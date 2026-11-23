
# LOSTRUCT ----------------------------------------------------------------
# make sure to do this on an HPC!!
# load libraries
require(pacman, devtools)
# remove comment to install lostruct from github
# devtools::install_github("petrelharp/local_pca/lostruct")
pacman::p_load(lostruct, 
               data.table, 
               tidyverse, 
               cluster, 
               zoo, 
               MetBrewer,
               ggokabeito,
               gtools)

# make sure in right working directory
# setwd("/u/home/h/hyangg/project-vlsork/10x_bcfs/")
setwd("/Users/heidiyang/inversions/lostruct_outlier_bcfs")

# set parameters
window_size <- 100 # todesco et al. used 100 SNP windows
k_kept <- 40 # Number of MDS kept
min_windows <- 4 # Minimum number of windows to call outlier regions
max_distance_between_outliers <- 10 # Max distance between outlier windows
n_permutations <- 1000 # Number of permutation for chromosomal clustering test
min_cor <- 0.8 # Correlation threshold for collapsing MDS
z_thresh <- 2 # threshold for defining MDS outlier windows

# get file name
args <- commandArgs(trailingOnly = TRUE) #if using it from terminal
bcf.file <- args[1]
# bcf.file <- "path/to/file"

# lostruct
sites <- vcf_positions(bcf.file)
win.fn.snp <- vcf_windower(bcf.file, size=window_size, type="snp", sites=sites) # partitioning VCF by window size
snp.pca <- eigen_windows(win.fn.snp, k=2, mc.cores=4) # Local PCA

# Distance matrix for MDS (NA removed)
pcdist <- pc_dist(snp.pca, mc.cores=4) 
na.wins <- is.na(snp.pca[,1])
pcdist <- pcdist[!na.wins, !na.wins]
nan.wins <- pcdist[,1]=="NaN"
pcdist <- pcdist[!nan.wins, !nan.wins]

# MDS of removed NA distance matrix
mds <- cmdscale(pcdist, eig=TRUE, k=k_kept)
mds.coords <- mds$points
colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))

# get start and end positions for each window
win.regions <- region(win.fn.snp)()
win.regions <- win.regions[!na.wins,][!nan.wins,]
win.regions |> mutate(mid=(start+end)/2) -> win.regions
for (k in 1:k_kept){
  name = paste("mds", str_pad(k, 2, pad = "0"), sep="")
  win.regions$tmp <- "NA"
  win.regions <- win.regions |> rename(!!name := tmp)
}
for (i in 1:k_kept){
  j = i + 4
  win.regions[,j] <- mds.coords[,i] # adding corresponding mds values
}
win.regions$n <- 1:nrow(win.regions)

# save file as RDS
saveRDS(win.regions, file= paste(bcf.file, "lostruct.win100.windows.rds"))


# LOSTRUCT OUTLIERS ----------------------------------------------------------------
        
# read RDS file names
RDS_files <- mixedsort(list.files(path = "/Users/heidiyang/inversions/lostruct_outlier_bcfs/",
                        pattern = ".rds")) 

# make empty tibble to fill with completed mds tibbles for each chr
chr_mds_outlier_tibble <- tibble(
  chr = 1:12,
  mds_tibble = vector("list", 12)
)

# collect outlier coordinate info 
outlier_coordinate_df <- data.frame(
  chr = character(),
  mds = numeric(),
  start_coord = numeric(),
  end_coord = numeric(),
  length = numeric(),
  stringsAsFactors = FALSE
)

# define necessary functions
# this one finds sequential runs of outlier windows within 10 windows of each other
find_sequential_runs <- function(df) {
  # Initialize list to store subsets of df
  ranges <- list()
  range_count <- 1
  start_index <- 1
  
  for (i in 2:nrow(df)) {
    if (df$n[i] - df$n[i-1] > 10) { # checks if the window numbers are greater than 10
      # Store completed range as a dataframe subset
      ranges[[range_count]] <- df[start_index:(i-1), ]
      range_count <- range_count + 1
      start_index <- i # this starts the index at the last non-consecutive number
    }
  }
  
  # Add the final range
  ranges[[range_count]] <- df[start_index:nrow(df), ]
  
  return(ranges)
}

for (chr in 1:length(RDS_files)) {
  # get the zscore outliers for all MDS axes for chromosome
  win.regions <- readRDS(RDS_files[chr]) # read file
  # make empty tibble for the chr
  mds_outlier_tibble <- tibble(
    mds = colnames(win.regions)[5:14],
    outlier_df = vector("list", 10),
    outlier_windows = vector("list", 10)
  )
  
  # chr colors for MDS plots
  okabe_ito_12 <- c(
    "#E69F00",  # Orange
    "#56B4E9",  # Sky Blue
    "#009E73",  # Bluish Green
    "#F0E442",  # Yellow
    "#0072B2",  # Blue
    "#D55E00",  # Vermillion
    "#CC79A7",  # Reddish Purple
    "#999999",  # Gray
    "#000000",  # Black
    "#88CCEE",  # Light Cyan
    "#AA4499",  # Purple
    "#DDCC77"   # Buff/Tan
  )

  # identify outliers and outlier regions for each MDS axis 
  for (i in 1:nrow(mds_outlier_tibble)) {
    mds_chosen <- mds_outlier_tibble[i, 1] |> pull()
    # identify outlier regions for MDS axis (zscore > z_thresh) # change Z-score threshold above
    sd_mds <- win.regions |> 
      dplyr::summarize(sd_mds = sd(.data[[mds_chosen]])) |> 
      pull(sd_mds)
    mean_mds <- win.regions |> 
      dplyr::summarize(mean_mds = mean(.data[[mds_chosen]])) |> 
      pull(mean_mds)
    zscore_mds <- win.regions |> 
      mutate(zscore = (.data[[mds_chosen]] - mean_mds) / sd_mds) 
    zscore_outlier_mds <- zscore_mds |> 
      filter(abs(zscore) > z_thresh) |> 
      arrange(start)
    
    # graph MDS coordinate and genomic position 
    # uncomment to make graphs
    # p <- ggplot(data = zscore_mds, 
    #        aes(x = mid / 1000000, 
    #            y = .data[[mds_chosen]],
    #            color = abs(zscore) > z_thresh)) +
    #   geom_point(alpha = 0.5) +
    #   scale_color_manual(values = c("gray", okabe_ito_12[chr]),
    #                      labels = c("≤ 2", "> 2"), # change zscore threshold here
    #                      name = "|Z-score|") +
    #   labs(x = "Genomic position (Mb)", 
    #        y = mds_chosen,
    #                  title = paste("Chr", chr, "MDS values of 100 SNP windows")) +
    #          theme_minimal() 
    # # save plot
    # ggsave(filename = paste("Chr_", chr, "_", mds_chosen, ".png", sep=""))
    
    # add outliers to chr tibble
    outliers <- zscore_outlier_mds |> 
      dplyr::select(chrom, start, end, mid, mds_chosen, n, zscore)
    mds_outlier_tibble$outlier_df[i] <- list(outliers)
    
    # find windows
    outlier_windows_list <- list(find_sequential_runs(outliers))
    mds_outlier_tibble$outlier_windows[i] <- outlier_windows_list
    
    # add outlier coordinates to larger df
    for (n in 1:length(outlier_windows_list[[1]])) { # for each window for this MDS
      window_df <- outlier_windows_list[[1]][[n]] # specify window
      if (nrow(window_df) > 1) { # if the window is more than 1 100-SNP window
        start_coord <- window_df$start[1]
        end_coord <- window_df$end[nrow(window_df)]
        length <- ((end_coord - start_coord) / 1000) # put in terms of Mb
        window_coords <- data.frame(
          chr = window_df$chr[1], 
          mds = mds_chosen, 
          start_coord = start_coord, 
          end_coord = end_coord,
          length = length)
        outlier_coordinate_df <- rbind(outlier_coordinate_df, window_coords)
      }
    }
    }
  # put the big tibble into the bigger chr tibble lol
  chr_mds_outlier_tibble$mds_tibble[[chr]] <- mds_outlier_tibble
}

# PCA ---------------------------------------------------

# load libraries
pacman::p_load(SNPRelate, gdsfmt, SeqArray, ggokabeito)

#setwd
setwd("/Users/heidiyang/inversions/lostruct_outlier_bcfs/") # need to make BCFs of each outlier region

# load and format data - do this only once since the file has already been made
bcf.fn <- "chr12_lostruct_outlier_pt2.bcf" # change for respective region
# format header
h <- seqVCF_Header(bcf.fn)
# there's some funky issues with the header format but this works
h$format$Number[h$format$ID=="QV"] <- "."
h$format$Number[h$format$ID=="TY"] <- "."
h$format$Number[h$format$ID=="CO"] <- "."
# save as GDS file
seqBCF2GDS(bcf.fn, "qlob.chr12_outlier_pos_pt2.gds")

# open file
genofile <- seqOpen("qlob.chr12_outlier_pos_pt2.gds", readonly=FALSE)
set.seed(1000)

# make PCA of individuals for the lostruct outlier regions
pca <- snpgdsPCA(genofile,
                 num.thread=2, 
                 autosome.only=FALSE)

# make pca df for plotting
pca_df <- data.frame(
  ind_ID = pca$sample.id,
  PC1 = pca$eigenvect[,1],    # the first eigenvector
  PC2 = pca$eigenvect[,2],    # the second eigenvector
  stringsAsFactors = FALSE)
# fix those yucky IDs
split_names <- sapply(pca_df$ind_ID, function (x) str_split(x, "[.]"))
pca_df$ind_ID <-t(as.data.frame(lapply(split_names, function (x) paste(x[2:3], collapse="."))))

# make PCA
ggplot(pca_df, aes(x=PC1, y=PC2)) + 
  geom_text(aes(x=PC1, y=PC2, label=ind_ID)) +
  theme_bw() + 
  ggtitle("Chr 12 Outlier Cluster PCA") + xlab("PC1") + ylab("PC2")

# if there are outliers (visual inspection) then identify them and re-run PCA
outlier_samples <- pca$sample.id[!pca$sample.id %in%
                                   c("Qlob.CD16.2.00F",
                                     "Qlob.SW.35.00F")]
pca <- snpgdsPCA(genofile, 
                 sample.id = outlier_samples, # exclude outliers
                 num.thread=2, 
                 autosome.only=FALSE)
# re-do code above
pca_df <- data.frame(
  ind_ID = pca$sample.id,
  PC1 = pca$eigenvect[,1],    # the first eigenvector
  PC2 = pca$eigenvect[,2],    # the second eigenvector
  stringsAsFactors = FALSE)

# fix IDs
split_names <- sapply(pca_df$ind_ID, function (x) str_split(x, "[.]"))
pca_df$ind_ID <-t(as.data.frame(lapply(split_names, function (x) paste(x[2:3], collapse="."))))

# make PCA
ggplot(pca_df, aes(x=PC1, y=PC2)) + 
  geom_text(aes(x=PC1, y=PC2, label=ind_ID)) +
  theme_bw() + 
  ggtitle("Chr 12 Outlier Cluster PCA") + xlab("PC1") + ylab("PC2")

# kmeans of PCA clusters
kmeans_pca <- as.data.frame(pca$eigenvect) |> 
  dplyr::select(!!1) # this selects the first PC

# function to get silhouette score for each k from 2-10
silhouette_score <- function(k) {
  km <- kmeans(kmeans_pca, centers = k, nstart=25) # performs kmeans clustering
  ss <- silhouette(km$cluster, dist(kmeans_pca)) # calculates silhouette score
  mean(ss[, 3]) # calculates average across iterations
}

# get avg silhouette score for each k - need this for silhouette score
k <- 2:10
avg_sil <- tibble(`k` = k,
                  `ss` = sapply(k, silhouette_score))

# plot the silhouette scores to find optimal k
ggplot(avg_sil,
       aes(k, ss)) +
  geom_line() +
  geom_point() +
  labs(
    x = '# clusters (k)',
    y = 'Average silhouette score',
    title = "Silhouette scores for k = 2 through k = 10"
  ) +
  theme_minimal()

# now that we know optimal k, we run kmeans again on the data with k = 2
# chr1: 8, chr2: 3, chr3: 3, chr4: 3, chr5: 2, chr6: 2, chr7: 3

km_final <- kmeans(kmeans_pca,centers = matrix(c(-0.1, 0.05, 0.2), nrow = 3))

# total within cluster sum of square
km_final$tot.withinss 

# cluster sizes
km_final$size 

pca_df <- pca_df |> 
  mutate(cluster = as.factor(km_final$cluster))

write.csv(pca_df, "chr12_outlier_cluster_df", row.names = FALSE)
#mutate(cluster = case_when(
#km_final$cluster == 1 ~ 2,
#km_final$cluster == 2 ~ 1,
#TRUE ~ km_final$cluster
#) |> as.factor())

# plot
ggplot(data = pca_df, 
       aes(x = PC1, 
           y = PC2, 
           color = cluster)) +
  geom_point(size = 2) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  #scale_color_okabe_ito() +
  labs(x = "PC1 (28.6%)", 
       y = "PC2 (5.0%)") +
  theme_minimal(base_family = "Avenir Next") +
  theme(text = element_text(size = 20))

win.regions <- win.regions |> 
  
  ggplot(data = win.regions, 
         aes(x = mid / 1000000, 
             y = mds01), 
         color = cluster) +
  geom_point(show.legend = FALSE) +
  scale_color_okabe_ito(alpha=0.5) +
  labs(x = "Chr12 position (Mb)", 
       y = "MDS Coordinate 1") +
  theme_minimal()


# HETEROZYGOSITY ----------------------------------------------------------

# heterozygosity (from PLINK or vcftools output from terminal) - opted to just use variants here
het_tbl <- read.table("chr12_het.txt", header=FALSE) # replace with correct heterozygosity file
het_tbl <- het_tbl |> 
  rename(ind = V1, 
         o_hom_count = V2, 
         total = V3, 
         het_prop = V4) |> 
  mutate(ind_ID = t(as.data.frame(lapply(split_names, function (x) paste(x[2:3], collapse=".")))),
         cluster = pca_df$cluster) |> 
  dplyr::select(-ind)

ggplot(data = het_tbl,
       aes(x = cluster,
           y = het_prop)) +
  geom_boxplot(aes(fill = cluster),
               show.legend = FALSE) + 
  #scale_fill_okabe_ito() +
  #scale_fill_manual(values = met.brewer(name = "Isfahan2", n = 3)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
  ylab("Proportion of heterozygous sites") +
  theme_minimal(base_family = "Avenir Next") +
  theme(text = element_text(size = 20))





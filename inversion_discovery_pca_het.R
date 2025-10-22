
# LOSTRUCT ----------------------------------------------------------------

# load libraries
require(pacman, devtools)
devtools::install_github("petrelharp/local_pca/lostruct")
pacman::p_load(lostruct, 
               data.table, 
               tidyverse, 
               cluster, 
               zoo, 
               MetBrewer,
               ggokabeito)

# make sure in right working directory
#setwd("/u/home/h/hyangg/project-vlsork/10x_bcfs/")
setwd("/Users/heidiyang/inversions/lostruct_outlier_bcfs")

# set parameters
window_size <- 100 # todesco et al. used 100 SNP windows
k_kept <- 40 # Number of MDS kept
min_windows <- 4 # Minimum number of windows to call outlier regions
max_distance_between_outliers <- 10 # Max distance between outlier windows
n_permutations <- 1000 # Number of permutation for chromosomal clustering test
min_cor <- 0.8 # Correlation threshold for collapsing MDS

# get file name
args <- commandArgs(trailingOnly = TRUE) #if using it from terminal
bcf.file <- args[1]

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
        
# open RDS file
win.regions <- readRDS("qlob_chr12.lostruct.win100.windows.rds")

## kmeans clustering to determine regions ##
# first let's get the fit2d tibble to just have the 2 MDS coordinates
kmeans_mds <- win.regions |> 
  dplyr::select(mds01, mds02)

# set le seed
#set.seed(123)

# function to get silhouette score for each k from 2-10
silhouette_score <- function(k) {
  km <- kmeans(kmeans_mds, centers = k, nstart=25) # performs kmeans clustering
  ss <- silhouette(km$cluster, dist(kmeans_mds)) # calculates silhouette score
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
km_final <- kmeans(kmeans_mds, 2)
# chr1: 8, chr2: 3, chr3: 3, chr4: 3, chr5: 2, chr6: 2, chr7: 3

# total within cluster sum of square
km_final$tot.withinss 

# cluster sizes
km_final$size 

# let's include the cluster number back to our initial dataset
win.regions$cluster <- as.factor(km_final$cluster)

# plot
ggplot(data = win.regions, 
       aes(x = mds01, 
           y = mds02, 
           color = cluster)) +
  geom_point(show.legend = FALSE) +
  scale_color_okabe_ito() +
  labs(x = "MDS Coordinate 1", 
       y = "MDS Coordinate 2",
       title = "MDS 2D plot of chromosome 7") +
  theme_minimal()


# this is my publication-ready figure
ggplot(data = win.regions, 
       aes(x = mid / 1000000, 
           y = mds06, 
           color = cluster)) +
  geom_point(show.legend = FALSE) +
  labs(x = "Genomic position (Mb)", 
       y = "MDS Coordinate 1",
       title = "Chr 12 MDS values of 100 SNP windows") +
  theme_minimal()

# now get the zscore outliers
mds_chosen <- "mds01"
sd_mds <- win.regions |> 
  dplyr::summarize(sd_mds = sd(.data[[mds_chosen]])) |> 
  pull(sd_mds)
mean_mds <- win.regions |> 
  dplyr::summarize(mean_mds = mean(.data[[mds_chosen]])) |> 
  pull(mean_mds)
zscore_mds <- win.regions |> # outliers have zscores > 1.5
  mutate(zscore = (.data[[mds_chosen]] - mean_mds) / sd_mds) |> 
  filter(abs(zscore) > 1.5) |> arrange(start)
neg_outliers <- zscore_mds |> 
  dplyr::select(chrom, start, end, mid, mds01, n, zscore, cluster) |> 
  filter(zscore < 0)
pos_outliers <- zscore_mds |> 
  dplyr::select(chrom, start, end, mid, mds01, n, zscore, cluster) |> 
  filter(zscore > 0)

# group by cluster
if (length(neg_outliers) > 0) {
  plot(neg_outliers$n, col=neg_outliers$cluster)
  plot(neg_outliers$mid, neg_outliers$zscore, col=neg_outliers$cluster)
  neg_outliers |> 
    group_by(cluster) |> 
    summarize(n = n()) 
}

clust_1 <- neg_outliers |> 
  select(cluster == 1) 

if (length(pos_outliers) > 0) {
  plot(pos_outliers$n, col=pos_outliers$cluster)
  plot(pos_outliers$mid, pos_outliers$zscore, col=pos_outliers$cluster)
  pos_outliers |> 
    group_by(cluster) |> 
    summarize(n=n())
}

# manually inspect the outlier DFs - take the first window position start 
# and end position of the last window

# chr 1: ends at 8185, then goes to 8272 - index # 628
# chr1:45286549-48090915

# chr 2:
# window number regions:7882-8167
# chr2:47787995-49014048

# chr 3:
# positive outliers
# window numbers: 4786-5027
# chr3:46959541-47980788
# other option: 5037-5095
# chr3:47991405-48225119
# altogether now: chr3:46959541-48225119
# negative outliers
# chr3:44995219-45302950

# chr4:56649791-59154308
chr4_outlier_length <- 59.154308 - 56.649791
# LD:56649000-59154000
# LD for graphing:50000000-65000000 

# chr5:52804869-53003411
# chr5:57131380-60969051

# chr6:29793595-30541449
# chr6:29.793595-30.248828
chr6_outlier_length <- 30.248828 - 29.793595
# for LD:29000000-31000000

# chr7:26414905-37730411
# chr7:38327482-39251068

# chr8:23888482-25636212

# chr9:33.232891-33.711931
chr9_outlier_length <- 33.711931 - 33.232891
# LD: 30000000-36000000


# chr10:35277492-37786975 # this is one, there's a big jump in basepairs

# chr11:30357747-31858058

# chr12:20929525-21008573

# chr12:18675632-20019679 - this is the one! (and it's actually negative on the MDS axis)
chr12_outlier_length <-20.019679 - 18.675632


# PCA ---------------------------------------------------

# load libraries
pacman::p_load(SNPRelate, gdsfmt, SeqArray, ggokabeito)

#setwd
setwd("/Users/heidiyang/inversions/lostruct_outlier_bcfs/")

# load and format data - do this only once since the file has already been made
bcf.fn <- "chr12_lostruct_outlier_pt2.bcf"
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

# fix those yucky IDs
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

# heterozygosity (from PLINK or vcftools output from terminal)
het_tbl <- read.table("chr12_het.txt", header=FALSE)
het_tbl <- het_tbl |> 
  rename(ind = V1, 
         o_hom_count = V2, 
         total = V3, 
         het_prop = V4) |> 
  filter(ind != "Qlob.CD16.2.00F" & ind != "Qlob.SW.35.00F") |> 
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


# ENVIRONMENTAL ASSOCIATION -----------------------------------------------

# load packages
pacman::p_load(lfmm, tidyverse, ggplot2, stringr, corrplot,
               Hmisc, PerformanceAnalytics, car, factoextra,
               geodata, terra, ggnewscale, ggpubr, prism, sf)
setwd("/Users/heidiyang/inversions/lostruct_outlier_bcfs")

# load data
clust_df <- read.csv("chr12_outlier_cluster_df", header = TRUE)
bioclim_df <- read.csv("bioclim_inversions_data.csv", header = TRUE)
comb_df <- left_join(clust_df, bioclim_df, by = join_by(ind_ID == sample_ID))

cors_df <- comb_df[, -c(seq(1:10))] |> 
  dplyr::select(-geometry)

# Create empty data frame to store results
anova_results <- data.frame(
  variable = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each variable in cors_df
for (i in 1:ncol(cors_df)) {
  var_name <- names(cors_df)[i]
  var_data <- cors_df[,i]
  
  # Run ANOVA
  anova_model <- aov(get(var_name) ~ cluster, data = comb_df)
  p_val <- summary(anova_model)[[1]][["Pr(>F)"]][1]
  
  # Store results
  anova_results <- rbind(anova_results, 
                         data.frame(variable = var_name, 
                                    p_value = p_val))
  

  # Create temporary data frame for plotting
  plot_data <- data.frame(
      variable_value = var_data,
      cluster = comb_df$cluster
    )
    
  # Create boxplot
  p <- ggplot(plot_data, aes(y = variable_value, x = as.factor(cluster))) +
    geom_boxplot() +
    labs(title = paste("Boxplot for", var_name),
          subtitle = paste("ANOVA p-value =", round(p_val, 4)),
          y = var_name,
          x = "Cluster") +
    theme_minimal()
    
    print(p)
}

anova_results <- anova_results |> 
  arrange(p_value)

# PCA of climate variables 
pca_result <- prcomp(cors_df, scale. = TRUE, center = TRUE)
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))
fviz_pca_var(pca_result, col.var = "contrib",
             gradient.cols = c("#00AFBB", "goldenrod1", "#FC4E07"))
pc_scores <- pca_result$x[, 1:4]

pca_data <- data.frame(
  PC1 = pc_scores[, 1],
  PC2 = pc_scores[, 2],
  PC3 = pc_scores[, 3],
  PC4 = pc_scores[, 4],
  cluster = comb_df$cluster
)

# ANOVA
oneway.test(PC2~cluster, data = pca_data) # chr12 for PC2

# graph
ggplot(pca_data, aes(y = PC2, 
                     x = as.factor(cluster), 
                     fill = as.factor(cluster))) +
  geom_boxplot(draw_quantiles = TRUE) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) + 
  guides(fill = "none") +
  #geom_boxplot(width=0.1, fill = "white") +
  stat_compare_means(comparisons = list(c("1", "2"), c("1", "3"), c("2", "3")),
                     method = "t.test",
                     label = "p.signif") +  # uses asterisks instead of p-values
  xlab("cluster") + ylab("Climate PC2 (19.5% of variance)") +
  theme_minimal(base_family = "Avenir Next") +
  theme(text = element_text(size = 20))

# look at loadings - check which variables
pc2_loadings <- pca_result$rotation[, 2]
print(pc2_loadings)

pc2_sorted <- sort(abs(pc2_loadings), decreasing = TRUE)
print(pc2_sorted)

# climate PCA
ggplot(pca_data,
       aes(x = PC1,
           y = PC2)) +
  geom_point(aes(color = as.factor(cluster))) +
  scale_color_okabe_ito() +
  theme_minimal()




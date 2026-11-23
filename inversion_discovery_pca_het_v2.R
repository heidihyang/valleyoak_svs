# load libraries
pacman::p_load(SNPRelate, gdsfmt, ggsignif, tidyverse, SeqArray, cluster, stringr, ggplot2)
#args <- commandArgs(trailingOnly = TRUE)

# make GDS file to use SNPRelate - only do this once
#h <- seqVCF_Header(args[1])
h <- seqVCF_Header("/Users/heidiyang/qlob_10x_snps.vcf.gz")
# # there's some funky issues with the header format but this works
h$format$Number[h$format$ID=="QV"] <- "."
h$format$Number[h$format$ID=="TY"] <- "."
h$format$Number[h$format$ID=="CO"] <- "."
# save as GDS file
seqBCF2GDS("/Users/heidiyang/qlob_10x_snps.vcf.gz")

# set working directory
setwd("/Users/heidiyang/Downloads/")

# open file
coord <- str_remove(args[1], ".bcf.gds") # if using HPC - include name of region as argument
filename <- args[1]
genofile <- seqOpen(filename, readonly=FALSE)

set.seed(1000)

# make PCA of individuals 
pca <- snpgdsPCA(genofile,
                 num.thread=2, 
                 autosome.only=FALSE)

# make pca df for plotting
pca_df <- data.frame(
  full_ID = pca$sample.id,
  PC1 = pca$eigenvect[,1],    # the first eigenvector
  PC2 = pca$eigenvect[,2],   # the second eigenvector
  stringsAsFactors = FALSE)

# fix ID formatting
split_names <- sapply(pca_df$full_ID, function (x) str_split(x, "[.]"))
pca_df$ind_ID <-t(as.data.frame(lapply(split_names, function (x) paste(x[2:3], collapse="."))))

# initial PCA
ggplot(pca_df, aes(x=PC1, y=PC2)) +
  geom_text(aes(x=PC1, y=PC2, label=ind_ID)) +
  theme_bw() +
  xlab(paste("PC1 (", round(pca$varprop[1],4)*100, "%)", sep="")) +
  ylab(paste("PC2 (", round(pca$varprop[2],4)*100, "%)", sep="")) +
  ggtitle(paste(coord, "PCA |", length(seqGetData(genofile, "variant.id")), "SNPs"))
# ggsave(paste(coord, "_initial_pca", ".png", sep=""))

# check if structure follows 3 clusters expected of inversions
# kmeans of PCA clusters
kmeans_pca <- as.data.frame(pca$eigenvect) |>
  dplyr::select(!!1) # this selects the first PC

silhouette_score <- function(k) {
  km <- kmeans(kmeans_pca, centers = k, nstart=25) # performs kmeans clustering
  ss <- silhouette(km$cluster, dist(kmeans_pca)) # calculates silhouette score
  mean(ss[, 3]) # calculates average across iterations
}

# get avg silhouette score for each k - need this for silhouette score
k <- 2:5 # specify range of # clusters to test
avg_sil <- tibble(`k` = k,
                  `ss` = sapply(k, silhouette_score))
max_k <- avg_sil |> 
  slice_max(ss) |> 
  pull(k)

if (max_k == 3) {
  # run kmeans
  km_final <- kmeans(kmeans_pca, centers = c(min(pca_df$PC1), mean(pca_df$PC1), max(pca_df$PC1)))
  
  # add cluster info to the pca_df
  pca_df <- pca_df |> 
    mutate(cluster = as.factor(km_final$cluster))
  
  # save the pca df
  write.table(pca_df, file=paste(coord,"_initial.tsv", sep=""), sep="\t")
  
  # axis limits with padding
  x_range <- range(pca_df$PC1, na.rm = TRUE)
  y_range <- range(pca_df$PC2, na.rm = TRUE)
  
  x_pad <- diff(x_range) * 0.1
  y_pad <- diff(y_range) * 0.1
  
  xlims <- c(x_range[1] - x_pad, x_range[2] + x_pad)
  ylims <- c(y_range[1] - y_pad, y_range[2] + y_pad)

  # plot PCA clusters
  ggplot(data = pca_df, 
         aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(size = 2, alpha = 1) +
    #geom_text(aes(label = ind_ID), show.legend = FALSE) +
    scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
    xlab(paste("PC1 (", round(pca$varprop[1],4)*100, "%)", sep="")) +
    ylab(paste("PC2 (", round(pca$varprop[2],4)*100, "%)", sep="")) +
    ggtitle(paste(coord, "PCA clusters")) +
    coord_cartesian(xlim = xlims, ylim = ylims) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_bw() +
    theme(panel.grid.minor = element_line(size = 0.5), 
          panel.grid.major = element_line(size = 1),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 14),
          legend.position = c(0.02, 0.02),
          legend.justification = c(0, 0),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA))
  #ggsave(paste(coord, "initial_pca_clusts.png", sep="_"))
  
  # heterozygosity (from PLINK or vcftools output from terminal)
  het_tbl <- read.table(paste(coord, ".bcf", "_het.txt", sep=""), header=FALSE)
  het_df <- het_tbl |> 
    dplyr::rename(ind = V1, 
           o_hom = V2, 
           e_hom = V3,
           total = V4, 
           f = V5,
           het_prop = V6) |> 
    mutate(ind_ID = t(as.data.frame(lapply(split_names, function (x) paste(x[2:3], collapse=".")))),
           cluster = pca_df$cluster) |> 
    dplyr::select(-ind)

  # plot heterozygosity per putative cluster
  ggplot(data = het_df,
         aes(x = cluster, y = het_prop)) +
    geom_boxplot(aes(fill = cluster),
                 show.legend = FALSE) +
    geom_signif(
      comparisons = list(c("1", "2"),
                         c("1", "3"),
                         c("2", "3")),
      test = "t.test",
      map_signif_level = TRUE
    ) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
    ylab("Proportion of heterozygous sites") +
    theme_bw() +
    theme(panel.grid.minor = element_line(size = 0.5), 
          panel.grid.major = element_line(size = 1),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 14),
          legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
          ) +
    ggtitle(paste(coord, "heterozygosity"))
  #ggsave(paste(coord, "_initial_het.png", sep=""))
}

# Mahalanobis distance outliers
# from ClassDiscovery package - modified to fit data
mahalanobisQC <- function(pca_object, pcs) {
  ss <- pca$eigenvect[, 1:pcs]
  maha <- sapply(1:nrow(ss), function(i) {
    v <- matrix(1/apply(ss[-i,], 2, var), ncol=1)
    as.vector(ss[i,]^2 %*% v)
  })
  pmaha <- 1-pchisq(maha, pcs) # prob of seeing test stat under X^2 dist with PC # degrees of freedom
  data.frame(statistic=maha, p.value=pmaha, ind_ID=pca_df$ind_ID)
}
for (i in c(2, 5, 10)) { # checking outlier individuals using different # PCs
  mah_samples <- mahalanobisQC(pca$eigenvect, i) |> 
    rownames_to_column("full_ID") |> 
    filter(p.value >= 0.01) |> pull(full_ID)
  outliers <- pca_df$full_ID[!(pca_df$full_ID %in% mah_samples)]
  
  # re-run without outliers
  pca <- snpgdsPCA(genofile, 
                   sample.id = mah_samples, # exclude outliers identified by mah outlier test
                   num.thread=2, 
                   autosome.only=FALSE)
  
  # re-do code above - only after manual visual inspection
  pca_df <- data.frame(
    full_ID = pca$sample.id,
    PC1 = pca$eigenvect[,1],    # the first eigenvector
    PC2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
  
  # fix IDs
  split_names <- sapply(pca_df$full_ID, function (x) str_split(x, "[.]"))
  pca_df$ind_ID <-t(as.data.frame(lapply(split_names, function (x) paste(x[2:3], collapse="."))))
  
  # plot PCA
  ggplot(pca_df, aes(x=PC1, y=PC2)) +
    geom_text(aes(x=PC1, y=PC2, label=ind_ID)) +
    theme_bw() +
    xlab(paste("PC1 (", round(pca$varprop[1],4)*100, "%)", sep="")) +
    ylab(paste("PC2 (", round(pca$varprop[2],4)*100, "%)", sep="")) +
    ggtitle(paste(coord, "PCA |", length(seqGetData(genofile, "variant.id")), "SNPs")) +
    theme_bw() +
    theme(panel.grid.minor = element_line(size = 0.5), 
          panel.grid.major = element_line(size = 1),
    ) 
  #ggsave(paste(coord, "pca_", i, ".png", sep=""))

  # this is a repeat of stuff - doing again after removing outlier individuals based on Mahalanobis distance
  # kmeans of PCA clusters
  kmeans_pca <- as.data.frame(pca$eigenvect) |>
    dplyr::select(!!1) # this selects the first PC
  
  silhouette_score <- function(k) {
    km <- kmeans(kmeans_pca, centers = k, nstart=25) # performs kmeans clustering
    ss <- silhouette(km$cluster, dist(kmeans_pca)) # calculates silhouette score
    mean(ss[, 3]) # calculates average across iterations
  }
  
  # get avg silhouette score for each k - need this for silhouette score
  k <- 2:5
  avg_sil <- tibble(`k` = k,
                    `ss` = sapply(k, silhouette_score))
  max_k <- avg_sil |> 
    slice_max(ss) |> 
    pull(k)
  
  if (max_k == 3) {
    # run kmeans
    km_final <- kmeans(kmeans_pca, centers = c(min(pca_df$PC1), mean(pca_df$PC1), max(pca_df$PC1)))
    
    # add cluster info to the pca_df
    pca_df <- pca_df |> 
      mutate(cluster = as.factor(km_final$cluster))
    
    ggplot(data = pca_df,
           aes(x = PC1,
               y = PC2,
               color = cluster)) +
      geom_point(size = 2) +
      scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
      xlab(paste("PC1 (", round(pca$varprop[1],4)*100, "%)", sep="")) +
      ylab(paste("PC2 (", round(pca$varprop[2],4)*100, "%)", sep="")) +
      theme_minimal() +
      ggtitle(paste(coord, "PCA clusters"))
    
    #Compute data-driven axis limits with padding
    x_range <- range(pca_df$PC1, na.rm = TRUE)
    y_range <- range(pca_df$PC2, na.rm = TRUE)
    
    x_pad <- diff(x_range) * 0.1
    y_pad <- diff(y_range) * 0.1
    
    xlims <- c(x_range[1] - x_pad, x_range[2] + x_pad)
    ylims <- c(y_range[1] - y_pad, y_range[2] + y_pad)
    
    p <- ggplot(data = pca_df, 
           aes(x = PC1, y = PC2, color = cluster)) +
      geom_point(size = 2, alpha = 1) +
      #geom_text(aes(label = ind_ID), show.legend = FALSE) +
      scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
      xlab(paste("PC1 (", round(pca$varprop[1],4)*100, "%)", sep="")) +
      ylab(paste("PC2 (", round(pca$varprop[2],4)*100, "%)", sep="")) +
      ggtitle(paste(coord, "PCA clusters")) +
      coord_cartesian(xlim = xlims, ylim = ylims) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
      theme_bw() +
      theme(panel.grid.minor = element_line(size = 0.5), 
            panel.grid.major = element_line(size = 1),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 14),
            legend.position = c(0.02, 0.02),
            legend.justification = c(0, 0),
            legend.background = element_rect(fill = alpha("white", 0.7), color = NA))
    print(p)
    #ggsave(paste(coord, "_pca_clusts_", i, ".png", sep=""))
    
    # save the pca df
    write.table(pca_df, file=paste(coord,"_", i,".tsv", sep=""), sep="\t")
    
    # heterozygosity (from PLINK or vcftools output from terminal)
    #het_tbl <- read.table(paste(args[1], "het.txt", sep="_"), header=FALSE)
    het_tbl <- read.table("chr5_32930681-33068794.bcf_het.txt", header=FALSE)
    het_df <- het_tbl |> 
      dplyr::rename(ind = V1, 
             o_hom = V2, 
             e_hom = V3,
             total = V4, 
             f = V5,
             het_prop = V6) |> 
      filter(ind %in% mah_samples) |> 
      mutate(ind_ID = t(as.data.frame(lapply(split_names, function (x) paste(x[2:3], collapse=".")))),
             cluster = pca_df$cluster) |> 
      dplyr::select(-ind)
    
    # Run pairwise t-tests
    t_1_2 <- t.test(het_df$het_prop[het_df$cluster == 2],
                    het_df$het_prop[het_df$cluster == 1])
    t_2_3 <- t.test(het_df$het_prop[het_df$cluster == 2],
                    het_df$het_prop[het_df$cluster == 3])
    
    if (mean(het_df$het_prop[het_df$cluster == 2]) > mean(het_df$het_prop[het_df$cluster == 1]) &&
        mean(het_df$het_prop[het_df$cluster == 2]) > mean(het_df$het_prop[het_df$cluster == 3]) &&
        t_1_2$p.value < 0.05 &&
        t_2_3$p.value < 0.05) {
      
      ggplot(data = het_df,
             aes(x = cluster, y = het_prop)) +
        geom_boxplot(aes(fill = cluster),
                     show.legend = FALSE) +
        geom_signif(
          comparisons = list(c("1", "2"),
                             c("1", "3"),
                             c("2", "3")),
          test = "t.test",
          map_signif_level = TRUE
        ) +
        scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
        ylab("Proportion of heterozygous sites") +
        theme_minimal() +
        ggtitle(paste(coord, "heterozygosity"))
      
      p <- ggplot(data = het_df,
             aes(x = cluster, y = het_prop)) +
        geom_boxplot(aes(fill = cluster),
                     show.legend = FALSE) +
        geom_signif(
          comparisons = list(c("1", "2"),
                             c("1", "3"),
                             c("2", "3")),
          test = "t.test",
          map_signif_level = TRUE
        ) +
        scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")) +
        ylab("Proportion of heterozygous sites") +
        theme_bw() +
        theme(panel.grid.minor = element_line(size = 0.5), 
              panel.grid.major = element_line(size = 1),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 12),
              plot.title = element_text(size = 14),
              legend.position = c(0.02, 0.98),
              legend.justification = c(, 1),
              legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
        ) +
        ggtitle(paste(coord, "heterozygosity"))
      print(p)
      #ggsave(paste(coord, "_het_", i, ".png", sep=""))
    }
  }
}

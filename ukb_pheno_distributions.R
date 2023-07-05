source("helper_functions.R")
library(ggridges)
library(tidyverse)

estimates_dir <- "" # Smoothed estimates after removing PCs
regressed_estimates_dir <- "" # Unsmoothed values after removing PCs
pca_path <- ""
gwas_path <- ""
img_dir <- ""
base_path <- ""
clustering_dir <- ""

clustering_file <- "hdbscan_labels_min25_EPS0.5_ukbb_pca_with_ids_UMAP_PC25_NC5_NN10_MD0.01_euclidean_20200214_065646.txt"

# Get the IDs in a data frame
pca_data <- import_pca_data(pca_path)
ukbb_ids <- pca_data[,c(1,2)]
gwas_data <- read_delim(gwas_path, delim = "\t")
base_data <- import_base_data(base_path)

clusters <- read_csv(file.path(clustering_dir, clustering_file), col_names = c("cluster"))#, col_types = cols(col_factor()))
clusters$cluster <- factor(clusters$cluster,
                           levels = min(clusters$cluster):max(clusters$cluster))

phenos <- c("SR4_FVC",
            "SR4_FEV1",
            "SR4_Standheight",
            "SR4_BMI",
            "SR4_Weight",
            "SR4_X30000_White_blood_cell_.leukocyte._count",
            "SR4_X30010_Red_blood_cell_.erythrocyte._count",
            "SR4_X30120_Lymphocyte_count",
            "SR4_X30130_Monocyte_count",
            "SR4_X30140_Neutrophill_count",
            "SR4_X30150_Eosinophill_count",
            "SR4_X30160_Basophill_count")

phenotype_captions <- create_phenotype_captions()

# This is the palette used in the main UKB figure (cluster labels with colours)
pal_ukb <- c(17,9,20,13,21,7,14,24,19,4,3,16,5,2,15,11,22,8,25,18,6,23,0,10,-1,12,1)
names(pal_ukb) <- c("#CB696D","#E17596","#E77325","#3E8E91","#3881AF","#586A9A","#4AA956",
                    "#C4902B","#58A057","#DCBD2E","#B53446","#B55E45","#AC6228","#FFEB2B","#F781BF",
                    "#806B87","#95519F","#CB6651","#6C856F","#FF8301","#F4EB31","#AF597D","#e4531a",
                    "#449C74","#864F70","#FFA60F","#FFC81D")
scales::show_col(names(pal_ukb))

# Specify the phenotype here and we can compare the pre- and post-PC adjusted values
pheno <- "SR4_X30140_Neutrophill_count"
pheno <- "SR4_FEV1"

cluster_df <- cbind.data.frame(ukbb_ids, base_data, clusters)
cluster_df <- left_join(cluster_df, gwas_data[c("FID","IID",phenos)], by = c("FID"="FID","IID"="IID"))

consensus_cluster_effects_df <- readRDS(paste(estimates_dir, "/", pheno, ".RDS",sep = ""))

# Unsmoothed values -- these are the values after regressing 40 PCs
# Results are stored in a lm() object, and IDs used are stored separately
unsmoothed_values <- readRDS(paste(regressed_estimates_dir, "/", pheno, "_full.RDS",sep = ""))
unsmoothed_ids <- readRDS(paste(regressed_estimates_dir, "/", pheno, "_full_ids.RDS",sep = ""))
unsmoothed_residuals <- cbind.data.frame(unsmoothed_ids, unsmoothed_values$residuals)
colnames(unsmoothed_residuals) <- c("FID","IID","residuals")

plotting_df <- left_join(cluster_df, consensus_cluster_effects_df, by = c("FID"="FID","IID"="IID"))
plotting_df <- left_join(cluster_df, unsmoothed_residuals, by = c("FID"="FID","IID"="IID"))

##### Ridge plots (multiple clusters) #####
# Note there are some parameters in ggplot to order and colour the densities correctly
# Clusters
cls <- c(17,18,19,22)

# image params
h <- 5
w <- 8

# Our data for plotting
temp <- subset(plotting_df, cluster %in% cls)
temp <- temp[c("FID","IID","cluster",pheno,"residuals")]
xlims <- c(-6,6) # x-axis range

# Calculate unadjusted means
plotting_means_prepc <- temp %>%
  filter(., !! sym(pheno) != -9) %>%
  group_by(cluster) %>%
  dplyr::summarise(mean = mean(!! sym(pheno)))

# Run some t-tests against cluster 17
for (cl in cls){
  if (cl !=17){
    temp1 <- subset(temp, cluster==17 & !is.na(residuals))
    temp2 <- subset(temp, cluster==cl & !is.na(residuals))
    print(paste(cl))
    print(t.test(temp1[[pheno]], temp2[[pheno]]))
  }
}

# Before PC adjustment
dp <- ggplot() +
  geom_density_ridges(data = temp %>%
                        filter(., !! sym(pheno) != -9),
                      aes(x = !! sym(pheno), y = as.factor(cluster), fill = as.factor(cluster)), alpha = 0.8) +
  scale_fill_manual("Populations",
                    values = names(sort(pal_ukb[pal_ukb %in% cls])),
                    labels = rev(rep(paste("Cluster",rev(cls)))),
                    guide = guide_legend(reverse = TRUE)) +
  geom_vline(xintercept = plotting_means_prepc$mean,
             linetype = "dashed",
             colour = names(sort(pal_ukb[pal_ukb %in% cls]))) +
  xlim(xlims) +
  xlab(paste("Residuals for", phenotype_captions[[pheno]])) +
  ylab("Cluster") +
  ggtitle("Distribution of non-PC-adjusted residuals", subtitle = phenotype_captions[[pheno]]) +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=14),
        plot.title = element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=12))

#ggsave(filename = paste(img_dir, "/",pheno,"_ridge_prepc.jpeg", sep=""), 
#       plot = dp, width = w, height = h,
#       device = "jpeg")

# After PC adjustment
# Calculate means
plotting_means_postpc <- temp %>%
  filter(., !! sym(pheno) != -9) %>%
  group_by(cluster) %>%
  dplyr::summarise(mean = mean(residuals))

# Run some t-tests against cluster 17
for (cl in cls){
  if (cl !=17){
    temp1 <- subset(temp, cluster==17 & !is.na(residuals))
    temp2 <- subset(temp, cluster==cl & !is.na(residuals))
    print(paste(cl))
    print(t.test(temp1$residuals, temp2$residuals))
  }
}

dp <- ggplot() +
  geom_density_ridges(data = temp %>%
                        filter(., !! sym(pheno) != -9),
                      aes(x = residuals, y = as.factor(cluster), fill = as.factor(cluster)), alpha = 0.8) +
  scale_fill_manual("Populations",
                    values = names(sort(pal_ukb[pal_ukb %in% cls])),
                    labels = rev(rep(paste("Cluster",rev(cls)))),
                    guide = guide_legend(reverse = TRUE)) +
  geom_vline(xintercept = plotting_means_postpc$mean,
             linetype = "dashed",
             colour = names(sort(pal_ukb[pal_ukb %in% cls]))) +
  xlim(xlims) +
  xlab(paste("Residuals for", phenotype_captions[[pheno]])) +
  ylab("Cluster") +
  ggtitle("Distribution of PC-adjusted residuals", subtitle = phenotype_captions[[pheno]]) +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=14),
        plot.title = element_text(size=20),
        legend.title = element_text(size=18),
        legend.text = element_text(size=12))

#ggsave(filename = paste(img_dir, "/",pheno,"_ridge_postpc.jpeg", sep=""), 
#       plot = dp, width = w, height = h,
#       device = "jpeg")

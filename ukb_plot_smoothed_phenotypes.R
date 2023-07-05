library(tidyverse)
source("helper_functions.R")

estimates_dir <- ""
pca_path <- ""
img_dir <- ""

# Get the IDs in a data frame
pca_data <- import_pca_data(pca_path)
ukbb_ids <- pca_data[,c(1,2)]

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

for (pheno in phenos){
  
  consensus_cluster_effects_df <- readRDS(paste(estimates_dir, "/", pheno, ".RDS",sep = ""))
  
  proj_path <- ""
  
  h <- 10
  w <- 8
  s <- 0.1
  a <- 0.6
  
  proj_array <- read.csv(proj_path, header=F, sep=" ", col.names = c("dim1","dim2"))
  proj_array <- cbind.data.frame(ukbb_ids, proj_array)
  
  plotting_data <- left_join(proj_array, consensus_cluster_effects_df)
  
  consensus_plot <- ggplot(plotting_data[!is.na(plotting_data[,"consensus_cluster_effect"]),]) +
    geom_point(size = s, alpha = a, aes(x = dim1, y = dim2, colour = consensus_cluster_effect), na.rm = TRUE) +
    scale_colour_gradient2(low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                           name = quote(sigma), limits = c(-0.5, 0.5)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "lightgrey",
                                          colour = "lightgrey"),
          axis.title=element_text(size=15),
          plot.title=element_text(size=20),
          legend.title=element_text(size=18),
          legend.text=element_text(size=12)) +
    ggtitle(paste("Phenotype distribution after smoothing\n(", phenotype_captions[[pheno]],")",sep="")) +
    xlab("UMAP1") + ylab("UMAP2")
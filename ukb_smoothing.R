# Script to run smoothing across clusters given phenotypes and PCs as arguments.

##### Import libraries #####
library(dplyr)
library(tidyverse)

##### IMPORT ARGS #####
n_pcs <- as.numeric(commandArgs()[which(commandArgs()=="--max_pcs") + 1])
phenos <- commandArgs()[which(commandArgs()=="--phenos") + 1]
pca_path <- commandArgs()[which(commandArgs()=="--pca_path") + 1]
gwas_path <- commandArgs()[which(commandArgs()=="--gwas_path") + 1]
clustering_dir <- commandArgs()[which(commandArgs()=="--clustering_dir") + 1]
residuals_dir <- commandArgs()[which(commandArgs()=="--residuals_dir") + 1]
estimate_dir <- commandArgs()[which(commandArgs()=="--estimate_dir") + 1]
training_path <- commandArgs()[which(commandArgs()=="--train") + 1]
testing_path <- commandArgs()[which(commandArgs()=="--test") + 1]

##### Script begins here #####
source("code/helper_functions.R") # Define some custom functions

# Create a vector of all of the UKB clusterings we've created
ukb_cluster_files <- list.files(clustering_dir, pattern = "ukbb")
ukb_cluster_paths <- paste(clustering_dir, ukb_cluster_files, sep = "/")

# Import PCA coordinates and phenotype data (originally used in a GWAS, hence the current name)
pca_data <- import_pca_data(pca_path)
gwas_data <- read_delim(gwas_path, delim = "\t")

# Get the IDs in a data frame
ukbb_ids <- pca_data[,c(1,2)]

# Import training and testing IDs
training_ids <- read.table(training_path, col.names = c("FID","IID"))
testing_ids <- read.table(testing_path, col.names = c("FID","IID"))

# Reduce the GWAS data frame to just the IDs and phenotype
# This is somewhat redundant but is imported from code where we looped through phenotypes
selected_gwas_variables <- c("FID","IID",phenos)
gwas_data <- gwas_data[,colnames(gwas_data) %in% selected_gwas_variables]

# Import the PCA data
pca_data <- select(pca_data, -used.in.pca.calculation)

# Create the data frame for our model. Attach the specified number of PCs and GWAS data.
# If n_pcs is "0", don't attach PCs and work with phenotype data directly
if (n_pcs > 0){
  pc_vars <- c(paste("PC", seq(1, n_pcs), sep = ""))
  model_df <- left_join(ukbb_ids, 
                        pca_data %>%
                          select(., FID, IID, pc_vars),
                        by = c("FID"="FID","IID"="IID")
  )
  model_df <- left_join(model_df, gwas_data, by = c("FID"="FID", "IID"="IID"))
} else if (n_pcs == 0){
  model_df <- left_join(ukbb_ids, gwas_data, by = c("FID"="FID", "IID"="IID"))
}

# Loop through all specified phenotypes (again from previous R script, for cluster job should just be one phenotype)
for (pheno in phenos){
  cat(paste("Beginning pheno", pheno, "file"))
  
  # Derive the valid IDs for our data
  # Drop NA, drop -9 (also NA), and reduce to ID variables
  valid_ids <- model_df %>%
    filter(., !(is.na(!! sym(pheno)))) %>% 
    filter(., (!! sym(pheno))!=-9) %>%
    select(., FID, IID)
  
  # Subset our model data to just the IDs from the previous step
  full_data <- subset(model_df, IID %in% valid_ids$IID)
  training_data <- subset(model_df, IID %in% intersect(valid_ids$IID, training_ids$IID))
  testing_data <- subset(model_df, IID %in% intersect(valid_ids$IID, testing_ids$IID))
  
  # Record the IDs used
  full_data_ids <- full_data[c("FID","IID")]
  
  # NEW: Split into testing and training here. Can still do full data since it's not very computationally intense.
  
  # If we specify some PC regression
  if (n_pcs > 0){
    # Define our model (Pheno ~ PC1 + PC2 + ... + PC[n_pcs])
    covars <- c(pc_vars)
    f <- generate_lm_formula(covars, pheno)
    
    # Get residuals from the full data set (because we want consensus estimates for each individual)
    # Use eval() since that prints out the model (shouldn't change results, but makes diagnostics easier)
    full_model <- eval(bquote(lm(.(f), data = full_data)))
    summary_full_model <- summary(full_model)
    print("Full model summary")
    print(summary_full_model)
    
    # Extract the residuals and the IDs
    pheno_residuals <- cbind.data.frame(full_data_ids, summary_full_model$residuals)
    
    # Save the residuals and summary data
    saveRDS(summary_full_model, paste(residuals_dir, "/", pheno, "_full_model_PC", as.character(n_pcs), ".RDS", sep = ""))
    saveRDS(pheno_residuals, paste(residuals_dir, "/", pheno, "_residuals_PC", as.character(n_pcs), ".RDS", sep = ""))
    
    ## Do the same, but now with training and testing data
    # Get the model
    # Actually, maybe we don't need this at all...
    training_model <- eval(bquote(lm(.(f), data = training_data)))
    summary_training_model <- summary(training_model)
    print("Training model summary")
    print(summary_training_model)
    
    pheno_residuals <- cbind.data.frame(training_data$FID, training_data$IID, summary_training_model$residuals)
    
    # Rename the columns
    colnames(pheno_residuals) <- c("FID","IID", pheno)
    
    saveRDS(summary_training_model, paste(residuals_dir, "/", pheno, "_training_model_PC", as.character(n_pcs), ".RDS", sep = ""))
    saveRDS(pheno_residuals, paste(residuals_dir, "/", pheno, "_training_residuals_PC", as.character(n_pcs), ".RDS", sep = ""))
    
    rm(summary_training_model)
  } else if (n_pcs == 0){
    # Get residuals relative to overall mean
    pheno_residuals <- full_data[c("FID","IID",pheno)]
    pheno_residuals[,pheno] <- pheno_residuals[,pheno] - mean(pheno_residuals[[pheno]])
    saveRDS(pheno_residuals, paste(residuals_dir, "/", pheno, "_residuals_PC", as.character(n_pcs), ".RDS", sep = ""))
  }
  
  # In case this data frame already exists, delete it (from a previous R script so probably unnecessary)
  if (exists("running_mean")){
    rm(running_mean)
  }
  
  # Initialize the data frame with IDs of observations with valid values
  running_mean <- full_data[,c("FID","IID")]
  running_mean$running_sum <- 0
  running_mean$running_count <- 0
  
  counter <- 1
  
  # Loop through the UKB clusterings
  for (cluster_path in ukb_cluster_paths){
    cluster_df <- create_clustering_df(cluster_path, ukbb_ids)
    cluster_name <- str_split(tail(str_split(cluster_path, "/")[[1]], 1), ".txt")[[1]][1]
    
    cat(paste("",as.character(counter))) # Print counter value
    
    # Join clusters to data
    cluster_model_joined_df <- left_join(pheno_residuals, cluster_df, by = c("FID"="FID", "IID"="IID"))
    full_data_joined_df <- left_join(full_data[c("FID","IID")], cluster_df, by = c("FID"="FID", "IID"="IID"))
    
    # Calculate the cluster means using training IDs
    cluster_means_training <- cluster_model_joined_df %>%
      filter(IID %in% training_ids$IID) %>%
      group_by(cluster) %>%
      summarize(mean_size = mean(!! sym(pheno))) %>%
      filter(cluster!=-1)
    
    # We assign means based on training data
    cluster_means <- cluster_means_training
    
    # Join the cluster means to the data frame
    #cluster_model_joined_df <- left_join(cluster_model_joined_df, cluster_means, by = c("cluster"="cluster"))
    cluster_model_joined_df <- left_join(full_data_joined_df, cluster_means, by = c("cluster"="cluster"))
    
    current_counts <- as.numeric(!is.na(cluster_model_joined_df$mean_size))
    
    running_mean$running_count <- running_mean$running_count + current_counts
    running_mean$running_sum <- rowSums(cbind(cluster_model_joined_df$mean_size, running_mean$running_sum), 
                                        na.rm = TRUE)
    
    counter <- counter + 1
  }
  
  running_mean$consensus_cluster_effect <- with(running_mean, running_sum/running_count)
  consensus_cluster_effects_df <- running_mean
  
  saveRDS(consensus_cluster_effects_df, paste(estimate_dir, "/", pheno, "_PC", as.character(n_pcs), ".RDS", sep = ""))
}

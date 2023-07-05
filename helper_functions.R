library(data.table)

#' plot_clusters
#'
#' Given a dataset with 2D coordinates, cluster labels, and plot parameters, generate the plot
#'
#' @param in_data Input data
#' @param cluster_file Not used right now
#' @param h # Plot height
#' @param w # Plot width
#' @param s # Point size
#' @param a # Point alpha [0-1]
#' @param param_fsize # Font size
#'
#' @return # Returns a ggplot object
#'
#' @examples
#' 
#' p <- plot_clusters(cluster_df, file.path(clustering_dir,clustering_file), 10, 6.49 , 0.01, 0.4)
plot_clusters <- function(in_data, cluster_file, h, w, s, a, param_fsize=10, label_alpha=1){
  # cluster membership for each point
  cluster_labels <- in_data$cluster
  
  # get average coordinates of each cluster for labeling
  cluster_coords <- in_data[,which(colnames(in_data) %in% c("cluster","dim1","dim2"))]
  
  coord_means_1 <- cluster_coords %>%
    group_by(cluster) %>%
    summarize(dim1 = mean(dim1))
  
  coord_means_2 <- cluster_coords %>%
    group_by(cluster) %>%
    summarize(dim2 = mean(dim2))
  
  coord_means <- left_join(coord_means_1, coord_means_2)
  
  # remove mean of unclustered (generally not informative)
  coord_means <- coord_means[coord_means$cluster!=-1,]
  coord_means$cluster_label <- coord_means$cluster
  
  
  # plot by cluster
  p_clusters <- ggplot(data = in_data, aes(x=dim1, y=dim2, colour=as.factor(cluster))) +
    geom_point(size=s, alpha=a) + theme(legend.position="none") +
    geom_point(data = in_data[in_data$cluster=="-1",], size=s, alpha=a, colour = "black",
               aes(x=dim1, y=dim2)) +
    if (label_alpha!=1){
      # If we specify a label alpha
      geom_label_repel(data = coord_means, aes(x=dim1, y=dim2, label=cluster_label, colour = as.factor(cluster)), size=5,
                       alpha=label_alpha)
    } else {
      geom_label_repel(data = coord_means, aes(x=dim1, y=dim2, label=cluster_label, colour = as.factor(cluster)), size=5)
    }
  
  # uncomment this to annotate plots with details of the parameters and results
  #annotate(geom = "text", label = paste("HDBSCAN clustering. Black = unclustered", cluster_params),
  #         x = -Inf, y = Inf, hjust = 0, vjust = 1, colour = "black", size=param_fsize)
  
  return(p_clusters)
}

get_cluster_var_freqs <- function(in_data, col_name, n_to_keep){
  # Take in base data (COB, ethnicity categories, clusters)
  # Output a table of the most frequency sub-categories
  
  # passing data frame column names as objects is tricky:
  # https://stackoverflow.com/questions/48062213/dplyr-using-column-names-as-function-arguments
  
  temp_cob_freqs <- in_data %>%
    #group_by(cluster, COB) %>% # analysis will be done on COB within cluster
    group_by(cluster, !! sym(col_name)) %>%
    summarize(n = n()) %>% # get counts of individuals by COB
    arrange(cluster, desc(n)) %>% # sort by descending
    mutate(p.freq = n/sum(n)) %>% # calculate percentage
    top_n(n_to_keep, wt = n)
  
  return(temp_cob_freqs)
}

get_var_freqs <- function(in_data, col_name1, col_name2, n_to_keep){
  # Take in base data (COB, ethnicity categories, clusters)
  # Column names must be strings, and are converted to column objects
  # n_to_keep: keep the N largest populations
  # Output a table of the most frequency sub-categories
  
  # passing data frame column names as objects is tricky:
  # https://stackoverflow.com/questions/48062213/dplyr-using-column-names-as-function-arguments
  
  temp_cob_freqs <- in_data %>%
    group_by(!! sym(col_name1), !! sym(col_name2)) %>%
    summarize(n = n()) %>% # get counts of individuals by COB
    arrange(!! sym(col_name1), desc(n)) %>% # sort by descending
    mutate(p.freq = n/sum(n)) %>% # calculate percentage
    top_n(n_to_keep, wt = n)
  
  return(temp_cob_freqs)
}

get_var_freqs_2 <- function(in_data, col_name1, col_name2, col_name3, n_to_keep){
  # Take in base data (COB, ethnicity categories, clusters)
  # Column names must be strings, and are converted to column objects
  # Output a table of the most frequency sub-categories
  # n_subpop is the count of col2 x col3 within col1
  # n_pop is the count of col2 within col1
  # e.g. inputting cluster, COB, eth gives:
  # (1) number of individuals per cluster (n_)
  # (2) number of individuals born within a country within a cluster (n_pop)
  # (3) number of individuals born within a country within a cluster by self-identified ethnicity (n_subpop)
  # proportions returned are (3)/(1) and (2)/(1)
  # table is sorted by (3) within cluster, so we get the largest sub-populations
  
  # counts of col_name2 x col_name3 within col_name1 (e.g. Eth x COB within cluster)
  temp1 <- in_data %>%
    group_by(!! sym(col_name1), !! sym(col_name2), !!sym(col_name3)) %>%
    summarize(n_var3 = n())
  
  # counts of col_name2 within col_name1 (e.g. Eth within cluster)
  temp2 <- in_data %>%
    group_by(!! sym(col_name1), !! sym(col_name2)) %>%
    summarize(n_var2 = n())
  
  # counts of col_name1 (e.g. cluster)
  temp3 <- in_data %>%
    group_by(!! sym(col_name1)) %>%
    summarize(n_var1 = n())
  
  temp4 <- left_join(temp1, temp2)
  temp5 <- left_join(temp4, temp3)
  
  temp5 <- temp5 %>%
    mutate(pct_var3_of_var1 = n_var3/n_var1, pct_var2_of_var1 = n_var2/n_var1) %>%
    arrange(!! sym(col_name1), desc(n_var3)) %>%
    top_n(n_to_keep, wt = n_var3)
  
  return(temp5)
}

get_cluster_2var_freqs <- function(in_data, col_name1, col_name2, n_to_keep){
  # Generate three tables and then combine them
  # Table 1: Count the number of appearances of a subpopulation within a cluster
  # Table 2: Count the number of appearances of a population within a cluster
  # Table 3: Count the number of observations in a cluster
  
  # Combine the three tables to generate proportions within the cluster and within the population
  # This gives us, e.g., the COB proportion by cluster as well as by population within a cluster
  
  temp1 <- in_data %>%
    group_by(cluster, !! sym(col_name1), !! sym(col_name2)) %>%
    summarize(n_subpop = n())
  
  temp2 <- in_data %>%
    group_by(cluster, !! sym(col_name1)) %>%
    summarize(n_pop = n())
  
  temp3 <- in_data %>%
    group_by(cluster) %>%
    summarize(n_cluster = n())
  
  temp4 <- left_join(temp1, temp2)
  temp5 <- left_join(temp4, temp3)
  
  temp5 <- temp5 %>%
    mutate(pct_cluster = n_subpop/n_cluster, pct_pop = n_subpop/n_cluster) %>%
    arrange(cluster, desc(n_subpop)) %>%
    top_n(n_to_keep, wt = n_subpop)
  
  return(temp5)
}

import_cluster_labels <- function(c_dir, c_file){
  temp_table <- read_csv(paste(c_dir, c_file, sep="/"), col_names = FALSE)
  colnames(temp_table) <- c("cluster")
  return(temp_table)
}

join_cluster_labels <- function(df, cluster_df){
  return(data.frame(cbind(df, cluster_df)))
}

extract_cluster_params <- function(fname){
  # Extract cluster parameters based on the filename
  params <- list()
  
  params[["filename"]] <- fname
  params[["mp"]] <- as.integer(str_split(str_split(fname, "_min")[[1]][2], "_")[[1]][1])
  params[["pcs"]] <- as.integer(str_split(str_split(fname, "_PC")[[1]][2], "_")[[1]][1])
  params[["nc"]] <- as.integer(str_split(str_split(fname, "_NC")[[1]][2], "_")[[1]][1])
  params[["md"]] <-  as.numeric(str_split(str_split(fname, "_MD")[[1]][2], "_")[[1]][1])
  params[["nn"]] <-  as.integer(str_split(str_split(fname, "_NN")[[1]][2], "_")[[1]][1])
  params[["eps"]] <- as.numeric(str_split(str_split(fname, "_EPS")[[1]][2], "_")[[1]][1])
  
  return(params)
  
}

extract_cluster_counts <- function(labels){
  # Count the number of clusters as well as the number of unclustered  
  cluster_counts <- list()
  
  cluster_counts[["clusters"]] <- length(unique(labels[labels$cluster!=-1,]$cluster))
  cluster_counts[["unclustered"]] <- sum(labels$cluster==-1)
  
  return(cluster_counts)
}

strip_filename <- function(fname){
  # Strip periods and file extension from filename
  temp_str <- str_split(fname, ".txt")[[1]][1]
  return(gsub(pattern = "\\.", replacement = "", x = temp_str))
}

import_eth_colours <- function(){
  eth_colours <- read_csv("eth_colours.csv", 
                          col_names = TRUE,
                          col_types = cols(
                            X1 = "i",
                            Ethnicity = "f",
                            Colour = "c"
                          ))
  eth_colours <- select(eth_colours, -X1)
}

import_gwas_data <- function(pheno_data_path, pheno){
  # Import our phenotype data used in the GWAS
  # Filter it to non-missing data
  gwas_data <- read_delim(pheno_data_path, delim = "\t")
  gwas_data <- select(gwas_data, FID, IID, Sex, Age, !! sym(pheno)) 
  #%>%
  #  subset(., !! sym(pheno)!=-9)
  
  gwas_data <- gwas_data[!gwas_data[,pheno]==-9,]
  
  return(gwas_data)
}

# "data/ukbb_aux_data.csv"
import_aux_data <- function(aux_data_path){
  pd <- read_csv(aux_data_path,
                 col_types = cols(
                   eth_txt = col_factor(),
                   COB = col_factor()
                 ))
  
  # drop some redundant clusters and format data
  pd <- select(pd, -X1, -eid_x, -eid_str, -eth1_str, -eth2_str, -eid_y, -genetic_grouping,
               -eid_str_geo, -northing_filled, -easting_filled, -FID_height, -IID_height, -eid_covar_x, 
               -EID_age, -eid_asthma, -eid_covar_y, -eid_blood, -northing_orig, -easting_orig, -eth_code,
               FEV1 = `3063_0_0`, -COB, -eth_txt)
  
  return(pd)
}

# "data/ukbb_cluster_base.csv"
import_base_data <- function(base_data_path){
  bd <- read_csv(base_data_path,
                 col_names = TRUE,
                 na = ".",
                 col_types = cols(
                   X1 = "i",
                   COB = "f",
                   ethnicity_full = "f",
                   ethnicity_parent_cat = "f"))[,2:4]
  
  return(bd)
}

import_pca_data <- function(pca_data_path){
  pd <- read_delim(pca_data_path, delim = " ", col_names=TRUE)
  
  return(pd)
}

import_aux_data <- function(aux_data_path){
  pd <- read_csv(aux_data_path,
                 col_types = cols(
                   eth_txt = col_factor(),
                   COB = col_factor()
                 ))
  
  # drop some redundant clusters and format data
  pd <- select(pd, -X1, -eid_x, -eid_str, -eth1_str, -eth2_str, -eid_y, -genetic_grouping,
               -eid_str_geo, -northing_filled, -easting_filled, -FID_height, -IID_height, -eid_covar_x, 
               -EID_age, -eid_asthma, -eid_covar_y, -eid_blood, -northing_orig, -easting_orig, -eth_code,
               FEV1 = `3063_0_0`, -COB, -eth_txt)
  
  return(pd)
}


create_clustering_df <- function(cluster_path, id_frame){
  temp_table <- fread(cluster_path, header = FALSE, 
                      colClasses = c("factor"))
  temp_table <- cbind.data.frame(id_frame, temp_table)
  colnames(temp_table) <- c("FID", "IID", "cluster")
  return(temp_table)
}

create_umap_df <- function(umap_data_path, id_frame){
  temp_table <- fread(umap_data_path, header = FALSE)
  temp_table <- cbind.data.frame(id_frame, temp_table)
  return(temp_table)
}

generate_lm_formula <- function(covars, phenotype){
  return(as.formula(paste(phenotype, " ~ ",
                          paste(covars, collapse = " + "),sep = "")
  ))
}

MSE <- function(v1, v2){
  return(mean((v1-v2)^2))
}

get_cluster_means <- function(in_data, dim1, dim2, cluster_var){
  # in_date: data frame
  # dim1: character name of column for dimension 1
  # dim2: character name of column for dimension 2
  # cluster_var: variable for which to calculate means
  cluster_coords <- in_data[,which(colnames(in_data) %in% c(cluster_var, dim1, dim2))]
  
  coord_means_1 <- cluster_coords %>%
    group_by(!! sym(cluster_var)) %>%
    summarize(x1 = mean(!! sym(dim1)))
  
  coord_means_2 <- cluster_coords %>%
    group_by(!! sym(cluster_var)) %>%
    summarize(x2 = mean(!! sym(dim2)))
  
  coord_means <- left_join(coord_means_1, coord_means_2, by = c(cluster_var))
  
  return(coord_means)
}

create_phenotype_captions <- function(){
  # Create a set of label pairs
  phenotype_captions <- c()
  phenotype_captions[["SR4_FVC"]] <- "FVC"
  phenotype_captions[["SR4_FEV1"]] <- "FEV1"
  phenotype_captions[["SR4_Standheight"]] <- "Standing height"
  phenotype_captions[["SR4_BMI"]] <- "BMI"
  phenotype_captions[["SR4_Weight"]] <- "Weight"
  phenotype_captions[["SR4_X30000_White_blood_cell_.leukocyte._count"]] <- "Leukocyte count"
  phenotype_captions[["SR4_X30010_Red_blood_cell_.erythrocyte._count"]] <- "Erythrocyte count"
  phenotype_captions[["SR4_X30120_Lymphocyte_count"]] <- "Lymphocyte count"
  phenotype_captions[["SR4_X30130_Monocyte_count"]] <- "Monocyte count"
  phenotype_captions[["SR4_X30140_Neutrophill_count"]] <- "Neutrophil count"
  phenotype_captions[["SR4_X30150_Eosinophill_count"]] <- "Eosinophil count"
  phenotype_captions[["SR4_X30160_Basophill_count"]] <- "Basophil count"
  
  return(phenotype_captions)
}

store_results_sieb <- function(pheno_, model_, corr_, mse_, sieb_, n_train_, n_test_){
  # Store the results of the phenotype model in a list to add to a data frame
  # Requires phenotype, model, correlation, and MSE
  temp_list <- list()
  
  temp_list[["phenotype"]] <- pheno_
  temp_list[["model"]] <- model_
  temp_list[["correlation"]] <- corr_
  temp_list[["mse"]] <- mse_
  temp_list[["sieb"]] <- sieb_
  temp_list[["n_train"]] <- n_train_
  temp_list[["n_test"]] <- n_test_
  
  return(temp_list)
}

make_word_cloud_generic <- function(df_, group_, words_, group_val_=NA, scale_=c(2,3), pal_="black"){
  # INPUTS
  # df_: data (must be data.table format to support column syntax, i.e. the part with "..")
  # group_: grouping variable (e.g. cluster), a column in df_
  # words_: the variable to show in the word cloud (e.g. country of birth), a column in df_
  # group_val_: specific group to do a wordcloud for (e.g. cluster==1)
  # scale_: wordcloud paramter specifying the range of the size of the words
  # pal_: colour palette (or colours to use)
  
  # OUTPUT
  # cloud_: a wordcloud object from the package wordcloud
  
  # If a specific group is provided, subset the data to that
  if (is.na(group_val_)){
    print("No subset provided, carrying on with entire data frame.")
  } else {
    df_ <- df_[c(df_[,..group_] == group_val_),]
  }
  
  word_counts <- data.frame(table(df_[,..words_,]))
  colnames(word_counts) <- c("label","freq")
  
  wordcloud(words = word_counts[,"label"], freq = word_counts[,"freq"], colors = pal_,
            random.order = F, rot.per = 0,scale=scale_)
}

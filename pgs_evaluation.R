library(ggplot2)
library(ggrepel)
library(stringr)

# Cluster description table
# Attach the modal region(s) of birth and UKB ethnic background(s) for my convenience
cluster_descriptions <- data.frame(Cluster=integer(),
                                   Description=character(),
                                   Description_2=character(),
                                   Modal_region=character(),
                                   Model_ethnicity=character(),
                                   stringsAsFactors = FALSE)

cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(0,"0","JPN","Japan","Other ethnic group"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(1,"1","ITA","Italy","White, Any other white background"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(2,"2","FIN","Finland","White, Any other white background"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(3,"3","ENG-WHI","England","White"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(4,"4","ENG-WHI","England","White"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(5,"5","MID","Middle East","Other ethnic group"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(6,"6","NAF","North Africa","Other ethnic group"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(7,"7","SAS","South Asia","Asian or Asian British"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(8,"8","SAS-IND","India, Africa","Asian or Asian British, Indian"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(9,"9","SAS-MIX","England","Mixed, White and Asian"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(10,"10","SOM","Somalia","Black or Black British, African"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(11,"11","ENG-MIX","England","Mixed, Any other mixed background"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(12,"12","AMR","Central/South America","Other ethnic group"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(13,"13","SEA","East/Southeast Asia", "Chinese"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(14,"14","ENG-EAS-MIX","England","Mixed"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(15,"15","NEP","Nepal","Asian or Asian British; Other ethnic group"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(16,"16","IBR","Spain and Portugal","White, Any other white background"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(17,"17","ENG-BRI","England","White, British"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(18,"18","SSAFR","Sub-Saharan Africa","Black or Black British, African"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(19,"19","WAFR-CAR","Caribbean, England, W. Africa","Black or Black British"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(20,"20","ENG-BRI","England","White, British"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(21,"21","ENG-BRI-OTH","England","White, British; Any other white background"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(22,"22","ENG-CAR-WAB","England","Mixed, White and Black; Any other mixed background"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(23,"23","HAFR","Greater Horn of Africa","Black or Black British"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(24,"24","ENG-BRI","England","White, British"), stringsAsFactors = F)
cluster_descriptions <- rbind.data.frame(cluster_descriptions, list(25,"25","ENG-AFR-CAR-MIX","England, South Africa, Mauritius, Caribbean","Other ethnic group; Mixed"), stringsAsFactors = F)
colnames(cluster_descriptions) <- c("Cluster","Description","Description_2","ROB","ETH")

# Import Shadi's PGS data from VIPRS
shadi_prs_dir <- ""
prs_tables <- list()

# Get the weighted mean FST plink values
mean_weighted_fst <- read.csv("",
                              header=F,
                              stringsAsFactors = F)

# Set up the variables to match the other file I'd created
mean_weighted_fst$POP1 <- paste("C",as.character(mean_weighted_fst$V1),sep="")
mean_weighted_fst$POP2 <- paste("C",as.character(mean_weighted_fst$V2),sep="")
mean_weighted_fst$MEAN_WEIGHTED_FST <- mean_weighted_fst$V3
mean_weighted_fst <- mean_weighted_fst[,c(4,5,6)]

mean_weighted_fst <- subset(mean_weighted_fst, POP1=="C17"|POP2=="C17")
mean_weighted_fst$Cluster <- ifelse(mean_weighted_fst$POP1=="C17",
                                    as.numeric(sub("C","",mean_weighted_fst$POP2, fixed=T)),
                                    as.numeric(sub("C","",mean_weighted_fst$POP1, fixed=T)))

# Rename
mean_fst <- mean_weighted_fst
colnames(mean_fst) <- c("POP1","POP2","MEAN_FST","Cluster")

for (f in list.files(shadi_prs_dir)){
  prs_tables[[str_split(f, ".csv")[[1]][1]]] <- read.csv(
    paste(shadi_prs_dir, f, sep="/"), header=T
  )
}

##### Plot R2 versus FST #####
# Can create some CIs using the variance equations from Momin et al (2023)
# This is from the function r2redux::r2_var
# We will adjust it since we are using a mean R2, so Var(R2_bar) = Sum(Var(R2))/n^2
# Note: May do some bootstrapping code here as the mean of the variance might not be straightforward for this distribution

get_ci <- function(R2, mv2){
  t100 <- 1/mv2 * (1-R2)^2
  lam <- R2/t100
  
  t100 <- (t100^2) * 2 * (mv2 + 2*lam)
  var1 <- t100
  
  dvr2 <- R2
  mv <- mv2
  
  uci<-qchisq(0.975,1,ncp=lam)
  uci<-(uci-lam-1)/(2*(mv+2*lam))^.5
  uci<-uci*var1^.5+dvr2
  
  lci<-qchisq(0.025,1,ncp=lam)
  lci<-(lci-lam-1)/(2*(mv+2*lam))^.5
  lci<-lci*var1^.5+dvr2
  
  return(list(lci, uci))
}

get_lam <- function(R2,mv2){
  t100 <- 1/mv2 * (1-R2)^2
  lam_ <- R2/t100
  
  return(lam_)
}

get_var <- function(R2, mv2){
  t100 <- 1/mv2 * (1-R2)^2
  lam <- R2/t100
  
  t100 <- (t100^2) * 2 * (mv2 + 2*lam)
  var1 <- t100
}

# Image parameters
out_dir <- ""
h_ <- 7.5 # height
w_ <- 12 # width
dpi_ <- 600 # dpi
ext <- "png"
text_angle <- 0 # label text angle
text_size <- 8 # label text size

axis_title_size <- 25
axis_text_size <- 15
plot_title_size <- 25
legend_title_size <- 18
legend_text_size <- 15

offset <- 0.007 # text offset for "Reference population" label

# Generate plots for each of the phenotype predictions
for (pheno in names(prs_tables)){
  temp <- subset(prs_tables[[pheno]], Model=="VIPRS") # Subset to VIPRS output
  
  # Get Variances of individual R^2 values
  # Need to sum these up and divide by squared number of folds
  temp$R2_var <- get_var(temp$R2, temp$Sample_size)
  temp$lam <- get_lam(temp$R2, temp$Sample_size)
  
  cis <- get_ci(temp$R2, temp$Sample_size)
  
  temp$lci <- pmax(rep(0,length(cis[[1]])), cis[[1]])
  temp$uci <- cis[[2]]
  
  #temp <- subset(temp, Sample_size > 50)
  prs_r2_means <- aggregate(R2 ~ Cluster, temp, FUN=mean)
  colnames(prs_r2_means) <- c("Cluster","R2_mean")
  
  # There is no closed-form of the Generalized Chi-Square so for now just use the median R^2 value's CI
  prs_r2_median <- aggregate(R2 ~ Cluster, temp, FUN=median)
  colnames(prs_r2_median) <- c("Cluster","R2_median")
  
  prs_r2_lci <- temp[temp$R2 %in% prs_r2_median$R2, c("Cluster","lci")]
  prs_r2_uci <- temp[temp$R2 %in% prs_r2_median$R2, c("Cluster","uci")]
  
  prs_r2_sdevs <- aggregate(R2 ~ Cluster, temp, FUN=sd)
  names(prs_r2_sdevs) <- c("Cluster","R2_sdev")
  
  compared <- merge(prs_r2_means, mean_fst, by = c("Cluster"))
  compared <- merge(compared, prs_r2_median, by = c("Cluster"))
  compared <- merge(compared, prs_r2_sdevs, by = c("Cluster"))
  compared <- merge(compared, prs_r2_lci, by = c("Cluster"))
  compared <- merge(compared, prs_r2_uci, by = c("Cluster"))
  compared <- merge(compared, cluster_descriptions, by = c("Cluster"))
  
  self_compare <- subset(prs_tables[[pheno]], Model=="VIPRS" & Cluster==17)
  
  # Plot with error bars
  p <- ggplot(compared, aes(y = R2_median, x = MEAN_FST, label = Description)) +
    geom_point(aes(y = R2_median, x = MEAN_FST)) +
    geom_smooth(aes(y=R2_median,x=MEAN_FST),method = "lm", se=F) +
    geom_errorbar(aes(ymin = lci, ymax = uci), alpha=0.2) +
    geom_text_repel(hjust = -0.1, vjust = 0, angle = text_angle, size = text_size) +
    geom_hline(yintercept = mean(self_compare$R2), linetype="dotted") + # Reference cluster
    xlim(0,0.13) + expand_limits(y = 0) +
    xlab(expression(plain("Mean "*F[ST])*plain(" from Cluster 17"))) +
    ylab(expression(plain("Median ")*R^2*plain(" of VIPRS folds"))) +
    ggtitle(paste("Phenotype:", pheno)) 
  
  print(p)
  
  ggsave(filename = paste(pheno,"_median_r2_vs_fst_error_bars",".",ext,sep=""),
         plot = p,
         dpi = dpi_,
         device = ext,
         path = out_dir, h = h_, w = w_)
  
  # Plot mean
  p <- ggplot(compared, aes(y = R2_mean, x = MEAN_FST, label = Description)) +
    geom_point(aes(y = R2_mean, x = MEAN_FST), size=2) +
    geom_text_repel(hjust = -0.1, vjust = 0, angle = text_angle, size = text_size, alpha=0.7) +
    geom_hline(yintercept = mean(self_compare$R2), linetype="dotted") + # Reference cluster
    annotate("text", x = 0.12, y = mean(self_compare$R2) + offset, label = expression(plain("Reference population ")*R^2), size=text_size) +
    xlim(0,0.13) + expand_limits(y = 0) +
    xlab(expression(plain("Mean "*F[ST])*plain(" from Cluster 17"))) +
    ylab(expression(plain("Mean ")*R^2*plain(" of VIPRS folds"))) +
    ggtitle(paste("Phenotype:", pheno)) +
    theme(
      axis.title=element_text(size=axis_title_size),
      axis.text=element_text(size=axis_text_size),
      plot.title=element_text(size=plot_title_size),
      legend.title=element_text(size=legend_title_size),
      legend.text=element_text(size=legend_text_size))
  
  print(p)
  
  ggsave(filename = paste(pheno,"_mean_r2_vs_fst",".",ext,sep=""),
         plot = p,
         dpi = dpi_,
         device = ext,
         path = out_dir, h = h_, w = w_)
  
  # Plot median
  p <- ggplot(compared, aes(y = R2_median, x = MEAN_FST, label = Description)) +
    geom_point(aes(y = R2_median, x = MEAN_FST), size=2) +
    geom_text_repel(hjust = -0.1, vjust = 0, angle = text_angle, size = text_size, alpha = 0.7) +
    geom_hline(yintercept = mean(self_compare$R2), linetype="dotted") + # Reference cluster
    annotate("text", x = 0.12, y = mean(self_compare$R2) + offset, label = expression(plain("Reference population ")*R^2), size=text_size) +
    xlim(0,0.13) + expand_limits(y = 0) +
    xlab(expression(plain("Mean "*F[ST])*plain(" from Cluster 17"))) +
    ylab(expression(plain("Median ")*R^2*plain(" of VIPRS folds"))) +
    ggtitle(paste("Phenotype:", pheno)) +
    theme(
      axis.title=element_text(size=axis_title_size),
      axis.text=element_text(size=axis_text_size),
      plot.title=element_text(size=plot_title_size),
      legend.title=element_text(size=legend_title_size),
      legend.text=element_text(size=legend_text_size))
  
  print(p)
  
  ggsave(filename = paste(pheno,"_median_r2_vs_fst",".",ext,sep=""),
         plot = p,
         dpi = dpi_,
         device = ext,
         path = out_dir, h = h_, w = w_)
  
  rm(compared, prs_r2_means, prs_r2_sdevs, prs_r2_median, prs_r2_lci, prs_r2_uci)
  
}

##### Check allele frequencies #####
# Import the effect sizes from VIPRS and look at MAF vs R2
top_snps <- read.table("data/allele_freqs/ldl_c_top_100_snps.keep")
top_snps_betas <- read.csv("data/cluster_analysis_alex/data_2/snp_estimates/ldl_top_100_effect_sizes.csv")
colnames(top_snps) <- c("SNP")

allele_freqs_strat <- read.table("data/allele_freqs/ldl_c_top100.frq.strat",header=T)
allele_freqs_7412 <- subset(allele_freqs_strat, SNP=="rs7412")

# See what the results are for each of the top SNPs
lm_results <- list()

ldl_prs <- prs_tables[["LDL"]]
temp <- aggregate(R2 ~ Cluster, ldl_prs, FUN=median)

for (snp in top_snps$SNP){
  temp2 <- subset(allele_freqs_strat, as.character(SNP)==as.character(snp))
  
  df_lm <- merge(temp2, temp, by.x = "CLST", by.y = "Cluster")
  df_lm <- merge(df_lm, cluster_descriptions, by.x = "CLST", by.y = c("Cluster"))
  
  lm_rs <- lm(R2 ~ MAF, data = df_lm)
  lm_results[[snp]] <- lm_rs
}

r2_vec <- c()
pval_vec <- c()

for (l in lm_results){
  r2_vec <- c(r2_vec, summary(l)$r.squared)
  pval_vec <- c(pval_vec, summary(l)$coefficients[2,4])
}

plot(r2_vec, main="R2 values of the top 100 SNPs by order")
plot(-log10(pval_vec), main="p-values of the top 100 SNPs by order")

# Two highest R2 values are for SNPs rs7412 and rs4420638
# Look at them more closely; set up some output parameters
out_dir <- ""

# keep these image parameters the same as above for now, change if needed
h_ <- 7.5 # height
w_ <- 12 # width
dpi_ <- 600 # dpi
ext <- "png"
text_angle <- 0 # label text angle
text_size <- 8 # label text size

axis_title_size <- 25
axis_text_size <- 15
plot_title_size <- 25
legend_title_size <- 18
legend_text_size <- 15

# rs4420638 has the highest R2
lm_rs4420638 <- lm_results[["rs4420638"]]
summary_rs4420638 <- summary(lm_rs4420638)
df_rs4420638 <- subset(allele_freqs_strat, SNP=="rs4420638" & CLST!=-1)

medR2 <- aggregate(R2 ~ Cluster, ldl_prs, FUN=median)

df_rs4420638 <- merge(df_rs4420638, medR2, by.x = "CLST", by.y = "Cluster")
df_rs4420638 <- merge(df_rs4420638, cluster_descriptions, by.x = "CLST", by.y = c("Cluster"))

p <- ggplot(df_rs4420638,aes(y=R2, x=MAF, label=Description)) +
  geom_point(size=2) +
  geom_smooth(method="lm", se = F, colour="black", size=0.5) +
  geom_text_repel(hjust = -0.1, vjust = 0, angle = text_angle, size = text_size) +
  ylab(expression(plain("R"^2))) +
  xlab("Minor allele frequency of rs4420638") +
  ggtitle(expression(plain("R"^2)*plain(" of VIPRS vs MAF of rs4420638"))) +
  theme(
    axis.title=element_text(size=axis_title_size),
    axis.text=element_text(size=axis_text_size),
    plot.title=element_text(size=plot_title_size),
    legend.title=element_text(size=legend_title_size),
    legend.text=element_text(size=legend_text_size))

print(p)

# rs7412 has the next-highest R2 but has by far the highest absolute effect size
lm_rs7412 <- lm_results[["rs7412"]]
summary_rs7412 <- summary(lm_rs7412)
df_rs7412 <- subset(allele_freqs_strat, SNP=="rs7412" & CLST!=-1)

medR2 <- aggregate(R2 ~ Cluster, ldl_prs, FUN=median)

df_rs7412 <- merge(df_rs7412, medR2, by.x = "CLST", by.y = "Cluster")
df_rs7412 <- merge(df_rs7412, cluster_descriptions, by.x = "CLST", by.y = c("Cluster"))

p <- ggplot(df_rs7412,aes(y=R2, x=MAF, label=Description)) +
  geom_point(size=2) +
  geom_smooth(method="lm", se = F, colour="black", size=0.5) +
  geom_text_repel(hjust = -0.1, vjust = 0, angle = text_angle, size = text_size) +
  ylab(expression(plain("R"^2))) +
  xlab("Minor allele frequency of rs7412") +
  ggtitle(expression(plain("R"^2)*plain(" of VIPRS vs MAF of rs7412"))) +
  theme(
    axis.title=element_text(size=axis_title_size),
    axis.text=element_text(size=axis_text_size),
    plot.title=element_text(size=plot_title_size),
    legend.title=element_text(size=legend_title_size),
    legend.text=element_text(size=legend_text_size))

print(p)
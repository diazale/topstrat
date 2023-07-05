require(RColorBrewer)
library(scales)
library(tidyverse)
library(Ternary) 

##### Import cluster IDs and labels ####
source("helper_functions.R")
ids <- read.csv(file.path("data","1KGP_ids.txt"), sep = " ", header=F, colClasses = c("character","character"))

colnames(ids) <- c("FID","IID")

cluster_file <- "hdbscan_labels_min25_EPS0.5_1000G_UMAP_PC16_NC5_NN50_MD0.01_euclidean_2019814225811"

clust_labels <- read.csv(file.path("data",cluster_file), header=F)
colnames(clust_labels) <- "cluster"
clust_labels$cluster <- factor(clust_labels$cluster)

cluster_df <- cbind.data.frame(ids, clust_labels)

# Rename clusters to match MS figures
cluster_df <- cluster_df %>% mutate(
  cluster_2 = case_when(
    cluster==6 ~ 0,
    cluster==8 ~ 1,
    cluster==18 ~ 2,
    cluster==4 ~ 3,
    cluster==17 ~ 4,
    cluster==12 ~ 5,
    cluster==9 ~ 6,
    cluster==0 ~ 7,
    cluster==5 ~ 8,
    cluster==19 ~ 9,
    cluster==3 ~ 10,
    cluster==16 ~ 11,
    cluster==13 ~ 12,
    cluster==1 ~ 13,
    cluster==20 ~ 14,
    cluster==11 ~ 15,
    cluster==7 ~ 16,
    cluster==15 ~ 17,
    cluster==2 ~ 18,
    cluster==14 ~ 19,
    cluster==10 ~ 20
  )
)

cluster_df$cluster <- cluster_df$cluster_2
cluster_df <- cluster_df[,c(1:3)]

##### Import admixture proportions for CSA 1KGP populations #####
admix_dir <- "data/admixture_1kgp"

admix_ids_path <- file.path(admix_dir, "csa_ids.txt") # IDs
admix_3.p_path <- file.path(admix_dir, "csa_1kgp.3.P") # By base pair
admix_3.q_path <- file.path(admix_dir, "csa_1kgp.3.Q") # These are the proportions for individuals

# Import data and set up columns
admix_ids <- read.delim(admix_ids_path, header = F, sep = " ", colClasses = c("character","character","factor"))
colnames(admix_ids) <- c("FID","IID","Pop")
admix_props <- read.delim(admix_3.q_path, header = F, sep = " ")
colnames(admix_props) <- c("Pop1","Pop2","Pop3")

# Combine into one data frame
admix_df <- cbind.data.frame(admix_ids, admix_props)

# Merge in cluster values
admix_df <- left_join(admix_df, cluster_df, by = c("FID"="FID","IID"="IID"))

# Split by population
data_points <- list()
for (pop in unique(admix_df$Pop)){
  print(pop)
  data_points[[pop]] <- subset(admix_df, Pop==pop)
}

# Match cluster labels and colours to Fig 1 in MS
colpal <- c("#47A265","#6A6089","#AA5E28","#54A552","#A73D52","#E8D430","#6E8372","#3F918B","#FFC51C","#FFBC2D","#E41A1C","#EF7718","#BE6255","#CB6651","#FF980A","#F781BF","#3880B2","#DA718A","#896191","#C9992C","#A7558A")

# These are the indices for colours (respectively) CLM, MXL, PEL, PUR
# 6, 16, 17, 19
# In the end we don't need this, but I did the work, so I'm keeping it here just in case....
clr <- c()
clr[["CLM"]] <- colpal[6]
clr[["MXL"]] <- colpal[16]
clr[["PEL"]] <- colpal[17]
clr[["PUR"]] <- colpal[19]

# table(admix_df$Pop, admix_df$cluster)
# There 8 classes of population-cluster when we do the crosstab above
# 4 pops, 6 clusters
# Clusters are 0, 5, 10, 15, 16, 18
# So colours from colpal vector are 1, 6, 16, 17, 19 (since we take the cluster position and +1 for 1-based indexing in R)

pch_pop <- c()
pch_pop[["PUR"]] <- 15
pch_pop[["PEL"]] <- 16
pch_pop[["MXL"]] <- 17
pch_pop[["CLM"]] <- 18
cx <- 1.7 # cex (font size)

s <- 8 # Save as a square with side length s

# mac-specific thing: fixing plot.new() error
# other OSes may need to use the x11 command
quartz(width=s, height=s)
#png("manuscripts/images/1kgp_ternary.png",
#    width = s, height = s, units = "in", res = 600)

TernaryPlot(alab="Amerindigenous ancestry (%)",blab="African ancestry (%)",clab="European ancestry (%)",
            lab.cex=2,
            axis.cex = 2,
            #grid.col = "#fcfcfc",
            grid.lines = 2 # get rid of grid lines
)
legend(-0.355, 1.05,#-0.2, # position
       cex=1.5, # font size
       pch=pch_pop, # shape vector
       legend=c("PUR","PEL","MXL","CLM"), # label vector
       bty="n", # box
       ncol = 4, # n columns
       y.intersp = 0.5, # vertical space b/w text
       x.intersp = 0.5, # horizontal space b/w text
       xpd = TRUE, # allow legend outside of plot
       text.width = 0.1 # legend column spacing
) 

th <- 1.5 # line thickness
th_2 <- 1 # additional thickness
a <- 0.35 # alpha
ba <- 0.6 # alpha of bolded points
xa <- 0.4 # alpha of X

# PUR
p <- "PUR"
cl <- 18
AddToTernary(points, subset(admix_df, Pop==p & cluster==cl)[4:6], col=alpha(colpal[cl+1],a), pch=pch_pop[[p]], cex = cx, lwd = th)
cl <- 0
AddToTernary(points, subset(admix_df, Pop==p & cluster==cl)[4:6], col=alpha("black",xa), pch=3, cex = cx, lwd = th+th_2)
AddToTernary(points, subset(admix_df, Pop==p & cluster==cl)[4:6], col=alpha(colpal[cl+1],ba), pch=pch_pop[[p]], cex = cx, lwd = th+th_2)

# PEL
p <- "PEL"
cl <- 16
AddToTernary(points, subset(admix_df, Pop==p & cluster==cl)[4:6], col=alpha(colpal[cl+1],a), pch=pch_pop[[p]], cex = cx, lwd = th)
cl <- 5
AddToTernary(points, subset(admix_df, Pop==p & cluster==cl)[4:6], col=alpha("black",xa), pch=3, cex = cx, lwd = th+th_2)
AddToTernary(points, subset(admix_df, Pop==p & cluster==cl)[4:6], col=alpha(colpal[cl+1],ba), pch=pch_pop[[p]], cex = cx, lwd = th+th_2)

# MXL
p <- "MXL"
cl <- 15
AddToTernary(points, subset(admix_df, Pop==p & cluster==cl)[4:6], col=alpha(colpal[cl+1],a), pch=pch_pop[[p]], cex = cx, lwd = th)
cl <- 10
AddToTernary(points, subset(admix_df, Pop==p & cluster==cl)[4:6], col=alpha("black",xa), pch=3, cex = cx, lwd = th+th_2)
AddToTernary(points, subset(admix_df, Pop==p & cluster==cl)[4:6], col=alpha(colpal[cl+1],ba), pch=pch_pop[[p]], cex = cx, lwd = th+th_2)


# CLM
p <- "CLM"
cl <- 5
AddToTernary(points, subset(admix_df, Pop==p & cluster==cl)[4:6], col=alpha(colpal[cl+1],a), pch=pch_pop[[p]], cex = cx, lwd = th)
cl <- 15
AddToTernary(points, subset(admix_df, Pop==p & cluster==cl)[4:6], col=alpha("black",xa), pch=3, cex = cx, lwd = th+th_2)
AddToTernary(points, subset(admix_df, Pop==p & cluster==cl)[4:6], col=alpha(colpal[cl+1],0.8), pch=pch_pop[[p]], cex = cx, lwd = th+th_2)

#dev.off()
graphics.off()
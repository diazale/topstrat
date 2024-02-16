#!/bin/bash

# Storage directories
# Specify paths to data, output directory, log directory
D_SET=data/pca_1000g_100
OUT_DIR=projections
LOG_DIR=umap_logs

# Specify whether the source file has headers
HEADER=F

# Specify how many ID columns
N_ID=0


python umap_script.py \
-dset $D_SET \
-pc 20 25 \
-nn 25 50 \
-md 0 0.0001 \
-nc 3 5 \
-n_id $N_ID \
-outdir $OUT_DIR \
-head $HEADER \
-log $LOG_DIR
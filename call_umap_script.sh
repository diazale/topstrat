#!/bin/bash

# Storage directories

# Array of UMAP parameters to loop over
PC_ARRAY=(20 25)
NN_ARRAY=(25 50)
MD_ARRAY=(0 0.0001)
NC_ARRAY=(3 5)

# Specify paths to data, output directory, log directory
D_SET=data/pca_1000g_100
OUT_DIR=projections
LOG_DIR=umap_logs

# Specify whether the source file has headers
HEADER=F

for pc in "${PC_ARRAY[@]}"
	do
	for nn in "${NN_ARRAY[@]}"
		do
		for md in "${MD_ARRAY[@]}"
			do
			for nc in "${NC_ARRAY[@]}"
				do
				python umap_script.py \
					-dset $D_SET \
					-pc $pc \
					-nn $nn \
					-md $md \
					-nc $nc \
					-outdir $OUT_DIR \
					-head $HEADER \
					-log $LOG_DIR
				done
			done
		done
	done

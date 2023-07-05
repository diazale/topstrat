#!/usr/bin/env bash

# Path of the clustering script
PGM_PATH=hdbscan_clustering.py

# File with a list of projections to use
PROJ_LIST=umap_list.txt

# HDBSCAN parameters
EPS_LIST=(0.3 0.5)
MP_LIST=(25 50)

# Directories to store results
OUT_DIR=hdbscan_clusters
LOG_DIR=hdbscan_logs
PROB_DIR=hdbscan_probs

# Flag for whether the input file has headers
HEAD=F

# Flag for whether to return cluster membership probabilities
PROBS=T

# Number of ID columns (not implemented, but keeping this here for now)
N_ID=0

# Get the list of UMAPs to run
declare -a UMAPS=($(cat ${PROJ_LIST} | tr '\n' ' '))

# Loop over each UMAP file and each HDBSCAN parametrization
for UMAP in ${UMAPS[@]}
do
    for EPS in ${EPS_LIST[@]}
    do
        for MP in ${MP_LIST[@]}
        do
            echo "Beginning file ${UMAP}"
            echo "Epsilon: ${EPS}, Minimum points: ${MP}"

            python ${PGM_PATH} \
                    -dset ${UMAP} \
                    -min_points ${MP} \
                    -eps ${EPS} \
                    -n ${N_ID} \
                    -head ${HEAD} \
                    -probs ${PROBS} \
                    -outdir ${OUT_DIR} \
                    -probdir ${PROB_DIR} \
                    -log ${LOG_DIR}
        done
    done
done
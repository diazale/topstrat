"""
A script used to cluster data using HDBSCAN.
Technically works with any input data, but we used it with UMAP.

Inputs:
-dset : Input file
-min_points : Minimum number of points in a cluster
-head : Boolean, whether the input file has headers
-eps : Epsilon value for HDBSCAN(e_hat)
-n_id : NOT IMPLEMENTED, UNUSED (Number of columns to skip if data has ID columns at the start)
-probs : Boolean, whether to return cluster membership probabilities

-outdir : Output directory
-probdir : Output directory for probabilities
-log : Log directory

At the moment this script only returns cluster labels and membership probabilities.
However, it stores the clusterer object and could be used to return over values like the MST.

"""

import argparse
import hdbscan
import numpy as np
import os
import sys
import time


def str2bool(v):
    """
    Converts string to boolean

    :param v: Input value
    :return: Returns
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


parser = argparse.ArgumentParser(description='Run HDBSCAN on specified datasets.')

parser.add_argument('-dset', metavar='dataset', type=str,
                    help='Input dataset')
parser.add_argument('-min_points', type=int,
                    help='Minimum number of points in a cluster')
parser.add_argument('-head', metavar='headers', type=str2bool,
                    help='Indicate whether the file has headers')
parser.add_argument('-eps', metavar='epsilon', type=float,
                    help='Epsilon value for unequal population sizes')
parser.add_argument('-n_id', metavar='number_ID_cols', type=int,
                    help='Number of ID columns in the file')
parser.add_argument('-probs', metavar='membership_probability', type=str2bool,
                    help='Select whether to return cluster membership probabilities')

# Directories
parser.add_argument('-outdir', type=str,
                    help='Output directory for cluster labels')
parser.add_argument('-probdir', type=str,
                    help='Output directory for cluster membership probabilities')
parser.add_argument('-log', type=str,
                    help='Log directory')

args = parser.parse_args()

# Import arguments
dset = args.dset
min_points = args.min_points
has_headers = args.head
eps = args.eps
n_id = args.n_id  # Note that this is unused so far
probs = args.probs

out_dir = args.outdir
prob_dir = args.probdir
log_dir = args.log

fname = dset.split('/')[-1]

proj = np.loadtxt(dset)

# time stamp to keep track of files after multiple runs
tstamp = time.strftime('%Y%m%d_%H%M%S',time.localtime(time.time()))

print('File:', fname)
print('Minimum points:', str(min_points))
print('Epsilon for selection:', str(eps))

# set up logging
log_name = 'log_clustering_MP' + str(min_points) + '_' + tstamp + '_' + \
           fname.split('.txt')[0] + '.txt'
log_file = os.path.join(log_dir, log_name)

orig_stdout = sys.stdout # print() statements
orig_stderr = sys.stderr # terminal statements
f = open(log_file, 'w')
sys.stdout = f
sys.stderr = f

print('File:', fname)

start_time = time.time()
clusterer = hdbscan.HDBSCAN(min_cluster_size=min_points,cluster_selection_epsilon=eps).fit(proj)
end_time = time.time()

print('Time to cluster:', str(end_time - start_time))
print('Number of clusters:', str(max(clusterer.labels_)))
print('Number unassigned:', str(np.sum(clusterer.labels_ == -1)))

fname_temp = 'hdbscan_labels_min' + str(min_points) + '_EPS' + \
str(eps) + '_' + fname.split('.txt')[0] + '.txt'

fname_probs_temp = 'hdbscan_probs_min' + str(min_points) + '_EPS' + \
str(eps) + '_' + fname.split('.txt')[0] + '.txt'

np.savetxt(os.path.join(out_dir, fname_temp), clusterer.labels_, fmt='%i')

# If the user specifies probabilities, create the file
# We'll use 6 digits for a total of 8 characters
if probs:
    np.savetxt(os.path.join(prob_dir, fname_probs_temp), np.round(clusterer.probabilities_,6), fmt='%1.6f')

# restore print statements to terminal
sys.stdout = orig_stdout
sys.stderr = orig_stderr
f.close()

print('Time to cluster:', str(end_time - start_time))
print('Number of clusters:', str(max(clusterer.labels_)))
print('Number unassigned:', str(np.sum(clusterer.labels_ == -1)))
print('Cluster label file:', str(fname_temp))
print('Cluster membership strength file:', str(fname_probs_temp))
print()
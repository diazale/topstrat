"""
A Python script to carry out UMAP on the top PCs of data.

Inputs:
-dset : Input file
-pc : Number of PCs to use (default 10)
-nn : Number of nearest neighbours for UMAP (default 15)
-md : Minimum distance in low-dimensional space for UMAP (default 0.1)
-nc : Number of components (i.e. dimensionality) to reduce to (default 2)
-met : Distance metric (default Euclidian)
-outdir : Output directory to store UMAP results
-head : Boolean flag whether the input file has headers
-n_id : Number of columns to skip, if you have ID columns (default 2)
-log : Log directory

Output:
Creates a UMAP file with the parameters used in the filename, separated by underscores:

[Input dataset name]_
UMAP_PC[Number of PCs used]_
NC[Number of components]_
NN[Number of neighbours in UMAP]_
MD[Minimum distance in UMAP]_
met[Distance metric used]_
[timestamp]

The file names can get quite long but they are very informative for downstream analysis.

The corresponding log file will have the same format

"""

import argparse
import numpy as np
import os
import sys
import time
import timeit

import umap


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


example_text = '''
Requires the following packages: argparse, numpy, logging, os, sys, time, timeit, umap.

File names are generated automatically based on the input file name and parameters.

Files are time-stamped to avoid being overwritten from new runs.

EXAMPLE USAGE

Dataset: my_pcs.txt
Run UMAP on the top 15 PCs and the top 25 PCs, using 10 neighbours, a minimum distance of 0.001, reducing to 3D, 
assuming the first row of the file my_pcs.txt contains headers and we have 0 ID columns.

python umap_script.py \\ 
-dset ~/Documents/my_umap_project/my_pcs.txt \\ 
-pc 15 25 \\ 
-nn 10 \\ 
-md 0.001 \\ 
-nc 3 \\ 
-n_id 0 \\
-head T \\ 
-outdir ~/Documents/my_umap_project/umap_projections \\ 
-log ~/Documents/my_umap_project/logs
'''

# Define parser
parser = argparse.ArgumentParser(description='Run UMAP on specified datasets.',
                                 epilog=example_text,
                                 formatter_class=argparse.RawTextHelpFormatter)

# Input variables
parser.add_argument('-dset', metavar='DSET', type=str,
                    help='Input dataset')
parser.add_argument('-pc', metavar='PC', type=int, nargs='+',
                    default=10,
                    help='List of number of top PCs to use (default 10)')
parser.add_argument('-nn', metavar='NN', type=int, nargs='+',
                    default=15,
                    help='List of number of neighbours for UMAP (default 15)')
parser.add_argument('-md', metavar='MD', type=float, nargs='+',
                    default=0.1,
                    help='List of minimum distances for UMAP (default 0.1)')
parser.add_argument('-nc', metavar='NC', type=int, nargs='+',
                    default=2,
                    help='List of low dimensional components to project to (default 2D)')
parser.add_argument('-met', metavar='MET', type=str,
                    default='euclidean',
                    help='Type of distance metric to use (default euclidean)')
parser.add_argument('-outdir', metavar='OUTDIR', type=str,
                    default='umap_data',
                    help='Output directory (default: umap_data)')
parser.add_argument('-head', metavar='HEADERS', type=str2bool,
                    help='Indicate whether the file has headers')
parser.add_argument('-n_id', metavar='NID', type=int,
                    default=2,
                    help='Number of ID columns (default 2)')
parser.add_argument('-log', metavar='LOG', type=str,
                    default='umap_logs',
                    help='Log directory (default: umap_logs)')

args = parser.parse_args()

# Timestamp to identify different runs of UMAP
tstamp = time.strftime('%Y%m%d_%H%M%S',time.localtime(time.time()))

# Import arguments
dset = args.dset
pcs = args.pc
nns = args.nn
mds = args.md
ncs = args.nc
met = args.met.lower()
out_dir = args.outdir
has_headers = args.head
n_id = args.n_id
log_dir = args.log

# Check if important parameters have been left empty
crit=False
if dset is None:
    print('ERROR: No input dataset specified.')
    crit=True

if has_headers is None:
    print('ERROR: Headers in file not specified.')
    crit=True

if crit:
    # Exit if there is a critical error.
    sys.exit(1)

# If the output and log directories don't exist, create them
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Import the full data set and subset to the specified PCs as required
try:
    print('Beginning import of data:', dset)
    with open(dset) as data:
        data_contents = data.readlines()

        pca_data = []

        # import top PCs
        if has_headers:
            # If it has headers, skip the first row
            for p in data_contents[1:]:
                pca_data.append(p.split()[n_id:len(p)])  # Skip the ID columns
        else:
            # If there are no headers, just strip the ID columns
            for p in data_contents:
                pca_data.append(p.split()[n_id:len(p)])  # Skip the ID columns

        pca_data_array = np.array(pca_data).astype(float)

        print(pca_data_array.shape)

        # Cleanup
        del pca_data
        del p
        del data_contents
except Exception as e:
    print(e)
    print('Error during data import')

    sys.exit(1)

for pc in pcs:
    for md in mds:
        for nn in nns:
            for nc in ncs:

                # Make sure the number of components is >= the number of PCs
                # UMAP will still work but this is not well characterized.
                if pc < nc:
                    print('WARNING: Number of PCs is less than request dimensions.')
                    #sys.exit(1)

                # Record parameters in a string
                param_str = dset.split('/')[-1].split('.txt')[0] + '_UMAP_PC' + str(pc) + '_NC' + str(nc) + '_NN' \
                            + str(nn) + '_MD' + str(md) + '_' + met

                # Set up the log file
                log_file = os.path.join(log_dir, 'log_umap_' + param_str + '_' + tstamp + '.txt')

                print('Parameters: ', '\n PCs:', str(pc), '\n NC:', str(nc), '\n NN:', str(nn), '\n MD:', str(md),
                      '\n Metric:', met, '\n Has headers:', str(has_headers))

                # set up logging
                orig_stdout = sys.stdout  # print() statements
                orig_stderr = sys.stderr  # terminal statements
                f = open(log_file, 'w')
                sys.stdout = f
                sys.stderr = f

                print('Parameters: ', '\n PCs:', str(pc), '\n NC:', str(nc), '\n NN:', str(nn), '\n MD:', str(md),
                      '\n Metric:', met, '\n Has headers:', str(has_headers))

                fname = param_str + '_' + tstamp + '.txt'

                # preamble for log
                print()
                print("Using UMAP version: " + umap.__version__)
                print("Reducing to " + str(nc) + " components")
                print("Using " + str(nn) + " neighbours")
                print("Using minimum distance of " + str(md))
                print("Using metric: " + met)
                print("Using " + str(pc) + " PCs")
                print()
                print("Input data shape: ", pca_data_array.shape)

                try:
                    # Carry out UMAP
                    start = timeit.default_timer()
                    umap_proj = umap.UMAP(n_components=nc, n_neighbors=nn,min_dist=md,metric=met,
                                          verbose=True).fit_transform(pca_data_array[:, :pc])
                    stop = timeit.default_timer()
                except Exception as e:
                    print(e)
                    print('Error during UMAP')

                    f.close()
                    print('Error during UMAP')
                    sys.exit(1)

                print()
                print("UMAP runtime: ", stop - start)

                out_file = os.path.join(out_dir, fname)

                print()
                print("Output file:", out_file)
                print("Output data shape:", umap_proj.shape)

                np.savetxt(out_file, umap_proj)

                # Cleanup
                del umap_proj

                # restore print statements to terminal
                sys.stdout = orig_stdout
                sys.stderr = orig_stderr
                f.close()

                # print runtime to terminal.
                print("Finished successfully! UMAP runtime: ", stop - start)

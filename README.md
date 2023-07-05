# Topological stratification of biobank data

This repo contains code and documentation related to the manuscript *Topological stratification of continuous genetic variation in large biobanks* by Diaz-Papkovich et al<sup>1</sup>. Our methodology uses UMAP<sup>2</sup> and the Malzer and Baum update to HDBSCAN<sup>3</sup>. 

## Dimensionality reduction and clustering code

We use the Python implementations of both methods. We have provided two Python scripts to carry out the dimensionality reduction and clustering:
1. `umap_script.py`
2. `hdbscan_clustering.py`

These scripts can be executed via command line.

The UMAP script uses 8 parameters and assumes the input is the PCA data:
* Input data
* Number of PCs<sup>^</sup> to use
* Number of neighbours
* Minimum distance in the low dimensional space
* Number of components (dimensions) to reduce to. For clustering, we recommend a value of at least 3.
* Flag to indicate if the input file has headers
* Output directory
* Log directory

To execute the script, you can run:

```
umap_script.py \
-dset [data path] \
-pc [number of PCs] \
-nn [number of neighbours] \
-md [minimum distance] \
-nc [number of components] \
-head [header flag] \
-outdir [output directory] \
-log [log directory]
```

The output filename contains the parameters, input file, and a timestamp:

`[input_file]_UMAP_PC[top_pcs]_NC[num_dimensions]_NN[num_neighbours]_MD[min_dist]_[distance]_[timestamp].txt`

The output format is a space-delimited text file generated by `numpy.savetxt`. It only outputs UMAP coordinates, not IDs.

For the HDBSCAN script, there are 8 parameters:
* Input data (read in using `numpy.loadtxt`)
* Minimum points in a cluster
* Epsilon value
* Flag to indicate if the input file has headers
* Flag to indicate whether to return cluster membership probabilities
* Output directory for cluster labels
* Membership probability directory
* Log directory

To execute the script, you can run:

```
python hdbscan_clustering.py \
-dset [UMAP data] \
-min_points [minimum points] \
-eps [epsilon] \
-head [header flag] \
-probs [probability flag] \
-outdir [output directory] \
-probdir [probability directory] \
-log [log directory]
```

The output filename contains the parameters and the input file:

`hdbscan_labels_min[min_points]_EPS[epsilon]_[input_file].txt`

The output format is a single column of numeric cluster labels. The label `-1` indicates a point was not place in a cluster.

Since UMAP and HDBSCAN are not very computationally intensive, we recommend running a grid search across parameters. We have included two shell scripts to do this:
* `call_umap_script.sh`
* `call_hdbscan_script.sh`

They are set up to run a demonstration on 1KGP data by default (see below) but can be modified easily.
The general idea:
1. Run UMAP on your data to reduce its dimensionality. For visualization, we find relatively higher values of minimum distance to be useful (0.3 to 0.5). For clustering, values near or equal to 0 work better.
2. Run HDBSCAN on your dimensionally-reduced UMAP data. For the epsilon parameter, we find 0.3 to 0.5 to be useful for biobank data.

There is some discussion on parameters in Diaz-Papkovich et al<sup>1</sup>. 

### Demonstration
The two driver scripts are set up to use the 1KGP data, which is freely available [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip). We have provided the top 100 PCs of the 1KGP data as part of this repo. 

Once you have run `call_umap_script.sh`, you will need to run the following command in terminal:

`ls projections/*.txt > umap_list.txt`

This will create a list of the UMAP files generated, which can be fed into the HDBSCAN driver script.


For convenience, we have included some of the data in this repo, which is listed below. We have also included some demo code that visualizes the data:
* `1kgp_interactive_demo.Rmd`: the R notebook
* `1kgp_interactive_demo.html`: the [output](https://diazale.github.io/topstrat/1kgp_interactive_demo.html) of the R notebook 

## Manuscript code
We have included the code used in our manuscript:
* `1kgp_figures.Rmd`: Visualizations of 1KGP data, except for the ternary plot
* `1kgp_ternary_plot.R`: Used to generate ternary plots
* `helper_functions.R`: A variety of background functions to import data, export data, visualize clusters, etc
* `pgs_evaluation.R`: Used to evaluate polygenic scores by clusters and look at potentially influential alleles
* `ukb_pheno_distributions.R`: Used to visualize phenotype distributions for selected populations in the UKB
* `ukb_plot_smoothed_phenotypes.R`: Used to plot the smoothed phenotypes from the UKB
* `ukb_smoothing.R`: The smoothing algorithm used to regularize data over the 604 parametrizations

## 1KGP data included

These can be found in the `data` directory:

* `1000G_UMAP_PC100_NC2_NN15_MD0.5_admixturenotebook`: A 2D UMAP of 1KGP data for visualization used in the manuscript
* `1KGP_ids.txt`: Columns of the IDs from the 1KGP
* `admixture_1kgp`: A directory of estimated admixture proportions for the CLM, MXL, PEL, and PUR populations of the 1KGP (K=3) used in the manuscript
* `affy_samples.20141118.panel`: Sampled population labels for the 1KGP
* `hdbscan_labels_min25_EPS0.5_1000G_UMAP_PC16_NC5_NN50_MD0.01_euclidean_2019814225811`: Clusters used in the manuscript
* `pca_1000g_100`: Top 100 PCs for the 1KGP 

## Citations

1. Diaz-Papkovich et al. "Topological stratification of continuous genetic variation in large biobanks." (2023)
2. McInnes, Leland, John Healy, and James Melville. "UMAP: Uniform manifold approximation and projection for dimension reduction." arXiv preprint arXiv:1802.03426 (2018).
3. Malzer, Claudia, and Marcus Baum. "A hybrid approach to hierarchical density-based cluster selection." 2020 IEEE international conference on multisensor fusion and integration for intelligent systems (MFI). IEEE, 2020.

## Manuscript code

This repo also contains code used for visualization and analyses in our manuscript.
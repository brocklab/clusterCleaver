# clusterCleaver
clusterCleaver is a computational pipeline for analyzing cluster data in single-cell RNA sequencing (scRNAseq). You can check out our preprint [here](https://www.biorxiv.org/content/10.1101/2024.05.28.596337v1). This is a scanpy-compatible package which leverages the Earth Mover's Distance (EMD) to find genes which can be used to distinguish between clusters of cells. 

clusterCleaver provides a number of advantages:

1. clusterCleaver is fast and compatible with scanpy, meaning you can iterate through multiple potential parameters.
2. clusterCleaver provides easy to interpret plots, giving you confidence in its predictions.
3. Most importantly, clusterCleaver is one of the few surface marker identification tools which has validated its predictions.

# Installation
You can currently install clusterCleaver using pip:

`pip install clusterCleaver`

# Processing assumptions
clusterCleaver is designed around data which is processed according to the current single-cell best [practices](https://www.sc-best-practices.org/preamble.html). 

clusterCleaver works best on data which:
1. Is not regressed.
2. Has gene expression data which is greater than or equal to zero.
3. Has gene expression data which does not have excessive values.

The last caveat is important as some normalization methods can yield results with excessively high gene expression values. The EMD has high sensitivity to outliers. For our purposes, this is a good thing. However, if gene expression values are too high, this will result in spurious results. 

# Usage
To get started with clusterCleaver, check out `tutorial.ipynb` in the home directory. 

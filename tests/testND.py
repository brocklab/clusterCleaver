# %%
import scanpy as sc
import pandas as pd
from scipy.stats import wasserstein_distance_nd, wasserstein_distance
import numpy as np
from tqdm import tqdm

from clusterCleaver import searchOptimal, visualization

# %%
adata = sc.read_h5ad('../data/adata231Process.h5ad')
surfaceGenes = pd.read_csv('../data/surfaceGenes.csv', skiprows=4)
surfaceGenes = surfaceGenes['Column2'].tolist()
# %%
topGenes = searchOptimal.searchExpressionDist1D(
    adata, surfaceGenes, modifier='remove0', label='leiden'
)
topGenes = topGenes['genes'].tolist()
# %%
import itertools

label = 'leiden'
clusterLabels = list(adata.obs[label].unique())

is0 = np.array(adata.obs[label] == clusterLabels[0]).astype(bool)
is1 = np.array(adata.obs[label] == clusterLabels[1]).astype(bool)
nTop = 200
geneCombos = list(itertools.combinations(topGenes[0:nTop], 2))
print(len(geneCombos))

for combo in tqdm(geneCombos):
    surfaceIdx = np.where(adata.var.index.isin(combo))[0]
    X = adata.X[:, surfaceIdx]
    X0 = X[is0, :]
    X1 = X[is1, :]
    break
# %%
def sliced_wasserstein(X, Y, num_proj = 1000):
    """
    Computes the average sliced wasserstein distance for two arrays

    Inputs:
    X, Y: Input arrays (must have same number o columns)
    num_proj: Number of samples to compute distances

    Outputs:
    Mean wasserstein distance

    Notes:
    This was originally taken from
    https://stats.stackexchange.com/questions/404775/calculate-earth-movers-distance-for-two-grayscale-images
    based on the python optimal transport (POT) package.
    """
    dim = X.shape[1]
    ests = []
    # sample uniformly from the unit sphere
    dir1 = np.random.randn(dim, num_proj)
    dir2 = np.divide(dir1, np.linalg.norm(dir1, axis = 0))

    X0_proj = np.matmul(X, dir2)
    X1_proj = np.matmul(Y, dir2)

    ests = []
    for i in range(num_proj):
        ests.append(wasserstein_distance(X0_proj[:, i], X1_proj[:, i]))
    return np.mean(ests)
# %%
sliced_wasserstein(X0, X1)
import numpy as np

from tqdm import tqdm
import pandas as pd
import itertools
import warnings
from scipy.sparse import issparse
from scipy.stats import wasserstein_distance, gaussian_kde


def searchExpressionDist1D(
    adata,
    surfaceGenes,
    modifier="remove0",
    metric="EMD",
    label="leiden",
    minCounts=100,
    scale=False,
    maxCombos=10000,
):
    """
    Computes a statistical distance metric on gene expression data.

    Inputs:
        - adata: Anndata object with .obs consisting of numeric leiden cluster column
        - surfaceGenes: List of surface genes to compare, must also be in adata.var.index
        - modifier: Will remove 0s for improved EMD score, else run on unmodified gene expression
        - label: Column name in .obs containing label identities
        - nGenes: Number of genes to search through
        - minCounts: Number of counts sufficient for gene to pass when using modifier "remove0"
        - maxCombos: Maximum number of combinations of genes to search for. Becomes relevant with larger numbers.
    Outputs:
        - dfScores: Modified surface genes dataframe with a new separation score
    """
    metricDict = {"EMD": wasserstein_distance, "bhat": bhattacharyyaHist}

    # Validate input data
    nGenes = 1
    assert modifier in ["remove0", "no0", None], 'Modifier must be "remove0" or None'
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    scGenes = np.array(adata.var.index)

    availableGenes = [gene for gene in surfaceGenes if gene in scGenes]
    assert (
        len(availableGenes) > 0
    ), "No surface genes match genes in adata.var.index, check both inputs"

    surfaceCombos = list(itertools.combinations(availableGenes, nGenes))

    if len(surfaceCombos) > maxCombos:
        warnings.warn(
            "The number of combos generated is beyond the set maximum number of combos. Was this intentional?"
        )
    print(f"Searching for {len(surfaceCombos)} combinations of {nGenes} gene(s)")

    comboScores = []
    expressedClusters = []
    adata.obs[label] = adata.obs[label].astype("string")

    clusterLabels = list(adata.obs[label].unique())
    assert (
        len(clusterLabels) == 2
    ), "Number of unique labels in adata.obs[label] must be 2"

    is0 = np.array(adata.obs[label] == clusterLabels[0]).astype(bool)
    is1 = np.array(adata.obs[label] == clusterLabels[1]).astype(bool)
    surfaceCombosWrite = []
    for combo in tqdm(surfaceCombos):
        surfaceIdx = np.where(adata.var.index.isin(combo))[0]
        X = adata.X[:, surfaceIdx]

        # Continue if no expression is found
        if sum(X) == 0:
            continue
        # Can scale data between 0 and 1 (not recommended)
        if scale:
            X = (X - min(X)) / (max(X) - min(X))

        X0 = X[is0, :]
        X1 = X[is1, :]

        # print(f"{combo} : {sum(X0)} \t {sum(X1)} \t {surfaceIdx}")
        if sum(X0) == 0 or sum(X1) == 0:
            continue
        if nGenes == 1:
            X0 = X0.ravel()
            X1 = X1.ravel()
            if np.mean(X0) > np.mean(X1):
                cluster = "0"
            else:
                cluster = "1"
        else:
            cluster = -1

        X0, X1 = modifyEMD(X0, X1, modifier, minCounts=minCounts)
        distFunc = metricDict[metric]
        if len(X0) < minCounts or len(X1) < minCounts:
            continue
        dist = distFunc(X0, X1)
        comboScores.append(dist)
        expressedClusters.append(cluster)
        surfaceCombosWrite.append(combo)

    dfScores = pd.DataFrame(
        {
            "genes": np.array(surfaceCombosWrite).ravel(),
            "scores": comboScores,
            "cluster": expressedClusters,
        }
    )

    dfScores["genes"] = dfScores["genes"].astype("string")
    return dfScores.sort_values(by="scores", ascending=False)


def modifyEMD(X0, X1, modifier="remove0", minCounts=100):
    """
    Selectively removes gene expression with counts of 0s for the higher expressing cluster.

    Inputs:
    - X0, X1: Numpy arrays of gene expression for two separate clusters
    - modifier: Selects how to remove genes
                Currently only available as "remove0"
    - minCounts: Prevents selection of genes with low cell numbers after removal

    Outputs:
    - X0New, X1New: Modified gene expression values
    """
    X0 = X0.copy()
    X1 = X1.copy()

    if modifier not in ["remove0", "no0"]:
        return X0, X1
    # Other checks

    if modifier == "remove0":
        if X1.mean() > X0.mean():
            X1New = X1[X1 > 0]
            X0New = X0
        elif X0.mean() > X1.mean():
            X0New = X0[X0 > 0]
            X1New = X1
        else:
            return X0, X1
        if X0New.shape[0] < minCounts or X1New.shape[0] < minCounts:
            return X0, X1
        else:
            return X0New, X1New
    elif modifier == "no0":
        X1New = X1[X1 > 0]
        X0New = X0[X0 > 0]

        # if len(X1New) < minCounts:
        #     X1New = X1
        # if len(X0New) < minCounts:
        #     X0New = X0
        return X0New, X1New


def calculateKDE(y, covarianceFactor=0.25):
    kernel = gaussian_kde(y)
    kernel.covariance_factor = lambda: covarianceFactor
    kernel._compute_covariance()
    return kernel


def bhattacharyya(p, q):
    """bhattacharyya score, or measure of overlap"""
    return np.sum(np.sqrt(p * q))


def calculateBhattacharyya(X0, X1, ptsEval=10000):
    """
    Calculates the Bhattacharyya score by:
    - Finding a kernel density estimate
    - Normalizing KDE
    - Computing score

    Inputs:
    - X0: First gene expression vector
    - X1: Second gene expression vector
    - ptsEval: Number of points to evaluate bhjattacharyya score
        (will decrease speed on increase of value)
    """
    kernel0 = calculateKDE(X0)
    kernel1 = calculateKDE(X1)

    expr = np.linspace(0, max(np.concatenate([X0, X1])), ptsEval)

    X0KDE = kernel0(expr)
    X1KDE = kernel1(expr)

    X0KDE /= sum(X0KDE)
    X1KDE /= sum(X1KDE)

    bScore = bhattacharyya(X0KDE, X1KDE)

    return bScore


def bhattacharyyaHist(p, q):
    """
    Very fast (vectorized) bhattacharyya coefficient calculator
    Inputs:
    - p: List of observed values
    - q: List of observed values
    Outputs:
    - bScore: Bhattacharyya coefficient
    """
    # Grab relevant information for later calculations
    full = np.concatenate([p, q])
    maxFull = np.max(full)
    minFull = np.min(full)
    # Calculate and normalize counts
    histRange = (minFull, maxFull)
    hist1, _ = np.histogram(p, bins="auto", range=histRange)
    hist2, _ = np.histogram(q, bins="auto", range=histRange)
    hist1 = hist1 / sum(hist1)
    hist2 = hist2 / sum(hist2)

    bScore = bhattacharyya(hist1, hist2)
    return bScore

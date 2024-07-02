# %%
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from pathlib import Path
import scipy.sparse

# %%
colors = ['#BB4E44', '#44B1BB', '#76BB44', '#8944BB']
fullPalette = list(colors + sns.color_palette('tab10'))
sns.set_palette(sns.color_palette(fullPalette))


# %%
def plotExpression(adata, genes, colorCol='leiden'):
    """
    Plots dual gene expression across two axes.

    Inputs:
    - adata: Anndata single-cell object
    - genes: list of two genes to plot
    - colorCol: How to color histograms by category

    Outputs:
    Plots dual gene expression across two axes.

    """
    assert len(genes) == 2, 'Must have two genes'
    X = adata[:, genes].X
    if scipy.sparse.issparse(X):
        X = X.toarray()

    dfExpr = pd.DataFrame([X[:, 0], X[:, 1], adata.obs[colorCol]]).T
    dfExpr.columns = [genes[0], genes[1], colorCol]
    sns.jointplot(data=dfExpr, x=genes[0], y=genes[1], hue=colorCol)


def plotHists(
    adata, gene, truncate0=False, colorCol='leiden', logScale=False, saveFig=''
):
    """
    Plots a histogram and stripplot of gene expression data.

    Inputs:
    - adata: Anndata single-cell object
    - gene: Gene of interest
    - truncate0: Removes 0 from histogram
    - colorCol: How to color histograms by category
    - logScale: Using logarithmic scale for plotting
    - saveFig: Save location of plot

    Outputs:
    Histogram and stripplot
    """
    assert (
        gene in adata.var.index
    ), 'Gene is not present in anndata object (adata.var.index)'
    surfaceIdx = np.where(adata.var.index.isin([gene]))[0][0]
    expression = adata.X[:, surfaceIdx]

    if scipy.sparse.issparse(expression):
        expression = expression.toarray()

    not0 = list(expression > 0)

    if truncate0:
        expression = expression[not0]
        colorVec = adata.obs[not0][colorCol].tolist()
    else:
        colorVec = adata.obs[colorCol]
    dfHist = pd.DataFrame(expression, colorVec).reset_index()
    dfHist.columns = [colorCol, 'expression']

    dfHist[colorCol] = dfHist[colorCol].astype('category')
    # , log_scale=(False, True)
    plt.figure(figsize=(7, 6))
    plt.subplot(211)
    sns.histplot(
        data=dfHist,
        x='expression',
        hue=colorCol,
        element='poly',
        # stat='proportion',
        log_scale=(False, logScale),
    ).set(xlabel='')

    plt.subplot(212)
    sns.stripplot(
        data=dfHist,
        x='expression',
        hue=colorCol,
        native_scale=True,
        legend=False,
        jitter = 0.45
        # jitter=True,
    ).set(xlabel=f'{gene} Expression')
    if len(saveFig) > 0:
        saveDirectory = Path(saveFig).parents[0]
        if saveDirectory.exists():
            plt.savefig(saveFig, dpi=500)
        else:
            print('Save path directory {saveDirectory} does not exist.')


def plotModifiedHists(x0, x1, gene='gene'):
    """
    Plot histograms of two vectors. Useful if custom modifications have been made

    Inputs:
    - x0, x1: Lists of expression values
    - gene: Gene name used for plotting purposes
    """
    colorCol = 'leiden'
    logScale = False
    label0 = np.repeat('0', len(x0))
    label1 = np.repeat('1', len(x1))
    colorVec = np.concatenate([label0, label1])
    dfHist = pd.DataFrame([colorVec, np.concatenate([x0, x1])]).T
    dfHist.columns = [colorCol, 'expression']

    dfHist[colorCol] = dfHist[colorCol].astype('category')
    # , log_scale=(False, True)
    plt.figure(figsize=(7, 6))
    plt.subplot(211)
    sns.histplot(
        data=dfHist,
        x='expression',
        hue=colorCol,
        element='poly',
        # stat='proportion',
        log_scale=(False, logScale),
    ).set(xlabel='')

    plt.subplot(212)
    sns.stripplot(
        data=dfHist,
        x='expression',
        hue=colorCol,
        native_scale=True,
        legend=False,
        alpha=0.5,
        # jitter = 0.45
        jitter=True,
    ).set(xlabel=f'{gene} Expression')

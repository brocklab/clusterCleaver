# %%
import scanpy as sc
import unittest

# %%
adata = sc.datasets.pbmc3k_processed()
adata.obs["bCells"] = 0
isBCell = adata.obs["louvain"].astype("string") == "B cells"
adata.obs.loc[isBCell, "bCells"] = 1

sc.pl.umap(adata, color="bCells")

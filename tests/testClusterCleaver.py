# %%
import scanpy as sc
import pandas as pd
import unittest

from clusterCleaver import searchOptimal, visualization

# %%
adata = sc.read_h5ad('../data/adata231Process.h5ad')
surfaceGenes = pd.read_csv('../data/surfaceGenes.csv', skiprows=4)
surfaceGenes = surfaceGenes['Column2'].tolist()
# %%
# topGenes = searchOptimal.searchExpressionDist1D(
#     adata, surfaceGenes, modifier='remove0', label='leiden'
# )


# %%
class testEMDSearch(unittest.TestCase):
    
    def testSearch1D(self):
        topGenes = searchOptimal.searchExpressionDist1D(
            adata, surfaceGenes, modifier='remove0', label='leiden'
        )

        esam = topGenes['genes'].tolist()[4]
        self.assertEqual(esam, 'ESAM', f'Gene (ESAM) not found {esam}')
if __name__ == '__main__':
    unittest.main()
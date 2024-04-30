import os
import sys
import umap
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import scvi
from scvi.model import TOTALVI
import matplotlib.pyplot as plt

sys.path.append('../utils_py')
import plotting_utils

out_dir = './out'
print(out_dir)

T_out_dir = os.path.join(out_dir, 'T_cells')
print(T_out_dir)

adata = ad.read_h5ad(os.path.join(T_out_dir, 'T_integrated.h5ad'))

sc.pl.umap(adata, color='cell_type_L2', frameon=False)

sc.pl.umap(
  adata,
  color=['CD3D', 'CD4', 'FOXP3', 'CTLA4', 'IL2RA', 'CD8A', 'TRDC', 'PASK', 'LRRN3', 'LEF1', 'TCF7', 'CCR7', 'IL7R', 'SELL', 'CD40LG', 'SLC4A10', 'KLRB1', 'GZMK', 'NKG7', 'EOMES', 'TBX21', 'PDCD1'],
	frameon=False,
	use_raw=True,
	legend_loc='none',
	colorbar_loc=None,
	vmax='p99',
	vmin='p1',
	ncols=5,
	wspace=0,
)

sc.pl.dotplot(
  adata,
  all_markers,
  groupby='cell_type_L2',
  dendrogram=False,
  cmap='Reds',
  swap_axes=False,
  use_raw=True,
  standard_scale='var',
)
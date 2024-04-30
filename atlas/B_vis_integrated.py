import os
import sys
import scvi
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

sys.path.append('../utils_py')
import plotting_utils

out_dir = './out'
print(out_dir)

B_out_dir = os.path.join(out_dir, 'B_cells')
print(B_out_dir)

adata = sc.read_h5ad(os.path.join(B_out_dir, 'batch_corrected_SCVI_B.h5ad'))

adata.uns['age_group_colors'] = np.array(['#00A087FF', '#00A087FF'])
datasets_palette = {
    'GSE157007': '#E64B35',
    'GSE213516': '#4DBBD5',
    'GSE214546': '#00A087',
    'HRA000203': '#3C5488',
		'HRA000624': '#F39B7F',
		'HRA003766': '#8491B4',
		'syn22255433': '#91D1C2'
}

sc.pl.umap(adata, color='dataset', groups=['GSE157007'], palette=datasets_palette, frameon=False)
sc.pl.umap(adata, color='dataset', groups=['GSE213516'], palette=datasets_palette, frameon=False)
sc.pl.umap(adata, color='dataset', groups=['GSE214546'], palette=datasets_palette, frameon=False)
sc.pl.umap(adata, color='dataset', groups=['HRA000203'], palette=datasets_palette, frameon=False)
sc.pl.umap(adata, color='dataset', groups=['HRA000624'], palette=datasets_palette, frameon=False)
sc.pl.umap(adata, color='dataset', groups=['HRA003766'], palette=datasets_palette, frameon=False)
sc.pl.umap(adata, color='dataset', groups=['syn22255433'], palette=datasets_palette, frameon=False)

sc.pl.umap(adata, color='age_group', groups=['young'], frameon=False)
sc.pl.umap(adata, color='age_group', groups=['old'], frameon=False)

sc.pl.umap(
  adata,
  color=['total_counts', 'pct_counts_mt'],
  frameon=False,
)

all_markers = ['CD79A', 'CD19', 'IGHM', 'IGHD', 'IGHA1', 'IGHA2', 'TBX21', 'ITGAX']

sc.pl.umap(
  adata,
  color='leiden_scVI_0.2',
  frameon=False,
  legend_loc='on data',
  size=3,
  palette='tab20c'
)

sc.pl.dotplot(
  adata,
  all_markers,
  groupby='leiden_scVI_0.2',
  dendrogram=False,
  cmap='Reds',
  swap_axes=False,
  use_raw=True,
  standard_scale='var',
)

adata.obs['cell_type_L2'] = 'other'
adata.obs['cell_type_L2'] = np.where(adata.obs['leiden_scVI_0.2'] == '1', 'B memory', adata.obs['cell_type_L2'])
adata.obs['cell_type_L2'] = np.where(adata.obs['leiden_scVI_0.2'] == '2', 'B atypical', adata.obs['cell_type_L2'])
adata.obs['cell_type_L2'] = np.where(adata.obs['leiden_scVI_0.2'] == '0', 'B naive', adata.obs['cell_type_L2'])

sc.pl.umap(
  adata,
  color=['cell_type_L2'],
  frameon=False,
  size=4,
  palette={
	  'B memory': '#A6DBBB',
    'B naive': '#80C39E',
		'B atypical': '#5AAC82',

		'other': '#ECD9F1'
  },
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

adata.obs.to_csv(os.path.join(T_out_dir, 'annotated_B_metadata.tsv'), sep='\t')
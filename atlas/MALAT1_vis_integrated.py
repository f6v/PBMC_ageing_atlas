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

sys.path.append('../utils_py')
import plotting_utils

out_dir = './out'
print(out_dir)

MALAT1_out_dir = os.path.join(out_dir, 'MALAT1_cells')
print(MALAT1_out_dir)

adata = sc.read_h5ad(os.path.join(MALAT1_out_dir, 'batch_corrected_SCVI_MALAT1.h5ad'))

adata.uns['age_group_colors'] = np.array(['#00A087FF', '#00A087FF'])

def highlight_dataset(adata, dataset_name):
	datasets_palette = {
    'GSE157007': '#E64B35',
    'GSE213516': '#4DBBD5',
    'GSE214546': '#00A087',
    'HRA000203': '#3C5488',
		'HRA000624': '#F39B7F',
		'HRA003766': '#8491B4',
		'syn22255433': '#91D1C2'
	}
	sc.pl.umap(
		adata,
		color='dataset',
		groups=dataset_name,
		palette=datasets_palette,
		frameon=False,
		size=10,
		legend_loc=None,
		save=f"{dataset_name}.png",
		title=''
	)

highlight_dataset(adata, 'GSE157007')
highlight_dataset(adata, 'GSE213516')
highlight_dataset(adata, 'GSE214546')
highlight_dataset(adata, 'HRA000203')
highlight_dataset(adata, 'HRA000624')
highlight_dataset(adata, 'HRA003766')
highlight_dataset(adata, 'syn22255433')

sc.pl.umap(adata, color='age_group', groups=['young'], frameon=False)
sc.pl.umap(adata, color='age_group', groups=['old'], frameon=False)

sc.pl.umap(
	adata,
  color=['total_counts', 'pct_counts_mt'],
  frameon=False,
	size=10
)

sc.pl.umap(
	adata,
  color=[
    'MALAT1', 'CD3D', 'CD4', 'CD8A',
		'IGHM', 'IGHA1', 'CD79A',
		'NKG7', 'KLRB1', 'NCAM1',
		'CCR7', 'PASK', 'LRRN3', 'EOMES', 'HAVCR2', 'PDCD1'	
  ],
  frameon=False,
  use_raw=True,
  legend_loc='none',
  colorbar_loc=None,
  vmax='p99',
  vmin='p1',
  ncols=4,
  wspace=0,
)

sc.tl.embedding_density(adata, basis='umap', groupby='age_group')
sc.pl.embedding_density(adata, basis='umap', key='umap_density_age_group')


import os
import sys
import scvi
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt

out_dir = './out'
print(out_dir)

adata = sc.read_h5ad(os.path.join(out_dir, 'batch_corrected_SCVI_cleaned.h5ad'))

celltypist_predictions = pd.read_csv(
  os.path.join(out_dir, 'celltypist_Immune_All_High.tsv'),
  sep='\t',
  index_col=0
)
adata.obs['celltypist_prediction'] = celltypist_predictions

def plot_celltypist(adata, ax, cell_type):
  sc.pl.umap(
    adata,
    ax=ax,
    color='celltypist_prediction',
    groups=[cell_type],
    legend_loc=None,
    frameon=False,
    show=False,
    palette=['#00A087FF']*24,
    size=4
  )

# n_cols = 6
# n_rows = 4

cell_types = list(set(adata.obs['celltypist_prediction']))
cell_types = np.reshape(
  np.resize(np.asarray(cell_types, dtype=str), 24),
  (-1,n_cols)
)
fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 3), constrained_layout=True)
  
for row_idx in range(0, n_rows):
  for col_idx in range(0, n_cols):
    cell_type_name = cell_types[row_idx, col_idx]
    current_ax = axs[row_idx, col_idx]
    current_ax.set_xlabel('xlabel', fontsize=8)
    current_ax.set_ylabel('ylabel', fontsize=8)
    plot_celltypist(adata, current_ax, cell_type_name)
    current_ax.set_title(cell_type_name, fontsize=10)

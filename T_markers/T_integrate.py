import os
import sys
import umap
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scvi
from scvi.model import TOTALVI
import matplotlib.pyplot as plt

sys.path.append('../utils_py')
import genes
import misc_utils

scvi.settings.num_threads = 6

out_dir = './out'
print(out_dir)

T_out_dir = os.path.join(out_dir, 'T_cells')
if not os.path.exists(T_out_dir):
  os.makedirs(T_out_dir)
print(T_out_dir)

model_path = os.path.join(T_out_dir, 'scvi_model_T')
print(model_path)

adata = sc.read_h5ad(os.path.join(T_out_dir, '3P_T.h5ad'))

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(
  adata,
  qc_vars=['mt'],
  percent_top=None,
  log1p=False,
  inplace=True
)

adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(
  adata,
  n_top_genes=3000,
  flavor='seurat_v3',
  batch_key='sample_name',
  subset=True,
  layer='counts',
)

TOTALVI.setup_anndata(
  adata,
  layer='counts',
  batch_key='sample_name',
  protein_expression_obsm_key='protein_expression',
	continuous_covariate_keys=['pct_counts_mt']
)
model = scvi.model.TOTALVI(adata)
model.train(max_epochs=500)
model.save(model_path, overwrite=True)

adata.obsm['X_totalVI'] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep='X_totalVI')
sc.tl.umap(adata)

for current_resolution in np.round(np.arange(0.2, 1.2, 0.1), 1):
  clustering_key = 'leiden_scVI_%0.1f' % current_resolution
  sc.tl.leiden(
    adata,
    key_added=clustering_key,
    resolution=current_resolution
  )

_, protein_means = model.get_normalized_expression(
    n_samples=25,
    include_protein_background=True,
    sample_protein_mixing=False,
    return_mean=True,
)
for i, p in enumerate(protein_means.columns):
  adata.obs[p] = protein_means.iloc[:, i]

adata.write_h5ad(os.path.join(T_out_dir, 'T_integrated.h5ad'))

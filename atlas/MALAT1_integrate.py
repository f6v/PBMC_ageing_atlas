import os
import sys
import scvi
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

out_dir = './out'
print(out_dir)

MALAT1_out_dir = os.path.join(out_dir, 'MALAT1_cells')
if not os.path.exists(MALAT1_out_dir):
  os.makedirs(MALAT1_out_dir)
print(MALAT1_out_dir)

model_path = os.path.join(MALAT1_out_dir, 'scvi_model_MALAT1')
print(model_path)

adata = sc.read_h5ad(os.path.join(out_dir, 'prepared_cleaned.h5ad'))

metadata_L1 = pd.read_csv(
  os.path.join(out_dir, 'annotated_metadata.tsv'),
  sep='\t',
  index_col=0
)
MALAT1_barcodes = metadata_L1.index[metadata_L1['cell_type_L1'].isin(['MALAT1+'])]

adata = adata[MALAT1_barcodes]

adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(
  adata,
  n_top_genes=3000,
  subset=True,
  layer='counts',
  flavor='seurat_v3',
)

scvi.model.SCVI.setup_anndata(
  adata,
  layer='counts',
  categorical_covariate_keys=['dataset', 'sample_name'],
  continuous_covariate_keys=['pct_counts_mt']
)

n_epochs = 400
model = scvi.model.SCVI(adata, n_hidden=64)
model.train(max_epochs=n_epochs, early_stopping=True)

latent = model.get_latent_representation()
adata.obsm['X_scVI'] = latent

sc.pp.neighbors(adata, use_rep='X_scVI')
sc.tl.umap(adata)

for current_resolution in np.round(np.arange(0.2, 0.5, 0.1), 1):
  clustering_key = 'leiden_scVI_%0.1f' % current_resolution
  sc.tl.leiden(
    adata,
    key_added=clustering_key,
    resolution=current_resolution
  )

adata.write(os.path.join(MALAT1_out_dir, 'batch_corrected_SCVI_MALAT1.h5ad'))

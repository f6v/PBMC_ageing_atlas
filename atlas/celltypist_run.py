import os
import scanpy as sc
import pandas as pd
import celltypist
from celltypist import models

out_dir = './out'
print(out_dir)

models.download_models(force_update = True)

adata = sc.read(os.path.join(out_dir, 'prepared_cleaned.h5ad'))

# metadata_L1 = pd.read_csv(
#   os.path.join(out_dir, 'annotated_metadata.tsv'),
#   sep='\t',
#   index_col=0
# )

# T_barcodes = metadata_L1.index[metadata_L1['cell_type_L1'].isin(['T CD4', 'T CD8'])]

# adata = adata[T_barcodes]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

predictions_immune_low = celltypist.annotate(
  adata,
  model = 'Immune_All_Low.pkl',
	# model = 'Immune_All_High.pkl',
  majority_voting = True
)
adata = predictions_immune_low.to_adata()

# adata.obs[['majority_voting']].to_csv(os.path.join(out_dir, 'celltypist_Immune_All_Low.tsv'), sep='\t', index_label='barcode')
# adata.obs[['majority_voting']].to_csv(os.path.join(out_dir, 'celltypist_Immune_All_High.tsv'), sep='\t', index_label='barcode')
adata.obs[['majority_voting']].to_csv(os.path.join(out_dir, 'celltypist_Immune_All_Low_T.tsv'), sep='\t', index_label='barcode')
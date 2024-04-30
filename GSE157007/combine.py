import os
import sys
import anndata as ad
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

sys.path.append('../utils_py')
import misc_utils
import plotting_utils

data_dir = '../../../data/immune_ageing/GSE157007/'
print(data_dir)

out_dir = './out'
print(out_dir)

metadata = pd.read_csv(
  os.path.join(out_dir, 'metadata.tsv'),
  sep='\t'
)
gex_metadata = metadata[
  metadata['assay'] == 'single-cell RNA'
].reset_index(drop=True)
all_samples = gex_metadata['sample_name'].tolist()

combined_adata = misc_utils.load_combined(out_dir, all_samples)
combined_adata = combined_adata[~combined_adata.obs.age_group.isin(['frail', 'cord_blood'])]
combined_adata.write(os.path.join(out_dir, 'combined.h5ad'))
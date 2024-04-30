import os
import sys
import anndata as ad
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

sys.path.append('../utils_py')
import misc_utils
import plotting_utils

metadata_dir = '../../../data/immune_ageing/syn22255433/'
print(metadata_dir)

out_dir = './out'
print(out_dir)

metadata = pd.read_csv(
  os.path.join(metadata_dir, 'metadata.csv'),
  sep=';'
)
all_samples = metadata['sample_name'].tolist()

combined_adata = misc_utils.load_combined(out_dir, all_samples)
combined_adata.write(os.path.join(out_dir, 'combined.h5ad'))
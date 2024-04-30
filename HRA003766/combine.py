import os
import sys
import anndata as ad
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

sys.path.append('../utils_py')
import misc_utils
import plotting_utils

data_dir = '../../../data/immune_ageing/HRA003766/'
print(data_dir)

out_dir = './out'
print(out_dir)

metadata = pd.read_csv(
  os.path.join(data_dir, 'metadata.csv'),
  sep=';'
)
all_samples = metadata['sample_name'].tolist()

combined_adata = misc_utils.load_combined(out_dir, all_samples)
combined_adata.write(os.path.join(out_dir, 'combined.h5ad'))
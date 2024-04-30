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

three_P_data_dir = '../../../data/immune_ageing/GSE164378/3P/'
print(three_P_data_dir)

out_dir = './out'
print(out_dir)

T_out_dir = os.path.join(out_dir, 'T_cells')
if not os.path.exists(T_out_dir):
  os.makedirs(T_out_dir)
print(T_out_dir)

def load_combined(data_dir):
	ab_adata = sc.read_10x_mtx(os.path.join(data_dir, 'ADT'))
	gex_adata = sc.read_10x_mtx(os.path.join(data_dir, 'GEX'))
	gex_adata.obsm['protein_expression'] = ab_adata.to_df()
	gex_adata.obsm['protein_expression'].rename(columns=lambda x: 'AB_' + x, inplace=True)
	
	metadata = pd.read_csv(os.path.join(data_dir, 'metadata.csv'), index_col=0)
	gex_adata.obs['cell_type_L2'] = metadata['celltype.l2']
	gex_adata.obs['batch'] = metadata['Batch']
	gex_adata.obs['sample_name'] = metadata['donor']

	T_cell_types = ['CD4 CTL', 'CD4 Naive', 'CD4 TCM', 'CD4 TEM', 'CD4 TEM', 'CD8 Naive', 'CD8 TCM', 'CD8 TEM', 'MAIT', 'Treg', 'gdT', 'dnT'] 
	gex_adata = gex_adata[gex_adata.obs['cell_type_L2'].isin(T_cell_types)]
	
	return gex_adata

combined_3P_data = load_combined(three_P_data_dir)

combined_3P_data.write_h5ad(os.path.join(T_out_dir, '3P_T.h5ad'))
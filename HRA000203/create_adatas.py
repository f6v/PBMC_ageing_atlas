import os
import anndata as ad
import pandas as pd
import scanpy as sc

base_data_dir = '../../../data/immune_ageing/HRA000203/'
print(base_data_dir)

cellranger_out_dir = os.path.join(base_data_dir, 'cellranger_out')
print(cellranger_out_dir)

out_dir = './out'
if not os.path.exists(out_dir):
  os.mkdir(out_dir)
print(out_dir)

metadata = pd.read_csv(
  os.path.join(base_data_dir, 'metadata.csv'),
  sep=';'
)

def create_adatas(data_dir, metadata, out_dir):
  for sample_idx in range(0, len(metadata)):
    sample_metadata = metadata.iloc[sample_idx]
    sample_name = sample_metadata['sample_name']
    print("%s" % (sample_name))

    sample_adata = sc.read_10x_mtx(
      os.path.join(data_dir, sample_name, 'count', 'sample_filtered_feature_bc_matrix'),
      cache=True
    )
    sample_adata.obs['sample_name'] = sample_name
    sample_adata.obs['age_group'] = sample_metadata['age_group']
    sample_adata.obs['sex'] = sample_metadata['sex']
    
    sample_adata.write(os.path.join(out_dir, '%s_raw.h5ad' % sample_name))

create_adatas(cellranger_out_dir, metadata, out_dir)
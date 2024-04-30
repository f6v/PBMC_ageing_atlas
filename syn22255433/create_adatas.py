import os
import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc

data_dir = '../../../data/immune_ageing/mogilenko_raw/'
print(data_dir)

metadata_dir = '../../../data/immune_ageing/Mogilenko_PBMC/'
print(metadata_dir)

out_dir = './out'
print(out_dir)

metadata = pd.read_csv(
  os.path.join(metadata_dir, 'metadata.csv'),
  sep=';'
)

def create_adatas(data_dir, metadata, out_dir):
  for sample_idx in range(0, len(metadata)):
    sample_metadata = metadata.iloc[sample_idx]
    sample_name = sample_metadata['sample_name']
    sample_dir = [x for x in os.listdir(data_dir) if sample_name in x]
    print("%s: %s" % (sample_name, sample_dir))
  
    if len(sample_dir) != 1:
      raise Exception('Should match only one dir')
    sample_adata = sc.read_10x_mtx(
      os.path.join(data_dir, sample_dir[0], 'outs', 'filtered_feature_bc_matrix'),
      cache=True
    )
    sample_adata.obs['sample_name'] = sample_name
    sample_adata.obs['age_group'] = sample_metadata['age_group']
    
    sample_adata.write(os.path.join(out_dir, '%s_raw.h5ad' % sample_name))

create_adatas(data_dir, metadata, out_dir)
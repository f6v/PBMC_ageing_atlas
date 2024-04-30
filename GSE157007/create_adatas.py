import os
import pandas as pd
import anndata as ad
import pandas as pd
import scanpy as sc

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
gex_metadata.head(3)

def create_adatas(data_dir, metadata, out_dir):
  for sample_idx in range(0, len(metadata)):
    sample_metadata = metadata.iloc[sample_idx]
    
    sample_name = sample_metadata['sample_name']
    sample_gsm = sample_metadata['gsm_name']
    print("%s: %s" % (sample_name, sample_gsm))

    sample_adata = sc.read_10x_mtx(
      os.path.join(data_dir, sample_gsm),
      cache=True
    )
    sample_adata.obs['sample_name'] = sample_name
    sample_adata.obs['age_group'] = sample_metadata['age_group']
    
    sample_adata.write(os.path.join(out_dir, '%s_raw.h5ad' % sample_name))

create_adatas(data_dir, gex_metadata, out_dir)
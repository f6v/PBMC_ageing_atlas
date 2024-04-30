import os, sys
import anndata as ad
import scanpy as sc

HRA000203_out_dir = '../HRA000203/out'
HRA000624_out_dir = '../HRA000624/out'
HRA003766_out_dir = '../HRA003766/out'
GSE157007_out_dir = '../GSE157007/out'
GSE213516_out_dir = '../GSE213516/out'
GSE214546_out_dir = '../GSE214546/out'
syn22255433_out_dir = '../syn22255433/out'

out_dir = './out'
if not os.path.exists(out_dir):
  os.makedirs(out_dir)
print(out_dir)

def load_filtered(out_dir):
  return sc.read_h5ad(os.path.join(out_dir, 'filtered.h5ad'))

adatas = {
  'HRA000203':   load_filtered(HRA000203_out_dir),
  'HRA000624':   load_filtered(HRA000624_out_dir),
	'HRA003766':   load_filtered(HRA003766_out_dir),
  'GSE157007':   load_filtered(GSE157007_out_dir),
  'GSE213516':   load_filtered(GSE213516_out_dir),
  'GSE214546':   load_filtered(GSE214546_out_dir),
  'syn22255433': load_filtered(syn22255433_out_dir)
}

combined_adata = ad.concat(adatas, label='dataset')
combined_adata.write(os.path.join(out_dir, 'old_young_1115k.h5ad'))
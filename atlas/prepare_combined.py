import os
import sys
import scanpy as sc

sys.path.append('../utils_py')
import genes
import misc_utils

out_dir = './out'
print(out_dir)

adata = sc.read_h5ad(os.path.join(out_dir, 'old_young_1115k.h5ad'))

adata = misc_utils.remove_barcodes(
  adata,
  os.path.join(out_dir, 'contaminating_barcodes.txt')
)

adata = misc_utils.remove_barcodes(
  adata,
  os.path.join(out_dir, 'doublet_barcodes.txt')
)

genes_to_remove = genes.IG_genes + genes.TCR_genes
adata = misc_utils.remove_genes(adata, genes_to_remove)

adata.write(os.path.join(out_dir, 'prepared_cleaned.h5ad'))
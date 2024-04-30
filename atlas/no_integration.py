import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import os
import sys
import numpy as np
import scanpy as sc
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="input_path")
parser.add_argument("-o", "--output", dest="output_dir")

args = parser.parse_args()
print(args)

adata = sc.read_h5ad(args.input_path)
print(adata)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
  adata,
  n_top_genes=3000,
  batch_key='dataset',
)
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'], n_jobs=16)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

adata.write(os.path.join(args.output_dir, 'not_integrated.h5ad'))

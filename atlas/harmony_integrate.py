import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import os
import sys
import numpy as np
import scanpy as sc
import scanpy.external as sce
from argparse import ArgumentParser

def run_integration(adata, out_dir):
	sce.pp.harmony_integrate(adata, 'dataset')

	adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
	sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
	sc.tl.umap(adata)

	current_out_path = os.path.join(out_dir, "harmony_integrated.h5ad")
	adata.write(current_out_path)

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="input_path")
parser.add_argument("-o", "--output", dest="output_dir")

args = parser.parse_args()
print(args)

adata = sc.read_h5ad(args.input_path)
print(adata)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=5000)
adata.raw = adata

adata = adata[:,adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

run_integration(adata.copy(), args.output_dir)

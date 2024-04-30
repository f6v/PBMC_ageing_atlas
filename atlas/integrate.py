import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import os
import sys
import scvi
import numpy as np
import scanpy as sc
from argparse import ArgumentParser

def run_integration(adata, n_hvg, out_dir):
	n_epochs = 400

	sc.pp.highly_variable_genes(
		adata,
		n_top_genes=n_hvg,
		subset=True,
		layer='counts',
		flavor='seurat_v3',
		batch_key="dataset"
	)

	scvi.model.SCVI.setup_anndata(
		adata,
		layer='counts',
		categorical_covariate_keys=['dataset', 'sample_name'],
		continuous_covariate_keys=['pct_counts_mt']
	)
	model = scvi.model.SCVI(adata)
	print(model)
	model.train(max_epochs=n_epochs, early_stopping=True)
	model_dir = os.path.join(out_dir, f"model_{n_hvg}")
	model.save(model_dir)

	latent = model.get_latent_representation()
	adata.obsm['X_scVI'] = latent

	sc.pp.neighbors(adata, use_rep='X_scVI')
	sc.tl.umap(adata)

	for current_resolution in np.round(np.arange(0.2, 0.8, 0.1), 1):
	  clustering_key = 'leiden_scVI_%0.1f' % current_resolution
	  sc.tl.leiden(
	    adata,
	    key_added=clustering_key,
	    resolution=current_resolution
	  )

	current_out_path = os.path.join(out_dir, f"integrated_{n_hvg}.h5ad")
	adata.write(current_out_path)

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="input_path")
parser.add_argument("-o", "--output", dest="output_dir")

args = parser.parse_args()
print(args)

adata = sc.read_h5ad(args.input_path)
print(adata)

adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

hvg_options = [5000]

for n_hvg in hvg_options:
	run_integration(adata.copy(), n_hvg, args.output_dir)

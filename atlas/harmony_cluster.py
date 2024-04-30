import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import os
import sys
import scvi
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

for current_resolution in np.round(np.arange(0.2, 0.9, 0.1), 1):
	clustering_key = 'leiden_harmony_%0.1f' % current_resolution
	sc.tl.leiden(
		adata,
		key_added=clustering_key,
		resolution=current_resolution
	)

current_out_path = os.path.join(args.output_dir, 'harmony_clustered.h5ad')
adata.write(current_out_path)

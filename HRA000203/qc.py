import os, sys
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('../utils_py')
import plotting_utils
import qc_utils

out_dir = './out'
print(out_dir)

combined_adata = ad.read_h5ad(os.path.join(out_dir, 'combined.h5ad'))

sc.pp.filter_genes(combined_adata, min_cells=10, inplace=True)

combined_adata.var['mt'] = combined_adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(
  combined_adata,
  qc_vars=['mt'],
  percent_top=None,
  log1p=False,
  inplace=True
)

qc_cutoffs = pd.DataFrame(
  {
    'sample_name': np.unique(combined_adata.obs['sample_name']),
    'total_counts_min': [700, 1500, 1000, 1500, 1500, 1300, 1500, 1300, 1300, 1300, 1500, 1500],
    'total_counts_max': [6000, 9000, 8000, 10000, 10000, 7000, 8000, 7000, 8000, 8000, 10000, 10000],
    'pct_counts_mt_min': [13, 9, 10, 10, 9, 10, 10, 10, 12, 10, 10, 10],
    'pct_counts_mt_max': [13, 9, 10, 10, 9, 10, 10, 10, 12, 10, 10, 10]
  }
)

plotting_utils.qc_plot(
	combined_adata,
  'total_counts',
  qc_cutoffs,
	'UMIs'
)

plotting_utils.qc_plot(
	combined_adata,
  'pct_counts_mt',
  qc_cutoffs,
	'% mt'
)

qc_df = qc_utils.get_qc_df(combined_adata, qc_cutoffs)
pd.options.display.float_format = '{:,.2f}'.format
qc_utils.format_qc_df(qc_df)

combined_adata = combined_adata[qc_df[qc_df['qc'] == 'pass']['barcode']].copy()
combined_adata.write(os.path.join(out_dir, 'filtered.h5ad'))
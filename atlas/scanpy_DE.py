import os
import sys
import re
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

out_dir = './out'
print(out_dir)

DE_out_dir = os.path.join(out_dir, 'scanpy_DE')
if not os.path.exists(DE_out_dir):
  os.mkdir(DE_out_dir)
print(DE_out_dir)

adata = sc.read_h5ad(os.path.join(out_dir, 'prepared_cleaned.h5ad'))

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

L2_metadata = pd.read_csv(
  os.path.join(out_dir, 'combined_metadata.tsv'),
  sep='\t',
  low_memory=False,
  index_col=0
)
L2_metadata.head(3)

adata.obs['cell_type_L2'] = L2_metadata['cell_type_L2']
adata.obs.head(3)

def save_to_xls(de_data, cell_type, out_path):
	save_mode = 'a' if os.path.exists(out_path) else 'w'
	with pd.ExcelWriter(out_path, mode=save_mode) as writer:  
	  de_data.to_excel(writer, sheet_name=cell_type, index=False)

def dataset_DE(adata, dataset, out_dir, min_counts=30):
	dataset_out_dir = os.path.join(out_dir, dataset)
	if not os.path.exists(dataset_out_dir):
		os.mkdir(dataset_out_dir)
	print(dataset_out_dir)

	dataset_adata = adata[adata.obs['dataset'] == dataset]
	print(f"Total cells: {len(dataset_adata)}")
	group_counts = dataset_adata.obs.groupby(['age_group', 'cell_type_L2']).size().reset_index(name='counts')
	cell_type_counts = pd.pivot(group_counts, index='cell_type_L2', columns='age_group', values='counts').reset_index()
	cell_type_counts['do_DE'] = cell_type_counts.apply(
    lambda row: row['old'] > min_counts and row['young'] > min_counts,
    axis=1
	)
	cell_types_for_DE = cell_type_counts[cell_type_counts['do_DE']]['cell_type_L2'].tolist()
  
	for cell_type in cell_types_for_DE:
		if cell_type in ['MALAT1+', 'Cycling', 'Tribo']:
			continue
		cell_type_adata = dataset_adata[dataset_adata.obs['cell_type_L2'] == cell_type].copy()
		print(f"+++ {cell_type}: {len(cell_type_adata)}")
		sc.tl.rank_genes_groups(
	    cell_type_adata,
	    groupby='age_group',
	    method='wilcoxon',
	    refence='young',
	    pts=True
	  )
		de_genes = sc.get.rank_genes_groups_df(cell_type_adata, group='old')
		de_genes = de_genes[de_genes['pvals_adj'] < 0.05]
		de_genes = de_genes[abs(de_genes['logfoldchanges']) > 0.25]
		de_genes = de_genes[de_genes['pct_nz_group'] > 0.1]
		de_genes.rename(
			columns={'names':'name', 'pvals':'pval', 'pvals_adj':'adj_pval', 'logfoldchanges':'lfc'},
			inplace=True
		)
		tsv_out_path = os.path.join(
			dataset_out_dir,
			f"{cell_type}_sig.tsv"
		)
		de_genes.to_csv(tsv_out_path, sep='\t', index=False)

dataset_DE(adata, 'GSE157007', DE_out_dir)
dataset_DE(adata, 'GSE213516', DE_out_dir)
dataset_DE(adata, 'GSE214546', DE_out_dir)
dataset_DE(adata, 'HRA000203', DE_out_dir)
dataset_DE(adata, 'HRA000624', DE_out_dir)
dataset_DE(adata, 'HRA003766', DE_out_dir)
dataset_DE(adata, 'syn22255433', DE_out_dir)

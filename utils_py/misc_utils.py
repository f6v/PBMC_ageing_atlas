import os
import scanpy as sc
import anndata as ad
import pandas as pd

def remove_genes(adata, gene_list):
  result = adata[:, ~adata.var_names.isin(list(gene_list))].copy()
  
  return result

def remove_barcodes(adata, barcode_path):
  with open(barcode_path) as f:
    barcode_list = [line.rstrip() for line in f]
    print(len(barcode_list))

  result = adata[~adata.obs.index.isin(barcode_list)].copy()
  
  return result

def load_combined(data_dir, all_samples):
  adatas = {}

  for current_sample in all_samples:
    print(current_sample)
    sample_adata = sc.read_h5ad(os.path.join(data_dir, '%s_raw.h5ad' % current_sample))

    doublet_calls = pd.read_csv(
      os.path.join(data_dir, '%s_doublets.tsv' % current_sample),
      sep='\t',
      index_col=0
    )['doublet_call']
    sample_adata.obs['doublet_class'] = doublet_calls

    adatas[current_sample] = sample_adata
  
    combined_adata = ad.concat(
      adatas,
      label='sample_name',
      index_unique='_'
    )
  
  return combined_adata

def get_cluster_markers(adata, selected_cluster):
  cluster_markers = sc.get.rank_genes_groups_df(adata, group=selected_cluster)
  cluster_markers = cluster_markers[abs(cluster_markers['logfoldchanges']) > 0.25]
  cluster_markers = cluster_markers[cluster_markers['pvals_adj'] < 0.05]
  cluster_markers = cluster_markers[cluster_markers['pct_nz_group'] > 0.1]
  #cluster_markers.sort_values('logfoldchanges', ascending=False, inplace=True)
  
  return cluster_markers
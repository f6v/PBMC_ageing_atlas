import numpy as np
import scanpy as sc

def qc_plot(adata, qc_column, qc_cutoffs, ylabel):
  ax = sc.pl.violin(
    adata,
    [qc_column],
    groupby='sample_name',
    rotation=90,
    stripplot=False,
		ylabel=ylabel,
    show=False,
		order=qc_cutoffs['sample_name'].tolist()
  )
  n_samples = len(np.unique(adata.obs['sample_name']))
  x_step = 1 / n_samples
  x_min = 0
  x_max = x_step
  ax.grid(False)

  for i in range(0, n_samples):
    y_min = qc_cutoffs.iloc[i][qc_column + '_min']
    ax.axhline(y=y_min, xmin=x_min, xmax=x_max, color='red')
    
    y_max = qc_cutoffs.iloc[i][qc_column + '_max']
    ax.axhline(y=y_max, xmin=x_min, xmax=x_max, color='red')
    ax.grid(False)

    x_min = x_min + x_step
    x_max = x_max + x_step
    
def plot_doublet_fractions(adata):
  doublet_counts = adata.obs\
    .groupby(['sample_name', 'doublet_class'])\
    .size()\
    .reset_index(name='counts')\
    .pivot(index=['sample_name'], columns='doublet_class',values='counts')\
    .reset_index('sample_name')\
    .sort_values(by='sample_name')
  doublet_counts['doublet_fraction'] = doublet_counts['doublet'] / (doublet_counts['doublet'] + doublet_counts['singlet'])
  
  doublet_counts.plot.bar(x='sample_name', y='doublet_fraction', color='#46732E99')
  
def plot_n_cells(adata):
  cell_counts = adata.obs['sample_name']\
    .value_counts()\
    .rename_axis('sample_name')\
    .reset_index(name='n_cells')\
    .sort_values(by='sample_name')
  
  cell_counts.plot.bar(x='sample_name', y='n_cells', color='#709AE199')
  
def plot_group_prop(df, col_to_group, col_to_count, x_lab):
  counts_df = df.groupby([col_to_group, col_to_count]).size().unstack()

  col_to_count_cats = list(set(df[col_to_count]))
  col_to_count_cats.sort()
  
  counts_df['total'] = counts_df.sum(axis=1)

  for current_cat in col_to_count_cats:
    counts_df[current_cat] = counts_df[current_cat] / counts_df['total']

  props_df = counts_df[col_to_count_cats]
  props_df.plot(
    kind='bar',
    stacked=True,
    ylabel='Fraction',
    xlabel=x_lab,
    fontsize=8,
    figsize=(9,3)
  ).legend(loc='center left', bbox_to_anchor=(1,0.5))
  
def subcluster_and_umap(adata, cluster_column, cluster_names, subcluster_res=0.1, pt_size=0.5):
  subcluster_column = f"{cluster_column}_{'_'.join(cluster_names)}"
  sc.tl.leiden(
    adata,
    restrict_to=(cluster_column, cluster_names),
    key_added=subcluster_column,
    resolution=subcluster_res
  )
  sc.pl.umap(
    adata,
    color=subcluster_column,
    frameon=False,
    legend_loc='on data',
    size=pt_size,
    palette='tab20'
  )
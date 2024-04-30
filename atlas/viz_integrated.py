import os
import sys
import scvi
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

sys.path.append('../utils_py')
import plotting_utils

out_dir = './out'
print(out_dir)

adata = sc.read_h5ad(os.path.join(out_dir, 'batch_corrected_SCVI_cleaned.h5ad'))

adata.uns['age_group_colors'] = np.array(['#00A087FF', '#00A087FF'])

datasets_palette = {
    'GSE157007': '#E64B35',
    'GSE213516': '#4DBBD5',
    'GSE214546': '#00A087',
    'HRA000203': '#3C5488',
		'HRA000624': '#F39B7F',
		'HRA003766': '#8491B4',
		'syn22255433': '#91D1C2'
}

# Datasets
sc.pl.umap(adata, color='dataset', groups=['GSE157007'], palette=datasets_palette, frameon=False)
sc.pl.umap(adata, color='dataset', groups=['GSE213516'], palette=datasets_palette, frameon=False)
sc.pl.umap(adata, color='dataset', groups=['GSE214546'], palette=datasets_palette, frameon=False)
sc.pl.umap(adata, color='dataset', groups=['HRA000203'], palette=datasets_palette, frameon=False)
sc.pl.umap(adata, color='dataset', groups=['HRA000624'], palette=datasets_palette, frameon=False)
sc.pl.umap(adata, color='dataset', groups=['HRA003766'], palette=datasets_palette, frameon=False)
sc.pl.umap(adata, color='dataset', groups=['syn22255433'], palette=datasets_palette, frameon=False)

# Age groups
sc.pl.umap(adata, color='age_group', groups=['young'], frameon=False)
sc.pl.umap(adata, color='age_group', groups=['old'], frameon=False)

# Doublets
sc.pl.umap(adata, color='doublet_class', groups=['doublet'], frameon=False)
sc.pl.umap(adata, color='doublet_class', groups=['singlet'], frameon=False)

sc.pl.umap(
  adata,
  color=['CMTM5', 'ITGA2B', 'PF4', 'HBA1'],
  frameon=False,
  use_raw=True,
  vmax='p99'
)

sc.pl.umap(
  adata,
  color=['PTPRC', 'CD3D', 'CD3G', 'CD3E', 'CD68', 'CD19', 'MKI67', 'TOP2A', 'NEAT1', 'MALAT1'],
  frameon=False,
  use_raw=True,
  vmax='p99'
)

sc.pl.umap(
  adata,
  color=['CD19', 'CD79A', 'MS4A1', 'JCHAIN', 'MZB1', 'ITGAX', 'IGHD', 'IGHM', 'IGHA1', 'IGHA2', 'IGHG1', 'IGHG2'],
  frameon=False,
  use_raw=True,
  vmax='p99'
)

sc.pl.umap(
  adata,
  color=['CD3D', 'CD3G', 'CD8A', 'CD8B', 'CD4', 'FOXP3', 'TRDC', 'TRGC1', 'CCR7',\
           'LEF1', 'IL7R', 'CD27', 'CD28', 'KLRB1', 'SLC4A10', 'PDCD1', 'HAVCR2', 'EOMES', 'CD34'],
  frameon=False,
  use_raw=True,
  vmax='p99'
)

sc.pl.umap(
  adata,
  color=['CD68', 'CD4', 'CD14', 'FCGR3A', 'CCR5'],
  frameon=False,
  use_raw=True,
	vmax='p99'
)

# ===== QC =====

sc.pl.umap(
  adata,
  color=['leiden_scVI_0.7'],
  frameon=False,
  legend_loc='on data',
  size=0.5,
  palette='tab20'
)

plotting_utils.subcluster_and_umap(adata, 'leiden_scVI_0.7', ['9'], 0.1)
plotting_utils.subcluster_and_umap(adata, 'leiden_scVI_0.7', ['11'], 0.1)
plotting_utils.subcluster_and_umap(adata, 'leiden_scVI_0.6', ['12'], 0.1)
plotting_utils.subcluster_and_umap(adata, 'leiden_scVI_0.7', ['4'], 0.2)

doublet_adata = adata[adata.obs['leiden_scVI_0.7_11'].isin(['9', '11,1', '18'])]
print(doublet_adata.obs.groupby('dataset').size())
doublet_barcodes = doublet_adata.obs.index.to_list()
print(len(doublet_barcodes))
with open(os.path.join(out_dir, 'doublet_barcodes.txt'), 'w') as f:
  for barcode in doublet_barcodes:
    f.write(f"{barcode}\n")
    
contaminating_adata = adata[adata.obs['leiden_scVI_0.7'].isin(['15', '16', '22'])]
print(contaminating_adata.obs[['dataset']].value_counts())
with open(os.path.join(out_dir, 'contaminating_barcodes.txt'), 'w') as f:
  for barcode in contaminating_adata.obs.index.to_list():
    f.write(f"{barcode}\n")

# ===== END QC =====


# ===== Cell Types =====

all_markers = ['PTPRC', 'CD3D', 'CD3G', 'CD8A', 'CD4',
               'IL2RA', 'FOXP3',
               'CTLA4', 'PDCD1',
               'CD27', 'CD28',
               'IL7R', 'CCR7', 'SELL', 'CD40LG', 'CD44', 'LEF1',
               'CD68', 'FCGR3A', 'CD14', 'CCR5',
               'CD19', 'CD79A', 'MS4A1', 'ITGAX', 'CD93', 
							 'CD24', 'CD38',
               'IGHD', 'IGHM', 'IGHA1', 'IGHA2',
               'JCHAIN', 'MZB1',
               'CD1C', 'FCER1A', 'CLEC10A',
               'GZMB', 'GZMK',
               'GNLY', 'XCL2', 'NKG7',
               'NCAM1', 'CCL5',
               'KLRB1', 'SLC4A10',
               'TRDC', 'TRGC1', 'TRGC2',
               'MALAT1', 'NEAT1',
               'MKI67', 'TOP2A',
               'PF4', 'CD34'
              ]
sc.pl.dotplot(
  adata,
  all_markers,
  groupby='leiden_scVI_0.7',
  dendrogram=False,
  cmap='Reds',
  swap_axes=False,
  use_raw=True,
  standard_scale='var',
)

plotting_utils.subcluster_and_umap(adata, 'leiden_scVI_0.6', ['2'], 0.1)
plotting_utils.subcluster_and_umap(adata, 'leiden_scVI_0.6', ['4'], 0.1)
plotting_utils.subcluster_and_umap(adata, 'leiden_scVI_0.6', ['12'], 0.1)
plotting_utils.subcluster_and_umap(adata, 'leiden_scVI_0.6', ['13'], 0.1)
plotting_utils.subcluster_and_umap(adata, 'leiden_scVI_0.6', ['15'], 0.1)

adata.obs['cell_type_L1'] = 'other'

adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6_15'] == '15,1', 'PC', adata.obs['cell_type_L1'])

adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6_13'] == '13,2', 'CD34+', adata.obs['cell_type_L1'])

adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6_15'] == '15,0', 'Cycling', adata.obs['cell_type_L1'])

adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6_2']  == '2,1', 'DC', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6']    == '14', 'pDC', adata.obs['cell_type_L1'])

adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6_13'] == '13,0', 'NK CD56-high', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6_13'] == '13,1', 'NK CD56-low', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6_4']  == '4,0', 'NK CD56-low', adata.obs['cell_type_L1'])

adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6']    == '9', 'B memory', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6']    == '6', 'B naive', adata.obs['cell_type_L1'])

adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6']    == '0', 'T CD4', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6']    == '1', 'T CD4', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6']    == '8', 'T CD4', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6_12'] == '12,0', 'T CD4', adata.obs['cell_type_L1'])

adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6']    == '5', 'T CD8', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6']    == '3', 'T CD8', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6']    == '7', 'T CD8', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6_4']  == '4,1', 'T CD8', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6_12'] == '12,1', 'T CD8', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6_12'] == '12,2', 'T CD8', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6_12'] == '12,3', 'T CD8', adata.obs['cell_type_L1'])

adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6']    == '11', 'MALAT1+', adata.obs['cell_type_L1'])

adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6']    == '10', 'Mono non-classical', adata.obs['cell_type_L1'])
adata.obs['cell_type_L1'] = np.where(adata.obs['leiden_scVI_0.6_2']  == '2,0', 'Mono classical', adata.obs['cell_type_L1'])

sc.pl.umap(
    adata,
    color=['cell_type_L1'],
    frameon=False,
    size=1,
    palette={
	    'B memory': '#A6DBBB',
      'B naive': '#80C39E',
			'B atypical': '#5AAC82',
      'PC': '#359566',
			
      'CD34+': '#C6DBEF',
      'T CD4': '#A6C5DD',
      'T CD8': '#86B0CB',

      'NK CD56-high': '#ECD9F1',
      'NK CD56-low': '#D6C1E8',

      'Mono classical': '#F3E55C',
      'Mono non-classical': '#EFB84C',
      'DC': '#EB8C3C',
      'pDC': '#E8602D',

      'MALAT1+': '#505050',
      'other': '#ECD9F1',
      'Cycling': '#A47786',
    },
    save='UMAP_cell_type_L1.png'
)

cell_type_markers =\
  ['CD3D', 'CD8A',
   'NCAM1', 'GZMK', 'NKG7',
	 'IL2RA',
   'MKI67',
   'CD79A', 'IGHD', 'IGHM', 'IGHA1', 'IGHA2',
   'JCHAIN', 'MZB1',
   'CD4',
   'CD68', 'FCGR3A', 'CD14',
   'CD1C', 'FCER1A',
	 'CD34',
	 'MALAT1'
  ]

sc.tl.dendrogram(adata, groupby='cell_type_L1', use_rep='X_scVI')

sc.pl.dotplot(
  adata,
  cell_type_markers,
  groupby='cell_type_L1',
  cmap='Reds',
  swap_axes=False,
  dendrogram=True,
  use_raw=True,
  standard_scale='var',
	save='cell_type_L1_dotplot.png'
)



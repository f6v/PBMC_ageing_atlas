import os
import pandas as pd
import anndata as ad
import scanpy as sc

data_dir = '../../../data/immune_ageing/GSE214546/'
print(data_dir)

out_dir = './out'
if not os.path.exists(out_dir):
	os.mkdir(out_dir)
print(out_dir)

metadata = pd.read_csv(
  os.path.join(data_dir, 'metadata.csv'),
  sep=';'
)

all_files = os.listdir(data_dir)
def create_adatas(data_dir, metadata, out_dir):
	for sample_idx in range(0, len(metadata)):
		sample_metadata = metadata.iloc[sample_idx]

		sample_name = sample_metadata['sample_name']
		sample_age = sample_metadata['age']
		sample_age_group = sample_metadata['age_group']
		sample_sex = sample_metadata['sex']

		sample_file = [file_name for file_name in all_files if sample_name in file_name]
		assert len(sample_file) == 1
		print("%s: %s" % (sample_name, sample_file[0]))

		if sample_age_group == 'old' and sample_age < 60:
			print("Skipping %s" % (sample_name))
			continue
		sample_adata = sc.read_10x_h5(os.path.join(data_dir, sample_file[0]))
		sample_adata.var_names_make_unique()

		sample_adata.obs['sample_name'] = sample_name
		sample_adata.obs['age_group'] = sample_age_group
		sample_adata.obs['sex'] = sample_sex

		sample_adata.write(os.path.join(out_dir, '%s_raw.h5ad' % sample_name))

create_adatas(data_dir, metadata, out_dir)
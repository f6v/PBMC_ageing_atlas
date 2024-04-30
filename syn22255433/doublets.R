suppressPackageStartupMessages({
library(scDblFinder)
library(DropletUtils)
library(tidyverse)
library(Matrix)
})

source('../doublet_utils.R')

metadata_dir <- '../../../data/immune_ageing/Mogilenko_PBMC/'
print(metadata_dir)

data_dir <- '../../../data/immune_ageing/mogilenko_raw'
print(data_dir)

out_dir <- './out'
print(out_dir)

metadata <- file.path(metadata_dir, 'metadata.csv') |>
  read_delim(show_col_types = F)

for(row_idx in 1:nrow(metadata)) {  
  current_sample <- metadata[[row_idx, 'sample_name']]
  sample_folder <- list.files(data_dir, pattern = current_sample)
  
  print(paste0(current_sample, ': ', sample_folder))
  flush.console()
  
  current_out_path <- file.path(out_dir, paste0(current_sample, '_doublets.tsv'))
  file.path(data_dir, sample_folder, 'outs', 'filtered_feature_bc_matrix') |>
    read10xCounts() |>
    process_sample(current_out_path)
}
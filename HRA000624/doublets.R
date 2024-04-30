suppressPackageStartupMessages({
library(scDblFinder)
library(DropletUtils)
library(tidyverse)
})
source('../doublet_utils.R')


data_dir <- '../../../data/immune_ageing/HRA000624'
print(data_dir)

cellranger_out_dir <- file.path(data_dir, 'cellranger_out')
print(cellranger_out_dir)

out_dir <- './out'
print(out_dir)

metadata <- file.path(data_dir, 'metadata.csv') |>
  read_delim(show_col_types = F)

for(row_idx in 1:nrow(metadata)) {  
  current_accession <- metadata[[row_idx, 'sample_accession']]
  current_sample <- metadata[[row_idx, 'sample_name']]
  
  print(paste(current_sample, '-', current_accession))
  flush.console()
  
  current_out_path <- file.path(out_dir, paste0(current_sample, '_doublets.tsv'))
  file.path(cellranger_out_dir, current_accession, 'outs', 'filtered_feature_bc_matrix') |>
    read10xCounts() |>
    process_sample(current_out_path)
}
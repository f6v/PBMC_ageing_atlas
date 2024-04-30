suppressPackageStartupMessages({
library(scDblFinder)
library(DropletUtils)
library(tidyverse)
})

source('../doublet_utils.R')

data_dir <- '../../../data/immune_ageing/GSE214546/'
print(data_dir)

out_dir <- './out'
print(out_dir)

metadata <- file.path(data_dir, 'metadata.csv') |>
  read_delim(show_col_types = F)

for(row_idx in 1:nrow(metadata)) {  
  current_sample <- metadata[[row_idx, 'sample_name']]

  print(current_sample)
  flush.console()
  
  current_out_path <- file.path(out_dir, paste0(current_sample, '_doublets.tsv'))
  file.path(data_dir, list.files(data_dir, pattern = current_sample)) |>
    read10xCounts() |>
    process_sample(current_out_path)
}
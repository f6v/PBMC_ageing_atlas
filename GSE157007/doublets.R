suppressPackageStartupMessages({
library(scDblFinder)
library(DropletUtils)
library(tidyverse)
})

source('../doublet_utils.R')

data_dir <- '../../../data/immune_ageing/GSE157007'
print(data_dir)

out_dir <- './out'
print(out_dir)

metadata <- file.path(out_dir, 'metadata.tsv') |>
  read_delim(show_col_types = F) |>
  filter(assay =='single-cell RNA')
metadata |> head(3)

for(row_idx in 1:nrow(metadata)) {  
  current_sample <- metadata[[row_idx, 'sample_name']]
  current_gsm_name <- metadata[[row_idx, 'gsm_name']]

  print(paste(current_sample, '-', current_gsm_name))
  flush.console()
  
  current_out_path <- file.path(out_dir, paste0(current_sample, '_doublets.tsv'))
  file.path(data_dir, current_gsm_name) |>
    read10xCounts() |>
    process_sample(current_out_path)
}
suppressPackageStartupMessages({
  library(tidyverse)
})

out_dir <- './out'
print(out_dir)

all_metadata <- file.path(out_dir, 'annotated_metadata.tsv') |>
  read_tsv(show_col_types = F) |>
  rename(barcode = `...1`) |>
  mutate(cell_type_L2 = cell_type_L1) |>
  select(barcode, sample_name, age_group, cell_type_L2, dataset)
all_metadata |> head(3)

T_metadata <- file.path(out_dir, 'annotated_T_metadata.tsv') |>
  read_tsv(show_col_types = F) |>
  rename(barcode = `...1`) |>
  select(barcode, sample_name, age_group, cell_type_L2, dataset)
  #mutate(cell_type_L3 = cell_type_L2)
T_metadata |> head(3)

all_except_T_metadata <- all_metadata |>
  filter(!(barcode %in% T_metadata$barcode))
nrow(all_except_T_metadata)

combined_metadata <- all_except_T_metadata |>
  bind_rows(T_metadata)

combined_metadata |>
  write_tsv(file.path(out_dir, 'combined_metadata.tsv'))

suppressPackageStartupMessages({
library(tidyverse)
library(openxlsx)
library(RColorBrewer)
})

out_dir <- './out'
print(out_dir)

metadata_dir <- '../../../data/immune_ageing/metadata'
print(metadata_dir)

load_metadata <- function(metadata_path) {
  metadata_path |>
    read_tsv(show_col_types = F) |>
    rename(barcode = `...1`) |>
    mutate(
      age_group = ifelse(age_group == 'young', 'Young', 'Old')
    )
}

save_to_xls <- function(data, sheet_name, out_path) {
  workbook <- createWorkbook()
  addWorksheet(workbook, sheetName = sheet_name, gridLines = T)
  writeDataTable(workbook, sheet = sheet_name, tableStyle = 'TableStyleLight1', data)
  saveWorkbook(workbook, out_path, overwrite = T)
  
  return(data)
}

cell_type_colors <- c('#A6DBBB', '#6DB890', '#359566', '#5AAC82', '#C6DBEF', '#9FC1D9', '#79A7C4', '#538DAE', '#2D7399', '#C1AADF', '#075A84', '#F3E55C',
         '#EFB84C', '#EB8C3C', '#E8602D', '#505050', '#ECD9F1', '#967BCE')
cell_type_names <- c('B memory', 'B naive', 'B atypical', 'PC', 'CD34+', 'T CD4', 'T CD8', 'NK CD56-high', 'NK CD56-low', 'NK CD25+','T CTLA4+',
                'Mono classical', 'Mono non-classical', 'DC', 'pDC', 'MALAT1+', 'other', 'Cycling')
names(cell_type_colors) <- cell_type_names

L1_scvi_metadata <- file.path(out_dir, 'annotated_metadata.tsv') |>
  load_metadata()
L1_harmony_metadata <- file.path(out_dir, 'harmony_annotated_metadata.tsv') |>
  load_metadata()


L1_scvi_metadata |>
  group_by(dataset) |>
  count(cell_type_L1) |>
  mutate(prop = n / sum(n)) |>
  mutate(cell_type_L1 = factor(cell_type_L1, levels = cell_type_names)) |>
  save_to_xls('PBMCs', file.path(out_dir, 'PBMC_dataset_props.xlsx')) |>
  ggplot(aes(x = dataset, y = prop, fill = cell_type_L1)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = cell_type_colors) +
  theme_minimal() +
  labs(x = '', y = 'Proportion', fill = 'Cell type') +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

all_samples <- L1_scvi_metadata |>
  filter(cell_type_L1 == 'B memory') |>
  distinct(sample_name, dataset) |>
  arrange(dataset) |>
  pull(sample_name)

L1_harmony_metadata |>
  group_by(dataset) |>
  count(cell_type_L1) |>
  mutate(prop = n / sum(n)) |>
  mutate(cell_type_L1 = factor(cell_type_L1, levels = cell_type_names)) |>
  save_to_xls('PBMCs', file.path(out_dir, 'PBMC_harmony_dataset_props.xlsx')) |>
  ggplot(aes(x = dataset, y = prop, fill = cell_type_L1)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = cell_type_colors) +
  theme_minimal() +
  labs(x = '', y = 'Proportion', fill = 'Cell type') +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

L2_T_metadata <- file.path(out_dir, 'annotated_T_metadata.tsv') |>
  load_metadata()

T_cell_type_colors <- c('#A6DBBB', '#5AAC82', '#80C39E', '#359566', '#C6DBEF', '#A6C5DD', '#86B0CB', '#D6C1E8', '#ECD9F1',
                        '#C1AADF', '#F3E55C', '#505050', '#C0C0C0')
T_cell_type_names <- c('CD4 Tn', 'CD4 Tcm', 'CD8 Tn', 'CD8 Tcm', 'CD4 CTL', 'CD4 Tem', 'CD8 Tem', 'MAIT', 'Tgd',
                       'NKT', 'Treg', 'Tdn', 'Tribo')
names(T_cell_type_colors) <- T_cell_type_names

L2_T_metadata |>
  group_by(dataset) |>
  count(cell_type_L2) |>
  mutate(prop = n / sum(n)) |>
  mutate(cell_type_L2 = factor(cell_type_L2, levels = T_cell_type_names)) |>
  save_to_xls('T cells', file.path(out_dir, 'T_dataset_props.xlsx')) |>
  ggplot(aes(x = dataset, y = prop, fill = cell_type_L2)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = T_cell_type_colors) +
  theme_minimal() +
  labs(x = '', y = 'Proportion', fill = 'Cell type') +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

B_cell_type_colors <- c('#5AAC82', '#A6DBBB', '#80C39E')
B_cell_type_names <- c('B atypical', 'B memory', 'B naive')
names(B_cell_type_colors) <- B_cell_type_names

L2_B_metadata <- file.path(out_dir, 'annotated_B_metadata.tsv') |>
  load_metadata()

L2_B_metadata |>
  group_by(dataset) |>
  count(cell_type_L2) |>
  mutate(prop = n / sum(n)) |>
  mutate(cell_type_L2 = factor(cell_type_L2, levels = B_cell_type_names)) |>
  save_to_xls('B cells', file.path(out_dir, 'B_dataset_props.xlsx')) |>
  ggplot(aes(x = dataset, y = prop, fill = cell_type_L2)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = B_cell_type_colors) +
  theme_minimal() +
  labs(x = '', y = 'Proportion', fill = 'Cell type') +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

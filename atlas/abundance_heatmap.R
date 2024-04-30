suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(ggsci)
  library(ComplexHeatmap)
  library(circlize)
  library(viridis)
})

out_dir <- './out'
print(out_dir)

cell_metadata <- file.path(out_dir, 'combined_metadata.tsv') |>
  read_tsv(show_col_types = F)

sample_metadata <- cell_metadata |>
  filter(!(cell_type_L2 %in% c('other'))) |>
  select(sample_name, age_group, dataset) |>
  distinct(sample_name, .keep_all = T) |>
  arrange(dataset, age_group)

props_df <- cell_metadata |>
  filter(!(cell_type_L2 %in% c('other'))) |>
  group_by(sample_name) |>
  dplyr::count(cell_type_L2) |>
  mutate(prop = n / sum(n)) |>
  select(-n) |>
  ungroup() |>
  pivot_wider(names_from = c('sample_name'), values_from = 'prop', values_fill = 0) |>
  as.data.frame()
rownames(props_df) <- props_df$cell_type_L2
props_df <- props_df |>
  select(-cell_type_L2)
# props_df
props_df <- props_df[, sample_metadata$sample_name]

props_matrix <- matrix(as.numeric(unlist(props_df)), nrow = nrow(props_df))
rownames(props_matrix) <- rownames(props_df)

options(repr.plot.width=16, repr.plot.height=8)

ha = HeatmapAnnotation(
  Dataset = sample_metadata$dataset,
  `Age group` = sample_metadata$age_group,
  col = list(
    Dataset = c(
      'GSE157007'   = '#E64B35',
      'GSE213516'   = '#4DBBD5',
      'GSE214546'   = '#00A087',
      'HRA000203'   = '#3C5488',
      'HRA000624'   = '#F39B7F',
      'HRA003766'   = '#8491B4',
      'syn22255433' = '#91D1C2'
    ),
    `Age group` = c('old' = '#F39B7F', 'young' = '#91D1C2')
  )
)

prop_heatmap <- Heatmap(
  props_matrix,
  cluster_rows = F,
  cluster_columns = F,
  column_split = sample_metadata$dataset,
  row_title_rot = 0,
  border = T,
  column_names_rot = 0,
  column_names_centered = T,
  top_annotation = ha,
  name = 'Fraction',
  heatmap_legend_param = list(
    legend_direction = 'horizontal',
    legend_width = unit(6, 'cm')
  )
)

draw(prop_heatmap, heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom', merge_legend = T)
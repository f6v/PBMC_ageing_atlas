cutoffs_plot <- function(seurat_obj, feature_column, feature_label, cutoffs_df, log = T) {
  base_plot <- seurat_obj[[]] |>
    ggplot(aes(x = 1, y = .data[[feature_column]], fill = .data[['sample_name']])) +
      geom_violin() +
      facet_grid(cols = vars(sample_name)) +
      theme_bw() +
      geom_hline(linetype = 'dashed', aes(yintercept = low), cutoffs_df) +
      geom_hline(linetype = 'dashed', aes(yintercept = high), cutoffs_df) +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none'
      ) +
      labs(y = feature_label)
  
  if(log) {
    base_plot <- base_plot + scale_y_log10()
  }
  return(base_plot)
}

feature_plot <- function(seurat_obj, features, ncol = 4) {
  FeaturePlot(
    seurat_obj,
    features = features,
    ncol = ncol,
    order = F,
    cols = plasma(11)
  ) & NoLegend()
}

plot_umap_dotplot <- function(seurat_obj, group_column, selected_features) {
  umap_plot <- DimPlot(
    seurat_obj,
    group.by = group_column,
    label = T
  ) + NoLegend()

  dotplot_plot <- DotPlot(
    seurat_obj,
    group.by = group_column,
    features = selected_features,
    cluster.idents = T
  ) +
  scale_color_gradientn(
    colours = colorRampPalette(rev(brewer.pal(8, 'Spectral')))(100)
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 16)
  )
  
  return(
    umap_plot + dotplot_plot + plot_layout(widths = c(1, 2)) & theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
  )
}

highlight_cells <- function(seurat_obj, column_name, column_value) {
  cells <- seurat_obj[[]] |>
    filter(!!sym(column_name) == !!column_value) |>
    rownames()
  n_cells <- length(cells)
  DimPlot(
    seurat_obj,
    cells.highlight = cells,
    cols.highlight = '#59A14F',
    sizes.highlight = 0.05,
    pt.size = 0.05
  ) +
    ggtitle(paste0(column_value, '(', n_cells, ')')) +
    NoLegend()
}


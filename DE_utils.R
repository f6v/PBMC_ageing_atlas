plot_volcano <- function(de_result) {
  is_de_mask <- de_result$p_val_adj < 0.05 & abs(de_result$avg_log2FC) > 0.25
  is_up_mask <- de_result$p_val_adj < 0.05 & de_result$avg_log2FC > 0.25
  is_down_mask <- de_result$p_val_adj < 0.05 & de_result$avg_log2FC < -0.25
  
  de_result$is_de <- 'Not changed'
  de_result[is_up_mask, 'is_de'] <- 'Up'
  de_result[is_down_mask, 'is_de'] <- 'Down'
  de_result$label <- NA
  de_result$label[is_de_mask] <- de_result$gene[is_de_mask]
  
  color_values <- c('#415484', 'grey60', '#D55640')
  names(color_values) <- c('Down', 'Not changed', 'Up')
  
  res_plot <- de_result|>
    ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), label = label, col = is_de)) +
    geom_point() +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = 'dashed', alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', alpha = 0.5) +
    theme_minimal() +
    theme(legend.position = 'top') +
    geom_text_repel(show.legend = F, max.overlaps = 20) +
    scale_color_manual(values = color_values) +
    guides(col = guide_legend(title = element_blank())) +
    labs(x = 'LFC', y = '-log10(p)')
  res_plot
}

plot_enrichment <- function(enrichment_df) {
  enrichment_df <- enrichment_df |>
    filter(source == 'GO:BP')
    
  if(nrow(enrichment_df) == 0) {
    return(plot_spacer()) 
  } else {
    return(
      enrichment_df |>
        mutate(minus_log_p_value = -log10(p_value), combined_title = paste(term_id, term_name)) |>
        slice_max(minus_log_p_value, n = 20) |>
        mutate(combined_title = str_trunc(combined_title, 50)) |>
        ggplot(aes(x = minus_log_p_value, y  = reorder(combined_title, minus_log_p_value))) +
        geom_bar(stat='identity', fill = '#E79F84') +
        scale_y_discrete(labels = function(x) sub('\\s', '\n', x) ) +
        ggtitle('GO:BP enriched in up DE') +
        theme_minimal() +
        labs(x = '-log10(p)', y = '') 
    )
  }
}

condition_de_test <- function(seurat_obj, cell_type_column, cell_type_value) {
  cell_type_barcodes <- seurat_obj[[]] |>
    filter(!!sym(cell_type_column) %in% cell_type_value) |>
    rownames()
  cell_type_seurat <- seurat_obj[, cell_type_barcodes]
  
  Idents(cell_type_seurat) <- 'age'

  all_markers <- FindMarkers(
    cell_type_seurat,
    ident.1 = 'O',
    #test.use = 'DESeq2',
    verbose = F
  ) |>
  drop_na()
  all_markers$gene <- rownames(all_markers)
  
  return(
    all_markers
  )
}

plot_de_test <- function(all_markers) {
  volcano_plot <- plot_volcano(all_markers)
  
  sig_markers <- all_markers |>
    filter(p_val_adj < 0.05 & avg_log2FC > 0)
  enrichment_res <- gost(query = sig_markers$gene, organism = 'mmusculus')
  go_bp_plot <- plot_enrichment(enrichment_res$result)

  return(
    volcano_plot + go_bp_plot + plot_layout(widths = c(3, 2))
  )
}
init_dir <- function(new_dir) {
  if (!dir.exists(new_dir)) {
	dir.create(new_dir)
  }
}

no_integration_workflow <- function(seurat_obj, n_pcs, clustering_resolutions, vars_to_regress = NULL) {
  seurat_obj |>
	SCTransform(method = 'glmGamPoi', vst.flavor = 'v2', vars.to.regress = vars_to_regress, verbose = F) |>
	RunPCA(verbose = F) |>
	RunUMAP(dims = 1:n_pcs, verbose = F) |>
	FindNeighbors(dims = 1:n_pcs, verbose = F) |>
	FindClusters(verbose = F, resolution = clustering_resolutions)
}

get_raw_groups <- function(seurat_obj, grouping_column, group_values) {
  DefaultAssay(seurat_obj) <- 'RNA'

  selected_barcodes <- seurat_obj[[]] |>
    filter(!!sym(grouping_column) %in% group_values) |>
    rownames()
  subset_seurat <- seurat_obj[, selected_barcodes]
  
  return(
    subset_seurat |>
      DietSeurat(assay = 'RNA')
  )
}


integration_workflow <- function(seurat_obj, split_column, n_pcs, clustering_resolutions,
                                 features_to_remove = NULL, vars_to_regress = NULL) {
  batch_objects <- seurat_obj |>
	SplitObject(split_column) |>
	lapply(function(batch_seurat) {
	  SCTransform(
		batch_seurat,
		vst.flavor = 'v2',
		method = 'glmGamPoi',
		vars.to.regress = vars_to_regress,
		verbose = F
	  ) |>
	  RunPCA(npcs = n_pcs, verbose = F)
	})

  features <- SelectIntegrationFeatures(batch_objects, nfeatures = 2000, verbose = F)
  if(!is.null(features)) {
    print('Filter features')
    print(length(features))
    features <- features[!(features %in% features_to_remove)]
    print(length(features))
  }
  batch_objects <- PrepSCTIntegration(batch_objects, anchor.features = features, verbose = F)
  anchors <- FindIntegrationAnchors(
    batch_objects, 
    normalization.method = 'SCT',
    reduction = 'rpca',
    anchor.features = features,
    verbose = F
  )
  seurat_integrated <- IntegrateData(anchors, normalization.method = 'SCT', verbose = F) |>
	RunPCA(verbose = F) |>
	RunUMAP(reduction = 'pca', dims = 1:n_pcs, verbose = F) |>
	FindNeighbors(reduction = 'pca', dims = 1:n_pcs, verbose = F) |>
	FindClusters(resolution = clustering_resolutions, verbose = F)
  
  return(seurat_integrated)
}

find_all_markers <- function(seurat_obj, selected_idents) {
  Idents(seurat_obj) <- selected_idents

  return(
	seurat_obj |>
	  FindAllMarkers(verbose = F)
  )
}

find_top_markers <- function(all_markers) {
  return(
	all_markers |>
      filter(p_val_adj < 0.05) |>
	  group_by(cluster) |>
	  slice_max(n = 20, order_by = avg_log2FC)
  )
}

add_metadata <- function(seurat_obj, column) {
  names(celltypist_predictions) <- celltypist_result$barcode
  seurat_obj <- seurat_obj |>
	AddMetaData(celltypist_predictions, 'celltypist_prediction')
  return(seurat_obj)
}
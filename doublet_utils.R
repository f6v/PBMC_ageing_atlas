make_final_calls <- function(all_results, n_runs) {
  all_results['doublet_call'] <- apply(all_results, 1, function(x) {
    all_calls <- x[paste0('run_', 1:n_runs)]
    n_doublets <- length(all_calls[all_calls == 'doublet'])

    return(
      ifelse(n_doublets > n_runs / 2, 'doublet', 'singlet')
    )
  })
  return(all_results)
}

process_sample <- function(current_sce, result_path, n_runs = 10) {
  run_results <- data.frame(barcode = colData(current_sce)$Barcode)

  for(run_idx in 1:n_runs) {
    doublets_sce <- scDblFinder(current_sce, verbose = F)
    run_results[paste0('run_', run_idx)] <- colData(doublets_sce)$scDblFinder.class
  }
  final_result <- make_final_calls(run_results, n_runs)
  final_result |>
    write_delim(result_path, delim = '\t')
}
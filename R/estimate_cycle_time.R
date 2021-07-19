estimate_cycle_time <- function(exprs, ensembl_ids=NULL,
                                handle_multiple_observations=c("remove", "mean", "median", "max"),
                                quiet=FALSE) {

  handle_multiple_observations <- match.arg(handle_multiple_observations)

  # TODO: Check input exprs has sample names as colnames

  # Get common genes
  if (is.null(ensembl_ids)) {
    common_genes <- intersect(rownames(coef_matrix), rownames(exprs))
    n_obs <- length(common_genes)
  } else {
    common_genes <- intersect(ensembl_ids, rownames(exprs))
    n_obs <- sum(rownames(exprs) %in% common_genes)
  }
  message(paste0("Using ", n_obs, "/", nrow(exprs), " observations from the input expression matrix."))

  # TODO: Check if provided gene IDs look like ensembl IDs

  # TODO: Handle provided ensembl IDs many-to-one observations

  # Process expression data
  normalised_exprs <- preprocessCore::normalize.quantiles.use.target(exprs[common_genes,], target=ref)
  colnames(normalised_exprs) <- colnames(exprs)
  rownames(normalised_exprs) <- common_genes

  # Get expected expression matrix
  expected_exprs <- coef_matrix %*% lp_matrix

  # Get the mean squared error
  mse <- sapply(1:ncol(exprs), function(i) {
    colMeans(sweep(x=expected_exprs[common_genes,], MARGIN=1, STATS=normalised_exprs[,i], FUN="-") ^ 2)
  })
  colnames(mse) <- colnames(exprs)

  # Estimate cycle time using the minimum MSE
  estimated_time <- rownames(mse)[apply(mse, 2, which.min)]
  estimated_time <- as.numeric(sub("^time_", "", estimated_time))
  names(estimated_time) <- colnames(mse)

  # TODO: Return an S4 object?
  return(estimated_time)
}

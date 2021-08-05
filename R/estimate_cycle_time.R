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
    common_genes <- intersect(ensembl_ids, rownames(coef_matrix))
    n_obs <- sum(rownames(exprs) %in% common_genes)
  }
  if (! quiet) {
    message(paste0("Using ", n_obs, "/", nrow(exprs), " observations from the input expression matrix."))
  }

  # TODO: Check if provided gene IDs look like ensembl IDs

  # TODO: Handle provided ensembl IDs many-to-one observations
  if (is.null(ensembl_ids)) {
    common_exprs <- exprs[common_genes,]
  } else if (! any(duplicated(ensembl_ids))) {
    rownames(exprs) <- ensembl_ids
    common_exprs <- exprs[common_genes,]
  } else {
    # TODO: Refactor this later
    probe_to_ensembl <- ensembl_ids
    rownames(probe_to_ensembl) <- rownames(exprs)

    # Handle one-to-one observations
    single_ensembl <- names(which(table(ensembl_ids) == 1))
    single_probes <- names(probe_to_ensembl %in% single_ensembl)
    single_exprs <- exprs[single_probes,]
    rownames(single_exprs) <- rownames(single_exprs)[probe_to_ensembl]

    # Consolidate ensembl IDs that have multiple probes
    multiple_ensembl <- names(which(table(ensembl_ids) > 1))
    consolidated_list <- list()
    for (ensembl_id in multiple_ensembl) {
      illumina_ids <- names(which(probe_to_ensembl == ensembl_id))
      # TODO: Currently only using means
      consolidated_list[[ensembl_id]] <- colMeans(exprs[illumina_ids,])
    }

    # Combine single and multiple matrices
    common_exprs  <- rbind(
      do.call(rbind, consolidated_list),
      single_exprs
    )
    common_exprs <- common_exprs[common_genes,]
  }

  # Process expression data
  normalised_exprs <- preprocessCore::normalize.quantiles.use.target(common_exprs[common_genes,], target=ref)
  colnames(normalised_exprs) <- colnames(common_exprs)
  rownames(normalised_exprs) <- common_genes

  # Get expected expression matrix
  expected_exprs <- coef_matrix %*% lp_matrix

  # Get the mean squared error
  mse <- sapply(1:ncol(common_exprs), function(i) {
    colMeans(sweep(x=expected_exprs[common_genes,], MARGIN=1, STATS=normalised_exprs[,i], FUN="-") ^ 2)
  })
  colnames(mse) <- colnames(common_exprs)

  # Estimate cycle time using the minimum MSE
  estimated_time <- rownames(mse)[apply(mse, 2, which.min)]
  estimated_time <- as.numeric(sub("^time_", "", estimated_time))
  names(estimated_time) <- colnames(mse)

  # TODO: Return an S4 object?
  return(list(estimated_time=estimated_time,
              mse=mse))
}

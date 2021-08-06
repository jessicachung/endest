#' Title
#'
#' @param exprs
#' @param ensembl_ids
#' @param entrez_ids
#' @param handle_multiple_observations
#' @param quiet
#'
#' @return
#' @export
#'
#' @examples
estimate_cycle_time <- function(exprs, ensembl_ids=NULL, entrez_ids=NULL,
                                handle_multiple_observations=c("remove", "mean", "median", "max"),
                                quiet=FALSE) {

  handle_multiple_observations <- match.arg(handle_multiple_observations)

  # Check input exprs has sample names as colnames
  if (is.null(colnames(exprs))) {
    stop("Error: Missing sample names. Expression matrix must have column names")
  }

  # TODO: Handle cases with missing rownames?

  # Determine what gene identifiers are used
  # TODO: Check if provided gene IDs look like ensembl/entrez IDs
  if (! is.null(ensembl_ids)) {
    gene_id_type <- "ensembl"
    gene_ids <- ensembl_ids
  } else if (! is.null(entrez_ids)) {
    gene_id_type <- "entrez"
    gene_ids <- as.character(entrez_ids)
  } else {
    gene_id_type <- identify_id_type(rownames(exprs))
    if (! quiet) {
      message(paste0("Using ", gene_id_type, " identifiers from rownames."))
    }
  }

  # Check if provided ID lengths are the same as the number of rows in the the expression matrix
  if (! is.null(ensembl_ids) | ! is.null(entrez_ids)) {
    if (length(gene_ids) != nrow(exprs)) {
      stop("Error: The number of provided gene IDs does not match the number of rows in the expression matrix.")
    }
  }

  if (gene_id_type == "entrez") {
    # Using entrez IDs means using fewer genes since not all ensembl IDs in the endest model have entrez IDs
    coef_matrix <- coef_matrix[names(ensembl_to_entrez),]
    rownames(coef_matrix) <- ensembl_to_entrez
  }

  # TODO: handle NAs in ensembl_ids and entrez_ids args
  # TODO: handle possible duplicate rownames (treat as probes?)

  # Get common genes
  if (is.null(ensembl_ids) & is.null(entrez_ids)) {
    common_genes <- intersect(rownames(coef_matrix), rownames(exprs))
    n_obs <- length(common_genes)
  } else {
    common_genes <- intersect(gene_ids, rownames(coef_matrix))
    n_obs <- sum(gene_ids %in% common_genes)
  }
  if (! quiet) {
    message(paste0("Using ", n_obs, "/", nrow(exprs), " observations from the input expression matrix."))
  }

  # Handle provided ensembl IDs many-to-one observations
  if (is.null(ensembl_ids) & is.null(entrez_ids)) {
    common_exprs <- exprs[common_genes,]
  } else if (! any(duplicated(gene_ids))) {
    rownames(exprs) <- gene_ids
    common_exprs <- exprs[common_genes,]
  } else {
    # TODO: Refactor this later
    probe_to_gene_id <- gene_ids
    names(probe_to_gene_id) <- rownames(exprs)

    # Handle one-to-one observations
    single_genes <- names(which(table(gene_ids) == 1))

    single_probes <- names(probe_to_gene_id)[probe_to_gene_id %in% single_genes]
    single_exprs <- exprs[single_probes,]
    rownames(single_exprs) <- probe_to_gene_id[rownames(single_exprs)]

    # Consolidate ensembl IDs that have multiple probes
    multiple_genes <- names(which(table(gene_ids) > 1))
    if (! quiet) {
      message(paste0("Consolidating ", length(multiple_genes), "genes that have multiple observations."))
    }
    consolidated_list <- list()
    for (gid in multiple_genes) {
      illumina_ids <- names(which(probe_to_gene_id == gid))
      # TODO: Currently only using means
      consolidated_list[[gid]] <- colMeans(exprs[illumina_ids,])
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

  # TODO: Return residuals?

  # TODO: Return an S4 object?
  return(list(estimated_time=estimated_time,
              mse=mse))
}

identify_id_type <- function(x) {
  if (all(grepl("^ENSG", x))) {
    return("ensembl")
  } else if (all(grepl("^\\d+$", x))) {
    return("entrez")
  } else {
    stop("Cannot determine gene IDs in expression matrix. Rownames should either be all Ensembl or ENTREZ identifiers, or manually supplied using the ensembl_ids or entrez_ids argument.")
  }
}

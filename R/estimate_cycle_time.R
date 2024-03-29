#' Estimate cycle time
#'
#' Estimate uterine cycle time of endometrium samples using expression data.
#'
#' @param exprs Expression matrix with rows as gene observations and columns as samples
#' @param ensembl_ids Optional vector of Ensembl gene identifiers corresponding to each row of the expression matrix \code{exprs}.
#' @param entrez_ids Optional vector of Entrez gene identifiers corresponding to each row of the expression matrix \code{exprs}.
#' @param handle_multiple_observations Method of handling multiple observations correspond to the same gene. \code{"remove"} will remove all observations of multiple observations, while \code{"mean"}, \code{"median"}, and \code{"max"} will consolidate multiple observations into one using the specified function.
#' @param return_full_residuals Return a list of residual matrices for all samples if \code{TRUE}. The matrices contain residual values for every time point for all genes.
#' @param model Use full model (recommended) or secretory model.
#' @param quiet Run quietly and don't output messages.
#'
#' @return Returns a list with the following items:
#' \describe{
#'   \item{estimated_time}{A vector of estimated cycle times from 0 to 99 for each sample.}
#'   \item{mse}{A matrix of mean squared error values for each sample at each timepoint.}
#'   \item{residuals}{A matrix of residual values after subtracting the cycle stage effect.}
#'   \item{full_residuals}{If \code{return_full_residuals=TRUE}, a list of matrices containing residual values for all timepoints.}
#' }
#' @export
#' @importFrom stats complete.cases
#'
#' @examples
#' data(example_rna_exprs)
#' results <- estimate_cycle_time(exprs=example_rna_exprs)
#' results$estimated_time
#'
#' data(example_array_exprs)
#' data(example_array_entrez_ids)
#' results <- estimate_cycle_time(exprs=example_array_exprs, entrez_ids=example_array_entrez_ids)
#' results$estimated_time

estimate_cycle_time <- function(exprs, ensembl_ids=NULL, entrez_ids=NULL,
                                handle_multiple_observations=c("mean", "median", "max", "remove"),
                                return_full_residuals=FALSE, model=c("full", "secretory"),
                                quiet=FALSE) {

  handle_multiple_observations <- match.arg(handle_multiple_observations)
  model <- match.arg(model)

  # Make sure exprs is a matrix, not a dataframe
  if (! "matrix" %in% class(exprs)) {
    stop("Error: exprs input must be a matrix.")
  }

  # Check there are no NAs in the input matrix
  if (! all(complete.cases(exprs))) {
    stop("Error: exprs input must not contain NA values.")
  }

  # Check input exprs has sample names as colnames
  if (is.null(colnames(exprs))) {
    stop("Error: Missing sample names. Expression matrix must have column names.")
  }

  # Check colnames are unique
  if (any(duplicated(colnames(exprs)))) {
    stop("Error: Column names in the expression matrix are not unique. Sample names must be unique.")
  }

  # Check gene IDs have been provided
  if (is.null(rownames(exprs)) & is.null(ensembl_ids) & is.null(entrez_ids)) {
    stop("Error: Missing gene identifiers. Ensembl or Entrez gene IDs must be set as rownames in the expression matrix or provided using the ensembl_ids/entrez_ids arguments.")
  }

  # Determine what gene identifiers are used
  if (! is.null(ensembl_ids)) {
    gene_id_type <- "ensembl"
    gene_ids <- ensembl_ids
  } else if (! is.null(entrez_ids)) {
    gene_id_type <- "entrez"
    gene_ids <- as.character(entrez_ids)
  } else {
    gene_id_type <- identify_id_type(rownames(exprs))
    gene_ids <- rownames(exprs)
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

  if (model == "full") {
    coef_matrix <- full_coef_matrix
    lp_matrix <- full_lp_matrix
  } else {
    coef_matrix <- sec_coef_matrix
    lp_matrix <- sec_lp_matrix
  }

  if (gene_id_type == "entrez") {
    # Using entrez IDs means using fewer genes since not all ensembl IDs in the endest model have entrez IDs
    coef_matrix <- coef_matrix[names(ensembl_to_entrez),]
    rownames(coef_matrix) <- ensembl_to_entrez
  }

  # Get common genes
  if (is.null(ensembl_ids) & is.null(entrez_ids)) {
    in_common <- intersect(rownames(coef_matrix), rownames(exprs))
    common_genes <- rownames(exprs)[rownames(exprs) %in% in_common]
    n_obs <- length(common_genes)
  } else {
    in_common <- intersect(gene_ids, rownames(coef_matrix))
    common_genes <- gene_ids[gene_ids %in% in_common]
    n_obs <- sum(gene_ids %in% common_genes)
  }
  if (n_obs == 0) {
    stop(paste0("Error: ", n_obs, " / ", nrow(exprs), " observations from the input expression matrix can be used. Check if your gene IDs are correct."))
  }
  if (! quiet) {
    message(paste0("Using ", n_obs, " / ", nrow(exprs), " observations from the input expression matrix."))
  }

  # Handle possible gene IDs many-to-one observations
  if (is.null(ensembl_ids) & is.null(entrez_ids) & ! any(duplicated(common_genes))) {
    # No genes in the rownames are duplicated
    common_exprs <- exprs[common_genes,]
  } else if (! any(duplicated(common_genes))) {
    # None of the provided genes are duplicated
    rownames(exprs) <- gene_ids
    common_exprs <- exprs[common_genes,]
  } else {
    # Change rownames of exprs to 1:nrow(exprs) since we can't always be sure unique rownames are provided
    probe_to_gene_id <- gene_ids
    rownames(exprs) <- seq_len(nrow(exprs))
    names(probe_to_gene_id) <- rownames(exprs)

    # Handle one-to-one observations
    single_genes <- names(which(table(gene_ids) == 1))
    single_probes <- names(probe_to_gene_id)[probe_to_gene_id %in% single_genes]
    single_exprs <- exprs[single_probes,]
    rownames(single_exprs) <- probe_to_gene_id[rownames(single_exprs)]

    # Consolidate ensembl IDs that have multiple probes
    multiple_genes <- names(which(table(gene_ids) > 1))
    if (handle_multiple_observations == "remove") {
      message(paste0("Removing ", length(multiple_genes), " genes that have multiple observations."))
      common_genes <- common_genes[common_genes %in% rownames(single_exprs)]
      common_exprs <- single_exprs[common_genes,]
    } else {
      if (! quiet) {
        message(paste0("Consolidating ", length(multiple_genes), " genes that have multiple observations using the ", handle_multiple_observations, "."))
      }
      multi_exprs <- consolidate_multiple_probes(multiple_genes, probe_to_gene_id, exprs=exprs, method=handle_multiple_observations)

      # Combine single and multiple matrices
      common_exprs  <- rbind(multi_exprs, single_exprs)
      common_exprs <- common_exprs[common_genes,]
    }
  }

  # Process expression data
  normalised_exprs <- preprocessCore::normalize.quantiles.use.target(common_exprs[common_genes,], target=ref)
  colnames(normalised_exprs) <- colnames(common_exprs)
  rownames(normalised_exprs) <- common_genes

  # Get expected expression matrix
  expected_exprs <- coef_matrix %*% lp_matrix

  # Get the mean squared error
  mse <- sapply(1:ncol(normalised_exprs), function(i) {
    colMeans(sweep(x=expected_exprs[common_genes,], MARGIN=1, STATS=normalised_exprs[,i], FUN="-") ^ 2)
  })
  colnames(mse) <- colnames(normalised_exprs)

  # Estimate cycle time using the minimum MSE
  estimated_time <- rownames(mse)[apply(mse, 2, which.min)]
  estimated_time <- as.numeric(sub("^time_|^pod_", "", estimated_time))
  names(estimated_time) <- colnames(mse)

  # Calculate residuals
  exp <- expected_exprs[common_genes,rownames(mse)[apply(mse, 2, which.min)]]
  residuals <- normalised_exprs - exp

  # Get full residuals if specified
  if (return_full_residuals) {
    full_residuals <- lapply(1:ncol(normalised_exprs), function(i) {
      sweep(x=expected_exprs[common_genes,] * -1, MARGIN=1, STATS=normalised_exprs[,i], FUN="+")
    })
    names(full_residuals) <- colnames(normalised_exprs)
  } else {
    full_residuals <- NULL
  }

  # TODO: Return an S4 object instead of a list?
  return(list(estimated_time=estimated_time,
              mse=mse,
              residuals=residuals,
              full_residuals=full_residuals))
}

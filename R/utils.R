identify_id_type <- function(x) {
  if (all(grepl("^ENSG", x))) {
    return("ensembl")
  } else if (all(grepl("^\\d+$", x))) {
    return("entrez")
  } else {
    stop("Cannot determine gene IDs in expression matrix. Rownames should either be all Ensembl or all Entrez identifiers, or manually supplied using the ensembl_ids or entrez_ids argument.")
  }
}

consolidate_multiple_probes <- function(multiple_genes, probe_to_gene_id, exprs, method="mean") {
  consolidated_list <- list()
  if (method == "mean") {
    f <- mean
  } else if (method == "median") {
    f <- stats::median
  } else {
    f <- max
  }
  for (gid in multiple_genes) {
    illumina_ids <- names(which(probe_to_gene_id == gid))
    consolidated_list[[gid]] <- apply(exprs[illumina_ids,], 2, f)
  }
  return(do.call(rbind, consolidated_list))
}


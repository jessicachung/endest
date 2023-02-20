test_that("Gene ID types can be identified", {
  expect_equal(identify_id_type(ensembl_to_entrez), "entrez")
  expect_equal(identify_id_type(1:100), "entrez")
  expect_equal(identify_id_type(names(ensembl_to_entrez)), "ensembl")
  expect_equal(identify_id_type(paste0("ENSG00000000", 1:100)), "ensembl")
})


test_that("Errors with invalid gene IDs", {
  expect_error(identify_id_type(LETTERS),
               "Cannot determine gene IDs")
  expect_error(identify_id_type(c(1:10, names(ensembl_to_entrez)[1:10])),
               "Cannot determine gene IDs")
})


test_that("Consolidate multiple probes", {
  multiple_genes <- c("A", "B", "C")
  m <- matrix(rep(c(1, 2, 1, 2, 1, 4, 9), times=3), nrow=7, ncol=3)
  colnames(m) <- paste0("sample_", 1:3)
  rownames(m) <- 1:7
  probe_to_gene_id <- c("A", "A", "B", "C", "B", "C", "C")
  names(probe_to_gene_id) <- rownames(m)

  res <- consolidate_multiple_probes(multiple_genes=multiple_genes,
                                     exprs=m,
                                     probe_to_gene_id=probe_to_gene_id)
  expect_true(all(res == matrix(rep(c(1.5, 1, 5), times=3), nrow=3)))

  res <- consolidate_multiple_probes(multiple_genes=multiple_genes,
                                     exprs=m,
                                     probe_to_gene_id=probe_to_gene_id,
                                     method="mean")
  expect_true(all(res == matrix(rep(c(1.5, 1, 5), times=3), nrow=3)))

  res <- consolidate_multiple_probes(multiple_genes=multiple_genes,
                                     exprs=m,
                                     probe_to_gene_id=probe_to_gene_id,
                                     method="median")
  expect_true(all(res == matrix(rep(c(1.5, 1, 4), times=3), nrow=3)))

  res <- consolidate_multiple_probes(multiple_genes=multiple_genes,
                                     exprs=m,
                                     probe_to_gene_id=probe_to_gene_id,
                                     method="max")
  expect_true(all(res == matrix(rep(c(2, 1, 9), times=3), nrow=3)))

})

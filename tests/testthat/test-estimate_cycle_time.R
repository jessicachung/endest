test_that("Errors with missing or duplicate sample names in expression matrix column names", {
  expect_error(estimate_cycle_time(exprs=matrix(0, nrow=5, ncol=5)),
               "Missing sample names")

  m <- matrix(0, nrow=5, ncol=5)
  colnames(m) <- paste0("sample", c(1,2,3,4,1))
  expect_error(estimate_cycle_time(exprs=m),
               "Column names in the expression matrix are not unique")
})

test_that("Errors with gene identifiers", {
  m <- matrix(0, nrow=5, ncol=5)
  colnames(m) <- paste0("sample", 1:5)
  expect_error(estimate_cycle_time(exprs=m),
               "Missing gene identifiers")

  expect_error(estimate_cycle_time(exprs=m, ensembl_ids=LETTERS[1:4]),
               "The number of provided gene IDs does not match the number of rows")

  expect_error(estimate_cycle_time(exprs=m, entrez_ids=1:10),
               "The number of provided gene IDs does not match the number of rows")
})

test_that("Estimation for RNA-seq", {
  data("example_rna_exprs")
  res <- estimate_cycle_time(example_rna_exprs)
  expect_equal(res$estimated_time, c("RS_1"=1, "RS_2"=31, "RS_3"=83))
  expect_equal(dim(res$mse), c(100, 3))
  expect_equal(dim(res$residuals), c(2000, 3))
})

test_that("Estimation for RNA-seq", {
  data("example_array_exprs")
  data("example_array_entrez_ids")
  res <- estimate_cycle_time(example_array_exprs, entrez_ids=example_array_entrez_ids)
  expect_equal(res$estimated_time, c("MA_1"=6, "MA_2"=55, "MA_3"=71))
  expect_equal(dim(res$mse), c(100, 3))
  expect_equal(dim(res$residuals), c(4685, 3))
})

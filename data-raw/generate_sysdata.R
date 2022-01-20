# Generate sysdata for endest

library(splines)
library(mgcv)
library(org.Hs.eg.db)  # version 3.14
library(dplyr)

combat_nc <- readRDS("/Users/jess/Work/2021_endo/endo_molecular_model/data/tidy_data/batch_normalised_exprs.rds")
combat_phenotype <- readRDS("/Users/jess/Work/2021_endo/endo_molecular_model/data/tidy_data/combat_phenotype_2021-03-23.rds")

fit_spline_model <- function(y, x, spline_k=6, spline_fx=FALSE, spline_bs="cc",
                             gamma=1, predict_range=seq(0,28,by=0.1),
                             weights=NULL, knots=NULL, return_coefs=FALSE,
                             return_fit=TRUE) {

  fit <- gam(y ~ s(x, bs=spline_bs, k=spline_k, fx=spline_fx),
             knots=knots, weights=weights, gamma=gamma)
  fit_summary <- summary(fit)

  pred <- predict(fit, newdata=data.frame(x=predict_range))
  pred <- data.frame(cycle_time=predict_range, pred=pred)

  results <- list(pred=pred,
                  edf=round(fit_summary$edf, 3),
                  R2=round(fit_summary$r.sq, 4),
                  dev_exp=round(fit_summary$dev.expl, 4),
                  s_table=fit_summary$s.table)
  if (return_coefs) results$coefs <- fit$coefficients
  if (return_fit) results$fit <- fit
  return(results)
}


# ------------------------------------------------------------
# Fit gene models

# Only study 1 samples
pheno <- combat_phenotype %>% dplyr::filter(! is.na(transformed_time))
exprs <- combat_nc[,pheno$sample_id]
cycle_time <- as.numeric(pheno$transformed_time)
genes <- rownames(exprs)
cycle_range <- seq(0, 100, by=1)

# Manually set knots
k_knots <- 30
spline_knots <- list(x=seq(0,100,length=k_knots))

gene_models <- lapply(genes, function(p) {
  fit_spline_model(y=exprs[p,], x=cycle_time, spline_k=k_knots, gamma=4, knots=spline_knots,
                   predict_range=cycle_range, return_coefs=TRUE)
})
names(gene_models) <- genes


# ------------------------------------------------------------
# Get coefficients and linear predictors

# Get one fit object in order to get s(x) values
fit <- fit_spline_model(y=exprs[1,], x=cycle_time, spline_k=k_knots, gamma=4, knots=spline_knots,
                        predict_range=cycle_range, return_fit=TRUE)$fit
lpmatrix <- predict(fit, newdata=data.frame(x=cycle_range), type="lpmatrix")

# Extract coefficients
coef_matrix <- do.call(rbind, lapply(gene_models, function(x) x$coefs))

# Multiply linear predictor matrix by coefs to check predictions are correct
exp_matrix_check <- coef_matrix %*% t(lpmatrix)
expected_exprs <- sapply(genes, function(p) gene_models[[p]]$pred$pred) %>%
  t %>% magrittr::set_colnames(paste0("time_", cycle_range))

# Sanity check all values are identical
stopifnot(all(exp_matrix_check == expected_exprs))

# Check that the first column is the same as the last column (time 0 is the same as time 100)
# stopifnot(all(expected_exprs[,1] == expected_exprs[,101]))
stopifnot(all(round(expected_exprs[,1], 10) == round(expected_exprs[,101], 10)))

# Remove time 100 from lpmatrix
lp_matrix <- t(lpmatrix)[,-101]
colnames(lp_matrix) <- paste0("time_", cycle_range[-101])


# ------------------------------------------------------------
# Get ref for normalisation

# Get one example sample, so we can use it as a reference for quantile normalisation
tmp <- apply(combat_nc, 2, sort)
tmp <- apply(tmp, 1, median)
ref <- tmp


# ------------------------------------------------------------
# Get Entrez to Ensembl gene identifiers

ensembl_to_entrez_list <- as.list(org.Hs.egENSEMBL2EG)
ensembl_ids <- rownames(coef_matrix)

table(ensembl_ids %in% names(ensembl_to_entrez_list))
# FALSE  TRUE
#  2818 17249

ensembl_to_entrez_list <- ensembl_to_entrez_list[ensembl_ids[ensembl_ids %in% names(ensembl_to_entrez_list)]]

table(sapply(ensembl_to_entrez_list, length))
# 1     2     3     4
# 17120   123     5     1

# Only get those with a one-to-one match
ensembl_to_entrez <- sapply(ensembl_to_entrez_list, function(x) {x[1]})
ensembl_to_entrez <- ensembl_to_entrez[sapply(ensembl_to_entrez_list, length) == 1]
length(ensembl_to_entrez)
# [1] 17120
ensembl_to_entrez[1:10]


# ------------------------------------------------------------
# Round coef_matrix to save space

coef_matrix <- round(coef_matrix, 6)

# ------------------------------------------------------------
# Save sysdata

save(coef_matrix, lp_matrix, ref, ensembl_to_entrez, file="R/sysdata.rda", version=2)

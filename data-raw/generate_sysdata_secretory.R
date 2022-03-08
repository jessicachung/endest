# Generate sysdata for endest

library(splines)
library(mgcv)
library(org.Hs.eg.db)  # version 3.14
library(dplyr)
library(stringr)

combat_nc <- readRDS("/Users/jess/Work/2021_endo/endo_molecular_model/data/tidy_data/batch_normalised_exprs.rds")
# combat_phenotype <- readRDS("/Users/jess/Work/2021_endo/endo_molecular_model/data/tidy_data/combat_phenotype_2021-03-23.rds")
sample_info <- readRDS("/Users/jess/Work/2021_endo/endo_molecular_model/data/tidy_data/sample_info.rds")


s1_secretory_pheno <- sample_info %>%
  filter(str_detect(sample_id, "^X")) %>%
  mutate(hr=as.numeric(original_histology_dating_pod_or_helene_rees),
         vo=as.numeric(vanessa_obers_pathology_pod)) %>%
  rowwise() %>%
  mutate(pod=mean(c(hr, vo), na.rm=TRUE)) %>%
  ungroup() %>%
  filter(is.finite(pod)) %>%
  select(sample_id, pod) %>%
  mutate(study="study_1")
s2_secretory_pheno <- sample_info %>%
  filter(! is.na(consensus_pod_urs_only)) %>%
  mutate(pod=as.numeric(consensus_pod_urs_only)) %>%
  select(sample_id, pod) %>%
  mutate(study="study_2")



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

subset_pheno <- rbind(s1_secretory_pheno, s2_secretory_pheno) %>%
  filter(sample_id %in% colnames(combat_nc))
table(subset_pheno$pod)

subset_cycle_time <- as.numeric(subset_pheno$pod)
subset_cycle_range <- seq(1, 14, by=0.1)

exprs <- combat_nc[,subset_pheno$sample_id]
probes <- rownames(combat_nc)
gene_models <- lapply(probes, function(p) {
  fit_spline_model(y=exprs[p,], x=subset_cycle_time,
                   spline_k=3, gamma=2, spline_bs="cr",
                   predict_range=subset_cycle_range, weights=NULL,
                   return_fit=FALSE, return_coefs=TRUE)
})
names(gene_models) <- probes


# ------------------------------------------------------------
# Get coefficients and linear predictors

# Get one fit object in order to get s(x) values
fit <- fit_spline_model(y=exprs[1,], x=subset_cycle_time, spline_k=3, gamma=2, spline_bs="cr",
                        predict_range=subset_cycle_range, weights=NULL, return_fit=TRUE)$fit
lpmatrix <- predict(fit, newdata=data.frame(x=subset_cycle_range), type="lpmatrix")

# Extract coefficients
coef_matrix <- do.call(rbind, lapply(gene_models, function(x) x$coefs))

# Multiply linear predictor matrix by coefs to check predictions are correct
exp_matrix_check <- coef_matrix %*% t(lpmatrix)
expected_exprs <- sapply(probes, function(p) gene_models[[p]]$pred$pred) %>%
  t %>% magrittr::set_colnames(paste0("time_", subset_cycle_range))

# Sanity check all values are identical
stopifnot(all(exp_matrix_check == expected_exprs))

lp_matrix <- t(lpmatrix)
colnames(lp_matrix) <- paste0("pod_", subset_cycle_range)

# ------------------------------------------------------------
# Get ref for normalisation

# Get one example sample, so we can use it as a reference for quantile normalisation
# tmp <- apply(combat_nc, 2, sort)
# tmp <- apply(tmp, 1, median)
# ref <- tmp


# ------------------------------------------------------------
# Get Entrez to Ensembl gene identifiers

# ensembl_to_entrez_list <- as.list(org.Hs.egENSEMBL2EG)
# ensembl_ids <- rownames(coef_matrix)
#
# table(ensembl_ids %in% names(ensembl_to_entrez_list))
# # FALSE  TRUE
# #  2818 17249
#
# ensembl_to_entrez_list <- ensembl_to_entrez_list[ensembl_ids[ensembl_ids %in% names(ensembl_to_entrez_list)]]
#
# table(sapply(ensembl_to_entrez_list, length))
# # 1     2     3     4
# # 17120   123     5     1
#
# # Only get those with a one-to-one match
# ensembl_to_entrez <- sapply(ensembl_to_entrez_list, function(x) {x[1]})
# ensembl_to_entrez <- ensembl_to_entrez[sapply(ensembl_to_entrez_list, length) == 1]
# length(ensembl_to_entrez)
# # [1] 17120
# ensembl_to_entrez[1:10]


# ------------------------------------------------------------
# Round coef_matrix to save space

coef_matrix <- round(coef_matrix, 6)

# ------------------------------------------------------------
# Save sysdata

sec_coef_matrix <- coef_matrix
sec_lp_matrix <- lp_matrix

save(sec_coef_matrix, sec_lp_matrix, file="tmp/sec_sysdata.rda", version=2)

# # Check size
# tools::checkRdaFiles("R/sysdata.rda")
# #                  size ASCII compress version
# # R/sysdata.rda 3644508 FALSE     gzip       2
#
# # Re-save with xz compression to reduce size
# tools::resaveRdaFiles("R/sysdata.rda", compress="xz", version=2)
#
# # Check size
# tools::checkRdaFiles("R/sysdata.rda")
# #                  size ASCII compress version
# # R/sysdata.rda 2003824 FALSE       xz       2

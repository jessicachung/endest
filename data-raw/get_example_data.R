# Get example data

devtools::load_all(".")
library(dplyr)

# RNA-seq expression data
combat_nc <- readRDS("/Users/jess/Work/2021_endo/endo_molecular_model/data/tidy_data/batch_normalised_exprs.rds")
combat_phenotype <- readRDS("/Users/jess/Work/2021_endo/endo_molecular_model/data/tidy_data/combat_phenotype_2021-03-23.rds")

# Microarray expression data
combat_array <- readRDS("/Users/jess/Work/2021_endo/endo_molecular_model/data/array_data/combat_exprs.rds")
probe_info <- readRDS("/Users/jess/Work/2021_endo/endo_molecular_model/data/array_data/illumina_v4_annotation.rds")

# Also load the coef_matrix
load("R/sysdata.rda")

# ------------------------------------------------------------
# Get a subset of RNA-seq expression data

# Get an example of a menstrual sample, a proliferative sample, and a secretory sample
set.seed(2021)
rna_samples <- c(
  combat_phenotype %>% filter(transformed_time %>% between(0,8)) %>% pull(sample_id) %>% sample(1),
  combat_phenotype %>% filter(transformed_time %>% between(8,58)) %>% pull(sample_id) %>% sample(1),
  combat_phenotype %>% filter(transformed_time %>% between(58,100)) %>% pull(sample_id) %>% sample(1)
)

rna_samples
# [1] "X210494" "X210164" "X210686"

# Check
res <- estimate_cycle_time(combat_nc[,rna_samples])
res$estimated_time
plot_mse(res, rna_samples[1])
plot_mse(res, rna_samples[2])
plot_mse(res, rna_samples[3])

# Reduce data size by subsampling genes
# Bias sampling by weighting curvy gene models higher
rel_prob <- rowSums(abs(coef_matrix))
subset_genes <- sample(rownames(coef_matrix), size=2000, prob=rel_prob)

res <- estimate_cycle_time(combat_nc[subset_genes,rna_samples])
res$estimated_time
plot_mse(res, rna_samples[1])
plot_mse(res, rna_samples[2])
plot_mse(res, rna_samples[3])

example_rna_exprs <- combat_nc[subset_genes,rna_samples]
colnames(example_rna_exprs) <- paste0("RS_", 1:3)


# ------------------------------------------------------------
# Get a subset of microarray expression data

# Get an example of a menstrual sample, a proliferative sample, and a secretory sample
# res <- estimate_cycle_time(combat_array[,1:100], entrez_ids=probe_info$EntrezReannotated)
# sort(res$estimated_time)
array_samples <- c("X210061", "X210210", "X210102")

# Check
res <- estimate_cycle_time(combat_array[,array_samples],
                           entrez_ids=probe_info$EntrezReannotated)
res$estimated_time
plot_mse(res, array_samples[1])
plot_mse(res, array_samples[2])
plot_mse(res, array_samples[3])

# Reduce data size by subsampling genes
# Bias sampling by weighting curvy gene models higher
rel_prob <- rowSums(abs(coef_matrix))
names(rel_prob) <- ensembl_to_entrez[names(rel_prob)]
rel_prob <- rel_prob[! is.na(names(rel_prob))]

probe_info <- probe_info %>%
  mutate(x=rel_prob[EntrezReannotated],
         x=ifelse(is.na(x), 1, x))
subset_probes <- sample(probe_info$IlluminaID, size=5000, prob=probe_info$x)


# Check
m <- match(subset_probes, probe_info$IlluminaID)
entrez <- probe_info$EntrezReannotated[m]
res <- estimate_cycle_time(combat_array[subset_probes,array_samples],
                           entrez_ids=entrez)
res$estimated_time
plot_mse(res, array_samples[1])
plot_mse(res, array_samples[2])
plot_mse(res, array_samples[3])

example_array_exprs <- combat_array[subset_probes,array_samples]
colnames(example_array_exprs) <- paste0("MA_", 1:3)
example_array_entrez_ids <- entrez
names(example_array_entrez_ids) <- rownames(example_array_exprs)


# ------------------------------------------------------------
# Save

usethis::use_data(example_rna_exprs, overwrite=TRUE)
usethis::use_data(example_array_exprs, overwrite=TRUE)
usethis::use_data(example_array_entrez_ids, overwrite=TRUE)

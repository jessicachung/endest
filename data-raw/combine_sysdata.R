load("tmp/full_sysdata.rda")
load("tmp/sec_sysdata.rda")

# Rename objects
full_coef_matrix <- coef_matrix
full_lp_matrix <- lp_matrix

# head(full_lp_matrix)
# head(full_coef_matrix)
# head(sec_lp_matrix)
# head(sec_coef_matrix)

save(full_lp_matrix,
     full_coef_matrix,
     sec_lp_matrix,
     sec_coef_matrix,
     ensembl_to_entrez,
     ref,
     file="R/sysdata.rda", version=2, compress="xz"
)

# 01_build_analysis_dataset/03_preprocess_metabolomics.R

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

if (!file.exists(files$analysis_cohort)) {
  stop("Missing analysis cohort file: ", files$analysis_cohort)
}

analysis_cohort <- readRDS(files$analysis_cohort)
metab_pd <- as.data.frame(analysis_cohort$metab_pd)
pd_ids <- as.character(analysis_cohort$sample_ids)

if (anyDuplicated(pd_ids)) {
  stop("Duplicated sample IDs in analysis cohort.")
}
if (!identical(rownames(metab_pd), pd_ids)) {
  stop("Metabolomics row order does not match sample_ids.")
}
if (anyNA(metab_pd)) {
  stop("Metabolomics matrix contains missing values before filtering.")
}

zero_var_keep <- apply(metab_pd, 2, var, na.rm = TRUE) > 0
metab_pd <- metab_pd[, zero_var_keep, drop = FALSE]

if (ncol(metab_pd) == 0) {
  stop("No metabolomics features remained after zero-variance filtering.")
}

metab_scaled <- scale(metab_pd)
metab_scaled <- as.matrix(metab_scaled)
rownames(metab_scaled) <- pd_ids
colnames(metab_scaled) <- colnames(metab_pd)

if (!identical(rownames(metab_scaled), pd_ids)) {
  stop("Scaled metabolomics row order does not match sample_ids.")
}
if (anyNA(metab_scaled)) {
  stop("Scaled metabolomics matrix contains NA values.")
}

metabolomics_processed <- list(
  metab_pd = metab_pd,
  metab_scaled = metab_scaled,
  sample_ids = pd_ids
)

saveRDS(metabolomics_processed, files$metabolomics_processed)

write.csv(
  data.frame(
    metric = c("n_final_samples", "n_features_after_zero_var"),
    value = c(nrow(metab_scaled), ncol(metab_scaled))
  ),
  file.path(paths$analysis_dataset, "metabolomics_preprocess_summary.csv"),
  row.names = FALSE
)
umap_input_clean <- data.frame(
  sample_id = pd_ids,
  metab_pd,
  check.names = FALSE,
  row.names = NULL
)
write.csv(
  umap_input_clean,
  file.path(paths$analysis_dataset, "umap_input_clean.csv"),
  row.names = FALSE
)

message("Metabolomics preprocessing complete.")
message("Patients: ", nrow(metab_scaled))
message("Features: ", ncol(metab_scaled))
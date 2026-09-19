# 01_build_analysis_dataset/03_preprocess_metabolomics.R

# This script preprocesses the metabolomics data from the analysis cohort.
# It checks sample IDs and missing values, removes zero-variance features,
# standardises the remaining features, and saves the processed data together
# with a preprocessing summary and a clean CSV input for UMAP.
#
# Before running:
# - Set `project_root` to the local project directory.
# - `files$analysis_cohort` must point to an existing analysis cohort created
#   by `01_build_clean_cohort.R`.
#
# Outputs:
# - `metabolomics_processed.rds`, containing the filtered and scaled
#   metabolomics data and sample IDs.
# - `metabolomics_preprocess_summary.csv`
# - `umap_input_clean.csv`

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

# Check that the analysis cohort created by the previous preprocessing step exists.
if (!file.exists(files$analysis_cohort)) {
  stop("Missing analysis cohort file: ", files$analysis_cohort)
}

analysis_cohort <- readRDS(files$analysis_cohort)
metab_pd <- as.data.frame(analysis_cohort$metab_pd)
pd_ids <- as.character(analysis_cohort$sample_ids)

# Check that sample IDs are unique and that the metabolomics rows are in the
# same order as the sample IDs stored in the analysis cohort.
if (anyDuplicated(pd_ids)) {
  stop("Duplicated sample IDs in analysis cohort.")
}
if (!identical(rownames(metab_pd), pd_ids)) {
  stop("Metabolomics row order does not match sample_ids.")
}

# Require the input metabolomics matrix to contain no missing values before filtering.
if (anyNA(metab_pd)) {
  stop("Metabolomics matrix contains missing values before filtering.")
}

# Remove features with zero variance across the analysis cohort.
zero_var_keep <- apply(metab_pd, 2, var, na.rm = TRUE) > 0
metab_pd <- metab_pd[, zero_var_keep, drop = FALSE]

if (ncol(metab_pd) == 0) {
  stop("No metabolomics features remained after zero-variance filtering.")
}

# Standardise the retained metabolomics features.
metab_scaled <- scale(metab_pd)
metab_scaled <- as.matrix(metab_scaled)
rownames(metab_scaled) <- pd_ids
colnames(metab_scaled) <- colnames(metab_pd)

# Check that scaling has preserved the expected sample order and introduced
# no missing values.
if (!identical(rownames(metab_scaled), pd_ids)) {
  stop("Scaled metabolomics row order does not match sample_ids.")
}
if (anyNA(metab_scaled)) {
  stop("Scaled metabolomics matrix contains NA values.")
}

# Store both the filtered unscaled data and the scaled data together with sample IDs.
metabolomics_processed <- list(
  metab_pd = metab_pd,
  metab_scaled = metab_scaled,
  sample_ids = pd_ids
)

# Save the processed metabolomics data.
saveRDS(metabolomics_processed, files$metabolomics_processed)

# Save a summary of the number of samples and features remaining after filtering.
write.csv(
  data.frame(
    metric = c("n_final_samples", "n_features_after_zero_var"),
    value = c(nrow(metab_scaled), ncol(metab_scaled))
  ),
  file.path(paths$analysis_dataset, "metabolomics_preprocess_summary.csv"),
  row.names = FALSE
)

# Create a CSV containing the filtered, unscaled metabolomics data with sample IDs
# as the first column for use as a clean input to the downstream UMAP script.
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

# 01_build_combined_clinical_affinity.R

# This script builds one combined clinical affinity matrix from UPDRS, MoCA,
# and LEDD for the representation-comparison sensitivity analysis.
#
# Fixed settings:
# - SNF K = `snf_params$K`
# - SNF sigma = `snf_params$sigma`
#
# Before running:
# - Set `project_root` to the local project directory.
# - `clinical_processed.rds` must already exist and contain UPDRS_adj,
#   MOCA_adj, LEDD_adj, and matching sample IDs.
#
# Outputs are saved under
# `05_sensitivity_analysis/outputs/representation_comparison/`, including
# the combined clinical matrices, affinity matrix, QC summary, and affinity
# distribution plot.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))
source(file.path(project_root, "05_sensitivity_analysis", "shared", "snf_helpers.R"))

out_dir <- file.path(project_root, "05_sensitivity_analysis", "outputs", "representation_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clinical_processed_file <- file.path(project_root, "01_build_analysis_dataset", "outputs", "clinical_processed.rds")
if (!file.exists(clinical_processed_file)) {
  stop("Missing clinical processed file: ", clinical_processed_file)
}

clinical_obj <- readRDS(clinical_processed_file)
pd_model <- as.data.frame(clinical_obj$pd_model)

required_cols <- c("UPDRS_adj", "MOCA_adj", "LEDD_adj")
missing_cols <- setdiff(required_cols, names(pd_model))
if (length(missing_cols) > 0) {
  stop("Missing required columns in clinical_processed.rds: ", paste(missing_cols, collapse = ", "))
}

sample_ids <- as.character(clinical_obj$sample_ids)
if (length(sample_ids) == 0L) {
  stop("clinical_processed.rds contains no sample IDs.")
}
if (is.null(rownames(pd_model))) {
  rownames(pd_model) <- sample_ids
}
if (!identical(rownames(pd_model), sample_ids)) {
  stop("Clinical processed row order does not match sample_ids.")
}

combined_clinical <- as.matrix(pd_model[, required_cols, drop = FALSE])
rownames(combined_clinical) <- sample_ids

if (anyNA(combined_clinical)) {
  stop("Combined clinical matrix contains NA values.")
}

combined_scaled <- scale(combined_clinical)
combined_scaled <- as.matrix(combined_scaled)
rownames(combined_scaled) <- sample_ids
colnames(combined_scaled) <- required_cols

if (anyNA(combined_scaled)) {
  stop("Scaled combined clinical matrix contains NA values.")
}

dist_combined <- SNFtool::dist2(combined_scaled, combined_scaled)
combined_clinical_affinity_mat <- SNFtool::affinityMatrix(
  dist_combined,
  K = snf_params$K,
  sigma = snf_params$sigma
)
rownames(combined_clinical_affinity_mat) <- sample_ids
colnames(combined_clinical_affinity_mat) <- sample_ids

check_affinity_matrix(combined_clinical_affinity_mat, "combined_clinical_affinity_mat")

saveRDS(combined_clinical, file.path(out_dir, "combined_clinical_raw.rds"))
saveRDS(combined_scaled, file.path(out_dir, "combined_clinical_scaled.rds"))
saveRDS(combined_clinical_affinity_mat, file.path(out_dir, "combined_clinical_affinity_mat.rds"))

offdiag <- combined_clinical_affinity_mat[upper.tri(combined_clinical_affinity_mat)]
write.csv(
  data.frame(
    metric = c("n_samples", "n_features", "min_affinity", "median_affinity", "mean_affinity", "max_affinity"),
    value = c(
      nrow(combined_clinical_affinity_mat),
      ncol(combined_scaled),
      min(offdiag, na.rm = TRUE),
      median(offdiag, na.rm = TRUE),
      mean(offdiag, na.rm = TRUE),
      max(offdiag, na.rm = TRUE)
    )
  ),
  file.path(out_dir, "combined_clinical_affinity_qc.csv"),
  row.names = FALSE
)

png(file.path(out_dir, "combined_clinical_affinity_histogram.png"), width = 2000, height = 1500, res = 300)
hist(offdiag, breaks = 50, main = "Combined clinical affinity distribution", xlab = "Affinity")
dev.off()

message("Combined clinical affinity complete.")
message("Saved outputs to: ", out_dir)

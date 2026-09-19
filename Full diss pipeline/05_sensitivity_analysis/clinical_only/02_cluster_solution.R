# 02_cluster_solution.R
# Uses the same eigengap + spectral clustering logic as the shared helpers
# 
# This script derives the cluster solution for the clinical-only sensitivity
# analysis from the fused three-network SNF matrix. It selects the number of
# clusters using the eigengap and then performs spectral clustering.
#
# Fixed settings:
# - Eigengap scan considers the first 10 candidate cluster numbers.
#
# Before running:
# - Set `project_root` to the local project directory.
# - `clinical_only_snf_results.rds` must already have been created by
#   `01_run_snf.R`.
#
# Outputs are saved under `05_sensitivity_analysis/outputs/clinical_only/`.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))
source(file.path(project_root, "05_sensitivity_analysis", "shared", "cluster_helpers.R"))

out_dir <- file.path(project_root, "05_sensitivity_analysis", "outputs", "clinical_only")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

snf_obj_file <- file.path(out_dir, "clinical_only_snf_results.rds")
if (!file.exists(snf_obj_file)) {
  stop("Missing SNF results file: ", snf_obj_file)
}

snf_obj <- readRDS(snf_obj_file)

if (!all(c("fused_network", "sample_ids") %in% names(snf_obj))) {
  stop("clinical_only_snf_results.rds must contain 'fused_network' and 'sample_ids'.")
}

W <- as.matrix(snf_obj$fused_network)
sample_ids <- as.character(snf_obj$sample_ids)

if (is.null(rownames(W)) || is.null(colnames(W))) {
  stop("Fused network must have row and column names.")
}
if (!identical(rownames(W), sample_ids)) {
  stop("Sample IDs do not match fused network row order.")
}
if (!isTRUE(all.equal(W, t(W)))) {
  stop("Fused network is not symmetric.")
}
if (anyNA(W)) {
  stop("Fused network contains NA values.")
}

scan <- compute_eigengap_scan(W, max_k = 10L)

chosen_k <- scan$best_k
if (is.na(chosen_k) || chosen_k < 2L) {
  stop("Chosen k is invalid.")
}

cluster_labels <- run_spectral_clustering(W, chosen_k)

save_cluster_outputs(
  out_dir = out_dir,
  prefix = "clinical_only",
  sample_ids = sample_ids,
  fused_network = W,
  scan = scan,
  chosen_k = chosen_k,
  cluster_labels = cluster_labels,
  model_name = "clinical_only_3network"
)

message("Clinical-only clustering complete.")
message("Chosen k: ", chosen_k)
message("Saved outputs to: ", out_dir)

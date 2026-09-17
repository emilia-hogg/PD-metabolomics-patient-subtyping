# Top-20 combined sensitivity analysis
# Cluster the fused 2-network solution using eigengap + spectral clustering.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))
source(file.path(project_root, "05_sensitivity_analysis", "shared", "top20_snf_helpers.R"))

branch_dir <- file.path(paths$sensitivity, "top20_combined")
out_dir <- file.path(branch_dir, "cluster_solution")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

snf_results_file <- file.path(branch_dir, "top20_combined_2network_snf_results.rds")
if (!file.exists(snf_results_file)) {
  stop("Missing SNF results file: ", snf_results_file)
}

snf_obj <- readRDS(snf_results_file)
if (!all(c("fused_network", "sample_ids") %in% names(snf_obj))) {
  stop("top20_combined_2network_snf_results.rds must contain fused_network and sample_ids.")
}

fused_network <- as.matrix(snf_obj$fused_network)
sample_ids <- as.character(snf_obj$sample_ids)

scan <- run_eigengap_scan(fused_network, max_k = 10L)
chosen_k <- scan$best_k
if (is.na(chosen_k) || chosen_k < 2L) {
  stop("Chosen k is invalid.")
}

cluster_labels <- run_spectral_clustering(fused_network, chosen_k)
second_cluster_labels <- NULL
if (!is.na(scan$second_best_k) && scan$second_best_k >= 2L) {
  second_cluster_labels <- run_spectral_clustering(fused_network, scan$second_best_k)
}

cluster_out <- save_cluster_outputs(
  out_dir = out_dir,
  prefix = "top20_combined_2network",
  sample_ids = sample_ids,
  fused_network = fused_network,
  scan = scan,
  chosen_k = chosen_k,
  cluster_labels = cluster_labels,
  model_name = "top20_combined_2network"
)

if (!is.null(second_cluster_labels)) {
  second_df <- data.frame(sample_id = sample_ids, cluster = as.integer(second_cluster_labels), stringsAsFactors = FALSE)
  write.csv(second_df, file.path(out_dir, "second_best_k_cluster_labels.csv"), row.names = FALSE)
  saveRDS(second_df, file.path(out_dir, "second_best_k_cluster_labels.rds"))
}

message("Top-20 combined clustering complete.")
message("Chosen k: ", chosen_k)
message("Saved outputs to: ", out_dir)

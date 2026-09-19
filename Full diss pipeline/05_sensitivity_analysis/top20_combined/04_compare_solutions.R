# This script compares the top-20 combined two-network SNF solution with the
# primary four-network SNF solution. It compares the fused networks and cluster
# assignments between the two solutions.
#
# Before running:
# - Set `project_root` to the local project directory.
# - The primary four-network SNF and cluster solution must already exist.
# - The top-20 combined SNF and cluster solution must already have been
#   generated.
#
# Outputs are saved under `05_sensitivity_analysis/outputs/top20_combined/comparison/`,
# including network similarity, cluster agreement, matched cluster labels,
# and the branch SNF results.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))
source(file.path(project_root, "05_sensitivity_analysis", "shared", "top20_snf_helpers.R"))

main_snf_dir <- file.path(paths$snf, "main_four_network")
main_cluster_dir <- file.path(main_snf_dir, "cluster_solution")
branch_dir <- file.path(paths$sensitivity, "top20_combined")
branch_cluster_dir <- file.path(branch_dir, "cluster_solution")
out_dir <- file.path(branch_dir, "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

main_fused_file <- file.path(main_snf_dir, "primary_four_network_fused_network.rds")
main_best_labels_file <- file.path(main_cluster_dir, "best_k_cluster_labels.rds")
branch_results_file <- file.path(branch_dir, "top20_combined_2network_snf_results.rds")
branch_best_labels_file <- file.path(branch_cluster_dir, "best_k_cluster_labels.rds")

required_files <- c(main_fused_file, main_best_labels_file, branch_results_file, branch_best_labels_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required file(s): ", paste(missing_files, collapse = "; "))
}

main_fused <- readRDS(main_fused_file)
main_labels <- readRDS(main_best_labels_file)
branch_obj <- readRDS(branch_results_file)
branch_labels <- readRDS(branch_best_labels_file)

if (!all(c("sample_id", "cluster") %in% names(main_labels))) {
  stop("Main best_k_cluster_labels.rds must contain sample_id and cluster.")
}
if (!all(c("sample_id", "cluster") %in% names(branch_labels))) {
  stop("Branch best_k_cluster_labels.rds must contain sample_id and cluster.")
}

main_labels$sample_id <- as.character(main_labels$sample_id)
main_labels$cluster <- as.integer(main_labels$cluster)
branch_labels$sample_id <- as.character(branch_labels$sample_id)
branch_labels$cluster <- as.integer(branch_labels$cluster)

branch_fused <- as.matrix(branch_obj$fused_network)
if (is.null(rownames(branch_fused)) || is.null(colnames(branch_fused))) {
  branch_ids <- as.character(branch_obj$sample_ids)
  rownames(branch_fused) <- branch_ids
  colnames(branch_fused) <- branch_ids
}

network_summary <- compare_networks_upper_triangle(main_fused, branch_fused)

main_cluster_vec <- setNames(main_labels$cluster, main_labels$sample_id)
branch_cluster_vec <- setNames(branch_labels$cluster, branch_labels$sample_id)
shared_ids <- intersect(names(main_cluster_vec), names(branch_cluster_vec))
if (length(shared_ids) < 2L) {
  stop("Not enough shared samples between main and branch cluster labels.")
}

main_cluster_vec <- main_cluster_vec[shared_ids]
branch_cluster_vec <- branch_cluster_vec[shared_ids]

merged_labels <- data.frame(
  sample_id = shared_ids,
  main_cluster = as.integer(main_cluster_vec),
  branch_cluster = as.integer(branch_cluster_vec),
  stringsAsFactors = FALSE
)

summary_obj <- list(
  pairwise_agreement = pairwise_co_clustering_agreement(merged_labels$main_cluster, merged_labels$branch_cluster),
  contingency = table(merged_labels$main_cluster, merged_labels$branch_cluster)
)

prefix <- "top20_combined_2network_vs_main_four_network"
write.csv(network_summary, file.path(out_dir, paste0(prefix, "_network_similarity.csv")), row.names = FALSE)
saveRDS(network_summary, file.path(out_dir, paste0(prefix, "_network_similarity.rds")))

save_representation_comparison_outputs(
  out_dir = out_dir,
  prefix = prefix,
  merged_labels = merged_labels,
  summary_obj = summary_obj
)

saveRDS(branch_obj, file.path(out_dir, paste0(prefix, "_branch_snf_object.rds")))

message("Top-20 combined comparison complete.")
message("Saved outputs to: ", out_dir)

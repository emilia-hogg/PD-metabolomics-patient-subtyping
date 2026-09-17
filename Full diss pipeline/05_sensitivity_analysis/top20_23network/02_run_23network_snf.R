# Top-20 sensitivity analysis
# 23-network SNF: 3 clinical networks + 20 separate metabolite networks

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))
source(file.path(project_root, "05_sensitivity_analysis", "shared", "top20_snf_helpers.R"))

out_dir <- file.path(paths$sensitivity, "top20_23network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

affinity_list_file <- file.path(out_dir, "top20_23network_affinity_list.rds")
if (!file.exists(affinity_list_file)) {
  stop("Missing affinity list file: ", affinity_list_file)
}

network_list <- readRDS(affinity_list_file)
if (!is.list(network_list) || length(network_list) != 23L) {
  stop("Expected 23 networks, but found: ", length(network_list))
}

network_names <- names(network_list)
if (is.null(network_names) || any(network_names == "")) {
  network_names <- paste0("network_", seq_along(network_list))
  names(network_list) <- network_names
}

sample_ids <- check_network_list(network_list)

K <- snf_params$K
T <- snf_params$t

fused_network <- run_snf_from_affinities(network_list = network_list, K = K, T = T)
qc_table <- summarise_fused_network(
  fused_network = fused_network,
  K = K,
  T = T,
  network_names = network_names
)

saveRDS(
  list(
    model_name = "top20_23network",
    sample_ids = sample_ids,
    K = K,
    T = T,
    network_names = network_names,
    inputs = network_list,
    fused_network = fused_network,
    qc = qc_table
  ),
  file.path(out_dir, "top20_23network_snf_results.rds")
)

saveRDS(fused_network, file.path(out_dir, "top20_23network_fused_network.rds"))
write.csv(qc_table, file.path(out_dir, "top20_23network_qc.csv"), row.names = FALSE)

message("Top-20 23-network SNF complete.")
message("Saved outputs to: ", out_dir)
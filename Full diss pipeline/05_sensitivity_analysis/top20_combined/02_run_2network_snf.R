# This script runs a two-network SNF using the combined clinical affinity
# matrix and the combined top-20 metabolite affinity matrix.
#
# Fixed settings:
# - K = `snf_params$K`
# - T = `snf_params$t`
#
# Before running:
# - Set `project_root` to the local project directory.
# - The combined clinical and top-20 metabolite affinity matrices must already
#   have been created by `01_build_top20_network.R`.
#
# Outputs are saved under `05_sensitivity_analysis/outputs/top20_combined/`,
# including the fused network, QC summary, and SNF results.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))
source(file.path(project_root, "05_sensitivity_analysis", "shared", "top20_snf_helpers.R"))

out_dir <- file.path(paths$sensitivity, "top20_combined")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clinical_aff_file <- file.path(out_dir, "combined_clinical_affinity_mat.rds")
metab_aff_file <- file.path(out_dir, "top20_metabolites_affinity_mat.rds")

required_files <- c(clinical_aff_file, metab_aff_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required affinity file(s): ", paste(missing_files, collapse = "; "))
}

combined_clinical_affinity_mat <- readRDS(clinical_aff_file)
top20_metabolites_affinity_mat <- readRDS(metab_aff_file)

network_list <- list(
  combined_clinical = combined_clinical_affinity_mat,
  top20_metabolites = top20_metabolites_affinity_mat
)
network_names <- names(network_list)

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
    model_name = "top20_combined_2network",
    sample_ids = sample_ids,
    K = K,
    T = T,
    network_names = network_names,
    inputs = network_list,
    fused_network = fused_network,
    qc = qc_table
  ),
  file.path(out_dir, "top20_combined_2network_snf_results.rds")
)

saveRDS(fused_network, file.path(out_dir, "top20_combined_2network_fused_network.rds"))
write.csv(qc_table, file.path(out_dir, "top20_combined_2network_qc.csv"), row.names = FALSE)

message("Top-20 combined 2-network SNF complete.")
message("Saved outputs to: ", out_dir)

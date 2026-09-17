# =========================================================
# 02_run_2network_snf.R
# Representation comparison branch
# 2-network SNF: Combined clinical + Metabolomics
# =========================================================

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))
source(file.path(project_root, "05_sensitivity_analysis", "shared", "snf_helpers.R"))

out_dir <- file.path(project_root, "05_sensitivity_analysis", "outputs", "representation_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

combined_file <- file.path(out_dir, "combined_clinical_affinity_mat.rds")
metabolomics_file <- file.path(project_root, "03_network_construction", "outputs", "metabolomics_network", "metab_affinity_mat.rds")

required_files <- c(combined_file, metabolomics_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required affinity file(s): ", paste(missing_files, collapse = "; "))
}

combined_clinical_affinity_mat <- readRDS(combined_file)
metab_affinity_mat <- readRDS(metabolomics_file)

network_list <- list(
  combined_clinical = combined_clinical_affinity_mat,
  metabolomics = metab_affinity_mat
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

save_fused_network_outputs(
  fused_network = fused_network,
  sample_ids = sample_ids,
  out_dir = out_dir,
  prefix = "representation_2network",
  qc_table = qc_table
)

saveRDS(
  list(
    model_name = "representation_2network",
    sample_ids = sample_ids,
    K = K,
    T = T,
    network_names = network_names,
    inputs = network_list,
    fused_network = fused_network,
    qc = qc_table
  ),
  file.path(out_dir, "representation_2network_snf_results.rds")
)

message("2-network SNF complete.")
message("Saved outputs to: ", out_dir)
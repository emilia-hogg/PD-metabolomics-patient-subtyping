# 01_run_snf.R
# This script runs the clinical-only sensitivity analysis using a three-network
# SNF model consisting of UPDRS, MoCA, and LEDD.
#
# Fixed settings:
# - K = `snf_params$K`
# - T = `snf_params$t`
#
# Before running:
# - Set `project_root` to the local project directory.
# - The UPDRS, MoCA, and LEDD affinity matrices must already exist in the
#   clinical network output directory.
#
# Outputs are saved under `05_sensitivity_analysis/outputs/clinical_only/`,
# including the fused network, QC summary, and SNF results.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))
source(file.path(project_root, "05_sensitivity_analysis", "shared", "snf_helpers.R"))

out_dir <- file.path(project_root, "05_sensitivity_analysis", "outputs", "clinical_only")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clinical_dir <- file.path(project_root, "03_network_construction", "outputs", "clinical_network")

required_files <- c(
  file.path(clinical_dir, "updrs_affinity_mat.rds"),
  file.path(clinical_dir, "moca_affinity_mat.rds"),
  file.path(clinical_dir, "ledd_affinity_mat.rds")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required clinical affinity file(s): ", paste(missing_files, collapse = "; "))
}

updrs_affinity_mat <- readRDS(file.path(clinical_dir, "updrs_affinity_mat.rds"))
moca_affinity_mat  <- readRDS(file.path(clinical_dir, "moca_affinity_mat.rds"))
ledd_affinity_mat  <- readRDS(file.path(clinical_dir, "ledd_affinity_mat.rds"))

network_list <- list(
  updrs = updrs_affinity_mat,
  moca = moca_affinity_mat,
  ledd = ledd_affinity_mat
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
  prefix = "clinical_only",
  qc_table = qc_table
)

saveRDS(
  list(
    model_name = "clinical_only_3network",
    sample_ids = sample_ids,
    K = K,
    T = T,
    network_names = network_names,
    inputs = network_list,
    fused_network = fused_network,
    qc = qc_table
  ),
  file.path(out_dir, "clinical_only_snf_results.rds")
)

message("Clinical-only SNF complete.")
message("Saved outputs to: ", out_dir)

# 00_config/config.R

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

paths <- list(
  root = project_root,
  config = file.path(project_root, "00_config"),
  analysis_dataset = file.path(project_root, "01_build_analysis_dataset", "outputs"),
  exploratory = file.path(project_root, "02_exploratory_analysis", "outputs"),
  clinical_networks = file.path(project_root, "03_network_construction", "outputs"),
  snf = file.path(project_root, "04_snf", "outputs"),
  sensitivity = file.path(project_root, "05_sensitivity_analysis", "outputs"),
  posthoc = file.path(project_root, "06_posthoc_analysis", "outputs"),
  dissertation_figures = file.path(project_root, "dissertation_figures"),
  dissertation_tables = file.path(project_root, "dissertation_tables")
)

files <- list(
  workspace_rdata = file.path(project_root, "restarted_pca_kmeans_working_session.RData"),
  clinical_rds = "/mnt/sde/buddhi/PROBAND_shared/Full_Patient_sample_information.rds",
  analysis_cohort = file.path(paths$analysis_dataset, "analysis_cohort.rds"),
  patient_metadata = file.path(paths$analysis_dataset, "patient_metadata.rds"),
  clinical_processed = file.path(paths$analysis_dataset, "clinical_processed.rds"),
  metabolomics_processed = file.path(paths$analysis_dataset, "metabolomics_processed.rds"),
  sample_ids = file.path(paths$analysis_dataset, "sample_ids.rds")
)

snf_params <- list(
  K = 30,
  sigma = 0.8,
  t = 20
)

analysis_flags <- list(
  run_duration_sensitivity = TRUE,
  run_no_ledd_sensitivity = TRUE,
  run_clinical_only_sensitivity = TRUE,
  run_representation_comparison = TRUE,
  run_metabolite_characterisation = TRUE,
  run_threshold_scans = TRUE
)

dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
invisible(lapply(paths[-1], dir.create, recursive = TRUE, showWarnings = FALSE))
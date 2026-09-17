# Top-20 sensitivity analysis (combined clinical + combined metabolite panel)
# Build one combined clinical affinity and one combined top-20 metabolite affinity.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))
source(file.path(project_root, "05_sensitivity_analysis", "shared", "top20_snf_helpers.R"))

out_dir <- file.path(paths$sensitivity, "top20_combined")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clinical_file <- files$clinical_processed
analysis_cohort_file <- files$analysis_cohort

if (!file.exists(clinical_file)) stop("Missing clinical processed file: ", clinical_file)
if (!file.exists(analysis_cohort_file)) stop("Missing analysis cohort file: ", analysis_cohort_file)

clinical_obj <- readRDS(clinical_file)
cohort_obj <- readRDS(analysis_cohort_file)

if (!all(c("pd_model", "sample_ids") %in% names(clinical_obj))) {
  stop("clinical_processed.rds must contain pd_model and sample_ids.")
}
if (!all(c("metab_pd", "sample_ids") %in% names(cohort_obj))) {
  stop("analysis_cohort.rds must contain metab_pd and sample_ids.")
}

pd_model <- as.data.frame(clinical_obj$pd_model)
clinical_ids <- as.character(clinical_obj$sample_ids)
metab_pd <- as.data.frame(cohort_obj$metab_pd)
metab_ids <- as.character(cohort_obj$sample_ids)

if (!identical(clinical_ids, metab_ids)) {
  stop("Clinical and metabolomics sample IDs do not match exactly.")
}
if (anyDuplicated(clinical_ids)) stop("Duplicated sample IDs in the cohort.")

if (is.null(rownames(pd_model))) rownames(pd_model) <- clinical_ids
if (!identical(rownames(pd_model), clinical_ids)) {
  stop("Clinical processed row order does not match sample_ids.")
}
if (is.null(rownames(metab_pd))) rownames(metab_pd) <- metab_ids
if (!identical(rownames(metab_pd), clinical_ids)) {
  stop("Metabolomics row order does not match sample_ids.")
}

clinical_cols <- c("UPDRS_adj", "MOCA_adj", "LEDD_adj")
missing_clinical <- setdiff(clinical_cols, names(pd_model))
if (length(missing_clinical) > 0) {
  stop("Missing required clinical columns: ", paste(missing_clinical, collapse = ", "))
}

resolved_metabs <- resolve_metabolite_columns(metab_pd, top20_metabolites)
if (length(resolved_metabs$missing) > 0) {
  stop(
    "Could not find the following requested metabolites in metabolomics data: ",
    paste(resolved_metabs$missing, collapse = ", ")
  )
}

top20_metab_raw <- metab_pd[, resolved_metabs$resolved, drop = FALSE]

combined_clinical_raw <- as.matrix(pd_model[, clinical_cols, drop = FALSE])
rownames(combined_clinical_raw) <- clinical_ids

combined_clinical_scaled <- scale(combined_clinical_raw)
combined_clinical_scaled <- as.matrix(combined_clinical_scaled)
rownames(combined_clinical_scaled) <- clinical_ids
colnames(combined_clinical_scaled) <- clinical_cols

if (anyNA(combined_clinical_scaled)) stop("Scaled combined clinical matrix contains NA values.")

combined_clinical_affinity_mat <- build_affinity_from_data(
  combined_clinical_scaled,
  sample_ids = clinical_ids,
  K = snf_params$K,
  sigma = snf_params$sigma
)

combined_top20_affinity_mat <- build_affinity_from_data(
  top20_metab_raw,
  sample_ids = clinical_ids,
  K = snf_params$K,
  sigma = snf_params$sigma
)

saveRDS(combined_clinical_raw, file.path(out_dir, "combined_clinical_raw.rds"))
saveRDS(combined_clinical_scaled, file.path(out_dir, "combined_clinical_scaled.rds"))
saveRDS(combined_clinical_affinity_mat, file.path(out_dir, "combined_clinical_affinity_mat.rds"))

saveRDS(top20_metab_raw, file.path(out_dir, "top20_metabolites_raw.rds"))
saveRDS(scale(top20_metab_raw), file.path(out_dir, "top20_metabolites_scaled.rds"))
saveRDS(combined_top20_affinity_mat, file.path(out_dir, "top20_metabolites_affinity_mat.rds"))

write.csv(
  data.frame(
    metric = c("n_samples", "n_clinical_features", "n_metabolites_requested", "n_metabolites_resolved"),
    value = c(nrow(combined_clinical_affinity_mat), ncol(combined_clinical_scaled), length(top20_metabolites), length(resolved_metabs$resolved)),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "top20_combined_build_summary.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    requested_name = top20_metabolites,
    resolved_name = resolved_metabs$resolved,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "top20_metabolite_resolution.csv"),
  row.names = FALSE
)

write.csv(
  make_top20_manifest(
    network_names = c("combined_clinical", "top20_metabolites"),
    network_types = c("clinical", "metabolite_panel"),
    source_columns = c(paste(clinical_cols, collapse = ";"), paste(resolved_metabs$resolved, collapse = ";")),
    output_files = c("combined_clinical_affinity_mat.rds", "top20_metabolites_affinity_mat.rds")
  ),
  file.path(out_dir, "network_manifest.csv"),
  row.names = FALSE
)

png(file.path(out_dir, "combined_clinical_affinity_histogram.png"), width = 2000, height = 1500, res = 300)
hist(combined_clinical_affinity_mat[upper.tri(combined_clinical_affinity_mat)], breaks = 50,
     main = "Combined clinical affinity distribution", xlab = "Affinity")
dev.off()

png(file.path(out_dir, "top20_metabolites_affinity_histogram.png"), width = 2000, height = 1500, res = 300)
hist(combined_top20_affinity_mat[upper.tri(combined_top20_affinity_mat)], breaks = 50,
     main = "Top-20 metabolite affinity distribution", xlab = "Affinity")
dev.off()

message("Top-20 combined affinity build complete.")
message("Outputs saved to: ", out_dir)
message("Resolved metabolites: ", length(resolved_metabs$resolved))

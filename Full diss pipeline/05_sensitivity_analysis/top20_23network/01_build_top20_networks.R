# Top-20 sensitivity analysis
# Build three clinical affinities plus 20 separate metabolite affinities.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))
source(file.path(project_root, "05_sensitivity_analysis", "shared", "top20_snf_helpers.R"))

out_dir <- file.path(paths$sensitivity, "top20_23network")
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

# Clinical networks
clinical_network_specs <- list(
  updrs = "UPDRS_adj",
  moca = "MOCA_adj",
  ledd = "LEDD_adj"
)

clinical_affinity_mats <- list()
clinical_manifest_rows <- list()
for (nm in names(clinical_network_specs)) {
  var_name <- clinical_network_specs[[nm]]
  if (anyNA(pd_model[[var_name]])) stop("Clinical variable contains NA values: ", var_name)

  one_col <- data.frame(value = pd_model[[var_name]], row.names = clinical_ids, check.names = FALSE)
  aff <- build_affinity_from_data(one_col, sample_ids = clinical_ids, K = snf_params$K, sigma = snf_params$sigma)
  clinical_affinity_mats[[nm]] <- aff
  saveRDS(aff, file.path(out_dir, paste0(nm, "_affinity_mat.rds")))

  clinical_manifest_rows[[length(clinical_manifest_rows) + 1L]] <- data.frame(
    network_name = nm,
    network_type = "clinical",
    source_column = var_name,
    output_file = paste0(nm, "_affinity_mat.rds"),
    stringsAsFactors = FALSE
  )
}

# Metabolite networks
metabolite_affinity_mats <- list()
metabolite_manifest_rows <- list()
for (met_name in resolved_metabs$resolved) {
  if (anyNA(metab_pd[[met_name]])) stop("Metabolite column contains NA values: ", met_name)

  slug <- make_slug(met_name)
  one_col <- data.frame(value = metab_pd[[met_name]], row.names = clinical_ids, check.names = FALSE)
  aff <- build_affinity_from_data(one_col, sample_ids = clinical_ids, K = snf_params$K, sigma = snf_params$sigma)
  metabolite_affinity_mats[[met_name]] <- aff
  saveRDS(aff, file.path(out_dir, paste0(slug, "_affinity_mat.rds")))

  metabolite_manifest_rows[[length(metabolite_manifest_rows) + 1L]] <- data.frame(
    network_name = met_name,
    network_type = "metabolite",
    source_column = met_name,
    output_file = paste0(slug, "_affinity_mat.rds"),
    stringsAsFactors = FALSE
  )
}

network_list <- c(clinical_affinity_mats, metabolite_affinity_mats)
network_names <- names(network_list)

saveRDS(network_list, file.path(out_dir, "top20_23network_affinity_list.rds"))
saveRDS(clinical_affinity_mats, file.path(out_dir, "clinical_affinity_mats.rds"))
saveRDS(metabolite_affinity_mats, file.path(out_dir, "metabolite_affinity_mats.rds"))

manifest <- do.call(rbind, c(clinical_manifest_rows, metabolite_manifest_rows))
write.csv(manifest, file.path(out_dir, "network_manifest.csv"), row.names = FALSE)

write.csv(
  data.frame(
    metric = c("n_samples", "n_clinical_networks", "n_metabolite_networks", "n_total_networks", "metabolites_requested", "metabolites_resolved"),
    value = c(nrow(pd_model), length(clinical_affinity_mats), length(metabolite_affinity_mats), length(network_list), length(top20_metabolites), length(resolved_metabs$resolved)),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "top20_23network_build_summary.csv"),
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

message("Top-20 23-network build complete.")
message("Outputs saved to: ", out_dir)
message("Resolved metabolites: ", length(resolved_metabs$resolved))

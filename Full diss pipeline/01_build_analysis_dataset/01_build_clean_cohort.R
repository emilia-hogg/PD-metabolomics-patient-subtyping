# 01_build_analysis_dataset/01_build_clean_cohort.R
# This script builds the analysis cohort by matching metabolomics data to clinical
# metadata, selecting patients with `PATIENT_CASE == "Diagnosis_3y"`, converting
# required clinical variables to numeric form, and applying complete-case
# filtering to the primary clinical variables.
#
# The resulting cohort contains the filtered clinical data, matching metabolomics
# data, and sample IDs. Cohort QC metrics are also saved as RDS and CSV files.
#
# Before running:
# - Set `project_root` to the local project directory.
# - The workspace file specified by `files$workspace_rdata` must exist and contain
#   the `imputation_2_data` object.
# - Clinical metadata must either contain `full_patient_info` in that workspace
#   or be available at `files$clinical_rds`.
#
# Required clinical variables include `Anonymised_sampleID`, `PATIENT_CASE`,
# `UPDRS_III`, `MOCA_total`, `LEDD_total`, `AGE`, and `GENDER`. If `age_onset`
# is available, it is used to calculate disease duration.
#
# Outputs:
# - `analysis_cohort.rds`
# - `patient_metadata.rds`
# - `sample_ids.rds`
# - `cohort_qc_summary.rds`
# - `cohort_qc_summary.csv`

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

# Check that the expected workspace file is available and load the metabolomics data.
if (!file.exists(files$workspace_rdata)) {
  stop("Missing workspace file: ", files$workspace_rdata)
}

load(files$workspace_rdata)

# The metabolomics input is expected to be stored as `imputation_2_data`.
if (!exists("imputation_2_data")) {
  stop("Object 'imputation_2_data' was not found in the loaded workspace.")
}

# Use clinical metadata from the loaded workspace when available; otherwise,
# read it from the clinical RDS file specified in the configuration.
if (exists("full_patient_info")) {
  clinical <- as.data.frame(full_patient_info)
} else if (file.exists(files$clinical_rds)) {
  clinical <- as.data.frame(readRDS(files$clinical_rds))
} else {
  stop("Could not find clinical metadata in the workspace or at: ", files$clinical_rds)
}

metab <- as.data.frame(imputation_2_data)

# Require valid sample IDs in the metabolomics row names so the two data sources
# can be matched.
if (is.null(rownames(metab)) || any(rownames(metab) == "")) {
  stop("Metabolomics matrix must have valid sample IDs in rownames.")
}

if (!"Anonymised_sampleID" %in% names(clinical)) {
  stop("Clinical metadata does not contain 'Anonymised_sampleID'.")
}

# Convert sample IDs to character before matching the clinical and metabolomics data.
clinical$Anonymised_sampleID <- as.character(clinical$Anonymised_sampleID)
rownames(metab) <- as.character(rownames(metab))

# Ensure that each clinical sample ID identifies at most one row.
if (anyDuplicated(clinical$Anonymised_sampleID)) {
  stop("Clinical metadata contains duplicated Anonymised_sampleID values.")
}

# Match metabolomics samples to the corresponding clinical metadata and stop if
# any metabolomics sample cannot be matched.
matched_idx <- match(rownames(metab), clinical$Anonymised_sampleID)
if (anyNA(matched_idx)) {
  missing_ids <- rownames(metab)[is.na(matched_idx)]
  stop("Some metabolomics IDs could not be matched to clinical metadata. First missing ID: ", missing_ids[1])
}

meta <- clinical[matched_idx, , drop = FALSE]
rownames(meta) <- rownames(metab)

# Confirm that the matched clinical metadata and metabolomics data are in the same order.
if (!identical(rownames(meta), rownames(metab))) {
  stop("Row order mismatch after matching clinical metadata to metabolomics.")
}

if (!"PATIENT_CASE" %in% names(meta)) {
  stop("Clinical metadata does not contain 'PATIENT_CASE'.")
}

# Restrict the matched data to rows labelled as `Diagnosis_3y`.
pd_meta <- meta[meta$PATIENT_CASE == "Diagnosis_3y", , drop = FALSE]

if (nrow(pd_meta) == 0) {
  stop("No PD cases found after filtering PATIENT_CASE == 'Diagnosis_3y'.")
}

# Check that the clinical variables required for the complete-case filter are present.
required_cols <- c("UPDRS_III", "MOCA_total", "LEDD_total", "AGE", "GENDER")
missing_cols <- setdiff(required_cols, names(pd_meta))
if (length(missing_cols) > 0) {
  stop("Missing required clinical columns: ", paste(missing_cols, collapse = ", "))
}

# Convert the clinical variables used in the cohort selection to numeric form.
pd_meta$AGE_num <- suppressWarnings(as.numeric(pd_meta$AGE))
pd_meta$UPDRS_III_num <- suppressWarnings(as.numeric(pd_meta$UPDRS_III))
pd_meta$MOCA_total_num <- suppressWarnings(as.numeric(pd_meta$MOCA_total))
pd_meta$LEDD_total_num <- suppressWarnings(as.numeric(pd_meta$LEDD_total))
pd_meta$age_onset_num <- if ("age_onset" %in% names(pd_meta)) {
  suppressWarnings(as.numeric(pd_meta$age_onset))
} else {
  NA_real_
}

# Calculate disease duration from age and age at onset. Negative values are
# treated as missing.
pd_meta$disease_duration <- pd_meta$AGE_num - pd_meta$age_onset_num
pd_meta$disease_duration[pd_meta$disease_duration < 0] <- NA

# Keep rows with complete values for the primary clinical variables and with
# GENDER recorded as either Male or Female.
primary_cols <- c("UPDRS_III_num", "MOCA_total_num", "LEDD_total_num", "AGE_num", "GENDER")
keep <- complete.cases(pd_meta[, primary_cols]) & pd_meta$GENDER %in% c("Male", "Female")
pd_model <- as.data.frame(pd_meta[keep, , drop = FALSE])

if (nrow(pd_model) == 0) {
  stop("No PD rows remained after complete-case filtering on primary clinical variables.")
}

# Use the final sample IDs to extract the corresponding metabolomics rows.
pd_ids <- as.character(pd_model$Anonymised_sampleID)
rownames(pd_model) <- pd_ids

pd_metab <- metab[pd_ids, , drop = FALSE]

# Confirm that the metabolomics rows correspond exactly to the final clinical cohort.
if (!identical(rownames(pd_metab), pd_ids)) {
  stop("Metabolomics row order does not match the final PD cohort.")
}

# Store the final clinical data, metabolomics data and sample IDs together.
analysis_cohort <- list(
  pd_model = pd_model,
  metab_pd = pd_metab,
  sample_ids = pd_ids
)

# Create a small QC summary describing the number of samples retained and
# basic checks on the final metabolomics data.
cohort_qc <- data.frame(
  metric = c(
    "n_total_matched",
    "n_pd_before_complete_case",
    "n_pd_after_complete_case",
    "n_metabolites",
    "any_duplicate_final_ids",
    "any_missing_metabolomics_values"
  ),
  value = c(
    nrow(meta),
    nrow(pd_meta),
    nrow(pd_model),
    ncol(pd_metab),
    anyDuplicated(pd_ids) > 0,
    anyNA(pd_metab)
  ),
  stringsAsFactors = FALSE
)

# Save the final cohort, supporting metadata, sample IDs and QC summary.
saveRDS(analysis_cohort, files$analysis_cohort)
saveRDS(pd_model, files$patient_metadata)
saveRDS(pd_ids, files$sample_ids)
saveRDS(cohort_qc, file.path(paths$analysis_dataset, "cohort_qc_summary.rds"))
write.csv(cohort_qc, file.path(paths$analysis_dataset, "cohort_qc_summary.csv"), row.names = FALSE)

message("Clean cohort saved to: ", files$analysis_cohort)
message("Final cohort size: ", nrow(pd_model))
message("Metabolomics features: ", ncol(pd_metab))

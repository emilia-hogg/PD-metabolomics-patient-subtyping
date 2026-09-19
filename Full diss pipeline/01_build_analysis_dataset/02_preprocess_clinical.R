# 01_build_analysis_dataset/02_preprocess_clinical.R
# This script preprocesses the clinical variables in the analysis cohort by
# residualising UPDRS III, MOCA, LEDD, and disease duration against age and
# gender. The processed clinical data and a summary of the retained samples
# are then saved.
#
# Before running:
# - Set `project_root` to the local project directory.
# - `files$analysis_cohort` must point to an existing analysis cohort created
#   by `01_build_clean_cohort.R`.
#
# Outputs:
# - `clinical_processed.rds`
# - `clinical_preprocess_summary.csv`

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

# Check that the analysis cohort created by the previous preprocessing step exists.
if (!file.exists(files$analysis_cohort)) {
  stop("Missing analysis cohort file: ", files$analysis_cohort)
}

analysis_cohort <- readRDS(files$analysis_cohort)
pd_model <- as.data.frame(analysis_cohort$pd_model)

# Check that all clinical variables required for the preprocessing are present.
required_cols <- c("AGE_num", "UPDRS_III_num", "MOCA_total_num", "LEDD_total_num", "disease_duration", "GENDER")
missing_cols <- setdiff(required_cols, names(pd_model))
if (length(missing_cols) > 0) {
  stop("Missing required columns in pd_model: ", paste(missing_cols, collapse = ", "))
}

# Residualise each clinical variable against age and gender.
pd_model$UPDRS_adj <- residuals(lm(UPDRS_III_num ~ AGE_num + GENDER, data = pd_model, na.action = na.exclude))
pd_model$MOCA_adj <- residuals(lm(MOCA_total_num ~ AGE_num + GENDER, data = pd_model, na.action = na.exclude))
pd_model$LEDD_adj <- residuals(lm(LEDD_total_num ~ AGE_num + GENDER, data = pd_model, na.action = na.exclude))
pd_model$duration_adj <- residuals(lm(disease_duration ~ AGE_num + GENDER, data = pd_model, na.action = na.exclude))

# Store the processed clinical data together with the sample IDs from the analysis cohort.
clinical_processed <- list(
  pd_model = pd_model,
  sample_ids = analysis_cohort$sample_ids
)

# Save the processed clinical dataset.
saveRDS(clinical_processed, files$clinical_processed)

# Save a summary of the number of retained samples and non-missing disease-duration values.
write.csv(
  data.frame(
    metric = c("n_final_samples", "n_nonmissing_duration"),
    value = c(nrow(pd_model), sum(!is.na(pd_model$disease_duration)))
  ),
  file.path(paths$analysis_dataset, "clinical_preprocess_summary.csv"),
  row.names = FALSE
)

message("Clinical preprocessing complete.")
message("Rows retained: ", nrow(pd_model))

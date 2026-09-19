# 01_build_analysis_dataset/04_baseline_characteristics.R

# This script generates baseline characteristics for the analysis cohort,
# including summaries of clinical variables, sex, and metabolomics data.
# Continuous variables are reported using mean (SD), median [IQR], and range,
# together with non-missing and missing counts. A publication-formatted table
# and a QC summary are also created.
#
# The analysis cohort is loaded from `files$analysis_cohort` and is not
# redefined by this script. If `files$metabolomics_processed` exists, the
# processed metabolomics data are used to report the number of retained
# metabolite features; otherwise this is derived from the raw cohort matrix.
#
# Before running:
# - Set `project_root` to the local project directory.
# - `files$analysis_cohort` must point to an existing analysis cohort.
# - `files$metabolomics_processed` is optional; if present, its sample IDs
#   must match those in the analysis cohort.
#
# Outputs:
# - `baseline_characteristics.csv`
# - `baseline_characteristics_publication.csv`
# - `baseline_characteristics.rds`
# - `baseline_characteristics_summary.txt`

# -----------------------------
# 1. Project setup
# -----------------------------
project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

config_file <- file.path(project_root, "00_config", "config.R")
if (!file.exists(config_file)) {
  stop("Missing config file: ", config_file)
}
source(config_file)

analysis_cohort_file <- files$analysis_cohort
if (!file.exists(analysis_cohort_file)) {
  stop("Missing analysis cohort file: ", analysis_cohort_file)
}

out_dir <- paths$analysis_dataset
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Helper functions
# -----------------------------
stop_if_missing <- function(x, object_name) {
  if (is.null(x)) {
    stop(object_name, " is NULL.")
  }
  invisible(TRUE)
}

# Return the first available variable name from a set of alternative names,
# stopping if none of the expected variables are present.
first_present <- function(x, candidates, object_name) {
  for (nm in candidates) {
    if (nm %in% names(x)) {
      return(nm)
    }
  }
  stop(
    "Could not find any of these variables in ", object_name, ": ",
    paste(candidates, collapse = ", ")
  )
}

clean_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[nchar(x) == 0] <- NA_character_
  x
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

fmt_num <- function(x, digits = 1) {
  if (!is.finite(x)) {
    return(NA_character_)
  }
  formatC(x, format = "f", digits = digits)
}

fmt_mean_sd <- function(x, digits = 1) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_character_)
  sd_x <- if (length(x) >= 2) stats::sd(x) else NA_real_
  paste0(fmt_num(mean(x), digits), " (", fmt_num(sd_x, digits), ")")
}

fmt_median_iqr <- function(x, digits = 1) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_character_)
  q <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  paste0(
    fmt_num(q[2], digits),
    " [",
    fmt_num(q[1], digits),
    ", ",
    fmt_num(q[3], digits),
    "]"
  )
}

fmt_min_max <- function(x, digits = 1) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_character_)
  paste0(
    fmt_num(min(x), digits),
    " to ",
    fmt_num(max(x), digits)
  )
}

pretty_label <- function(var_name) {
  switch(
    var_name,
    AGE_num = "Age (years)",
    AGE = "Age (years)",
    UPDRS_III_num = "UPDRS III",
    UPDRS_III = "UPDRS III",
    MOCA_total_num = "MoCA total",
    MOCA_total = "MoCA total",
    LEDD_total_num = "LEDD",
    LEDD_total = "LEDD",
    disease_duration = "Disease duration (years)",
    age_onset_num = "Age at onset (years)",
    age_onset = "Age at onset (years)",
    var_name
  )
}

summarise_continuous <- function(df, var_name, section = "Clinical characteristics", digits = 1) {
  if (!var_name %in% names(df)) {
    return(NULL)
  }

  raw_x <- safe_numeric(df[[var_name]])
  finite_x <- raw_x[is.finite(raw_x)]

  data.frame(
    section = section,
    characteristic = pretty_label(var_name),
    level = "",
    n_nonmissing = sum(is.finite(raw_x)),
    n_missing = sum(!is.finite(raw_x)),
    percent = NA_real_,
    mean_sd = fmt_mean_sd(finite_x, digits),
    median_iqr = fmt_median_iqr(finite_x, digits),
    min_max = fmt_min_max(finite_x, digits),
    value = NA_character_,
    notes = if (length(finite_x) == 0) "No finite values available" else "",
    stringsAsFactors = FALSE
  )
}

summarise_binary_categorical <- function(df, var_name, label = "Sex", level_order = c("Male", "Female")) {
  if (!var_name %in% names(df)) {
    return(NULL)
  }

  x <- clean_ids(df[[var_name]])
  total_n <- length(x)
  x_nonmissing <- x[!is.na(x)]

  if (length(x_nonmissing) == 0) {
    return(NULL)
  }

  levels_all <- unique(c(level_order, sort(unique(x_nonmissing))))
  tab <- table(factor(x_nonmissing, levels = levels_all))

  out <- lapply(names(tab), function(lv) {
    n <- as.integer(tab[[lv]])
    data.frame(
      section = "Sex",
      characteristic = label,
      level = lv,
      n_nonmissing = n,
      n_missing = sum(is.na(x)),
      percent = 100 * n / total_n,
      mean_sd = NA_character_,
      median_iqr = NA_character_,
      min_max = NA_character_,
      value = as.character(n),
      notes = "",
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

# Summarise metabolomics sample and feature counts, including zero-variance
# filtering and the number of features retained after preprocessing.
summarise_metabolomics <- function(raw_metab, processed_metab = NULL, sample_ids = NULL) {
  raw_metab <- as.data.frame(raw_metab)

  if (is.null(rownames(raw_metab)) && !is.null(sample_ids)) {
    rownames(raw_metab) <- sample_ids
  }

  raw_mat <- as.matrix(raw_metab)
  raw_n_patients <- nrow(raw_mat)
  raw_n_features <- ncol(raw_mat)
  raw_missing <- sum(!is.finite(raw_mat))

  zero_var_keep <- apply(raw_mat, 2, function(x) {
    x <- safe_numeric(x)
    x <- x[is.finite(x)]
    if (length(x) < 2) {
      return(FALSE)
    }
    stats::var(x) > 0
  })

  n_zero_var_removed <- sum(!zero_var_keep)
  n_retained_from_raw <- sum(zero_var_keep)

  if (!is.null(processed_metab)) {
    processed_metab <- as.data.frame(processed_metab)
    proc_mat <- as.matrix(processed_metab)
    if (nrow(proc_mat) != raw_n_patients) {
      stop("Processed metabolomics row count does not match raw metabolomics row count.")
    }
    processed_n_features <- ncol(proc_mat)
  } else {
    processed_n_features <- n_retained_from_raw
  }

  rows <- list(
    data.frame(
      section = "Metabolomics",
      characteristic = "Patients with metabolomics data",
      level = "",
      n_nonmissing = raw_n_patients,
      n_missing = 0L,
      percent = NA_real_,
      mean_sd = NA_character_,
      median_iqr = NA_character_,
      min_max = NA_character_,
      value = as.character(raw_n_patients),
      notes = "Matched to the analysis cohort",
      stringsAsFactors = FALSE
    ),
    data.frame(
      section = "Metabolomics",
      characteristic = "Metabolites before zero-variance filtering",
      level = "",
      n_nonmissing = raw_n_features,
      n_missing = 0L,
      percent = NA_real_,
      mean_sd = NA_character_,
      median_iqr = NA_character_,
      min_max = NA_character_,
      value = as.character(raw_n_features),
      notes = "",
      stringsAsFactors = FALSE
    ),
    data.frame(
      section = "Metabolomics",
      characteristic = "Zero-variance metabolites removed",
      level = "",
      n_nonmissing = n_zero_var_removed,
      n_missing = 0L,
      percent = NA_real_,
      mean_sd = NA_character_,
      median_iqr = NA_character_,
      min_max = NA_character_,
      value = as.character(n_zero_var_removed),
      notes = "",
      stringsAsFactors = FALSE
    ),
    data.frame(
      section = "Metabolomics",
      characteristic = "Metabolites retained after filtering",
      level = "",
      n_nonmissing = processed_n_features,
      n_missing = 0L,
      percent = NA_real_,
      mean_sd = NA_character_,
      median_iqr = NA_character_,
      min_max = NA_character_,
      value = as.character(processed_n_features),
      notes = if (is.null(processed_metab)) {
        "Derived from the raw cohort matrix"
      } else {
        "Taken from metabolomics_processed.rds"
      },
      stringsAsFactors = FALSE
    ),
    data.frame(
      section = "Metabolomics",
      characteristic = "Scaling applied",
      level = "",
      n_nonmissing = NA_integer_,
      n_missing = NA_integer_,
      percent = NA_real_,
      mean_sd = NA_character_,
      median_iqr = NA_character_,
      min_max = NA_character_,
      value = "Yes",
      notes = "Standardised feature-wise after filtering",
      stringsAsFactors = FALSE
    ),
    data.frame(
      section = "Metabolomics",
      characteristic = "Missing values in raw metabolomics",
      level = "",
      n_nonmissing = raw_missing,
      n_missing = 0L,
      percent = if (length(raw_mat) > 0) 100 * raw_missing / length(raw_mat) else NA_real_,
      mean_sd = NA_character_,
      median_iqr = NA_character_,
      min_max = NA_character_,
      value = as.character(raw_missing),
      notes = "",
      stringsAsFactors = FALSE
    )
  )

  do.call(rbind, rows)
}

# Create a publication-formatted version of the baseline table with selected
# columns, reader-friendly column names, and formatted percentages.
make_publication_table <- function(df) {
  pub <- df[, c(
    "section",
    "characteristic",
    "level",
    "n_nonmissing",
    "percent",
    "mean_sd",
    "median_iqr",
    "min_max",
    "n_missing",
    "notes"
  )]

  names(pub) <- c(
    "Section",
    "Characteristic",
    "Level",
    "N",
    "Percent",
    "Mean (SD)",
    "Median [IQR]",
    "Min to Max",
    "Missing",
    "Notes"
  )

  pub$Percent <- ifelse(
    is.na(pub$Percent),
    NA_character_,
    paste0(formatC(pub$Percent, format = "f", digits = 1), "%")
  )

  pub$N <- ifelse(is.na(pub$N), NA_character_, as.character(pub$N))
  pub$Missing <- ifelse(is.na(pub$Missing), NA_character_, as.character(pub$Missing))

  pub
}

# -----------------------------
# 3. Load analysis cohort
# -----------------------------
analysis_cohort <- readRDS(analysis_cohort_file)
stop_if_missing(analysis_cohort, "analysis_cohort")

required_objects <- c("pd_model", "metab_pd", "sample_ids")
missing_objects <- setdiff(required_objects, names(analysis_cohort))
if (length(missing_objects) > 0) {
  stop(
    "analysis_cohort.rds is missing required elements: ",
    paste(missing_objects, collapse = ", ")
  )
}

pd_model <- as.data.frame(analysis_cohort$pd_model)
metab_pd_raw <- as.data.frame(analysis_cohort$metab_pd)
sample_ids <- clean_ids(analysis_cohort$sample_ids)

if (length(sample_ids) == 0) {
  stop("sample_ids is empty.")
}
if (anyNA(sample_ids)) {
  stop("sample_ids contains missing values.")
}
if (anyDuplicated(sample_ids)) {
  stop("sample_ids contains duplicates.")
}

if (nrow(pd_model) != length(sample_ids)) {
  stop("pd_model row count does not match sample_ids length.")
}

if (is.null(rownames(pd_model))) {
  rownames(pd_model) <- sample_ids
} else {
  rownames(pd_model) <- clean_ids(rownames(pd_model))
}

if (is.null(rownames(metab_pd_raw))) {
  rownames(metab_pd_raw) <- sample_ids
}

if (!identical(clean_ids(rownames(metab_pd_raw)), sample_ids)) {
  stop("metab_pd row order does not match sample_ids.")
}

if (!identical(clean_ids(rownames(pd_model)), sample_ids)) {
  rownames(pd_model) <- sample_ids
}

# -----------------------------
# 4. Identify variables to summarise
# -----------------------------
age_var <- first_present(pd_model, c("AGE_num", "AGE"), "pd_model")
updrs_var <- first_present(pd_model, c("UPDRS_III_num", "UPDRS_III"), "pd_model")
moca_var <- first_present(pd_model, c("MOCA_total_num", "MOCA_total"), "pd_model")
ledd_var <- first_present(pd_model, c("LEDD_total_num", "LEDD_total"), "pd_model")
gender_var <- first_present(pd_model, c("GENDER", "Sex", "sex"), "pd_model")

duration_var <- if ("disease_duration" %in% names(pd_model)) "disease_duration" else NA_character_
age_onset_var <- if ("age_onset_num" %in% names(pd_model)) {
  "age_onset_num"
} else if ("age_onset" %in% names(pd_model)) {
  "age_onset"
} else {
  NA_character_
}

# -----------------------------
# 5. Optional metabolomics processed object
# -----------------------------
metabolomics_processed_file <- files$metabolomics_processed
processed_metab <- NULL

if (file.exists(metabolomics_processed_file)) {
  metabolomics_processed <- readRDS(metabolomics_processed_file)
  if (!all(c("metab_pd", "metab_scaled", "sample_ids") %in% names(metabolomics_processed))) {
    stop(
      "metabolomics_processed.rds must contain metab_pd, metab_scaled, and sample_ids."
    )
  }

  proc_ids <- clean_ids(metabolomics_processed$sample_ids)
  if (!identical(proc_ids, sample_ids)) {
    stop("Sample IDs in metabolomics_processed.rds do not match the analysis cohort.")
  }

  processed_metab <- metabolomics_processed$metab_pd
  if (!identical(clean_ids(rownames(processed_metab)), sample_ids)) {
    stop("Processed metabolomics row order does not match sample_ids.")
  }
}

# -----------------------------
# 6. Build the baseline table
# -----------------------------
baseline_rows <- list()

baseline_rows[[length(baseline_rows) + 1L]] <- data.frame(
  section = "Cohort",
  characteristic = "Total analysis cohort",
  level = "",
  n_nonmissing = nrow(pd_model),
  n_missing = 0L,
  percent = 100,
  mean_sd = NA_character_,
  median_iqr = NA_character_,
  min_max = NA_character_,
  value = as.character(nrow(pd_model)),
  notes = "Complete-case cohort for UPDRS III, MoCA, LEDD, age, and sex",
  stringsAsFactors = FALSE
)

continuous_candidates <- c(age_var, updrs_var, moca_var, ledd_var, duration_var, age_onset_var)
continuous_candidates <- continuous_candidates[!is.na(continuous_candidates)]
continuous_candidates <- unique(continuous_candidates)

for (v in continuous_candidates) {
  row <- summarise_continuous(pd_model, v)
  if (!is.null(row)) {
    baseline_rows[[length(baseline_rows) + 1L]] <- row
  }
}

sex_row <- summarise_binary_categorical(
  pd_model,
  gender_var,
  label = "Sex",
  level_order = c("Male", "Female")
)
if (!is.null(sex_row)) {
  baseline_rows[[length(baseline_rows) + 1L]] <- sex_row
}

metab_row <- summarise_metabolomics(
  raw_metab = metab_pd_raw,
  processed_metab = processed_metab,
  sample_ids = sample_ids
)
baseline_rows[[length(baseline_rows) + 1L]] <- metab_row

baseline_table <- do.call(rbind, baseline_rows)

section_order <- c("Cohort", "Clinical characteristics", "Sex", "Metabolomics")
baseline_table$section <- factor(baseline_table$section, levels = section_order)
baseline_table <- baseline_table[order(baseline_table$section, baseline_table$characteristic, baseline_table$level), ]
baseline_table$section <- as.character(baseline_table$section)

publication_table <- make_publication_table(baseline_table)

# -----------------------------
# 7. QC summary object
# -----------------------------
qc_summary <- data.frame(
  metric = c(
    "cohort_size",
    "continuous_variables_reported",
    "sex_levels_reported",
    "raw_metabolomics_features",
    "metabolomics_features_retained",
    "zero_variance_metabolites_removed",
    "metabolomics_processed_file_found",
    "raw_metabolomics_missing_values"
  ),
  value = c(
    nrow(pd_model),
    length(continuous_candidates),
    if (!is.null(sex_row)) length(unique(sex_row$level)) else 0L,
    ncol(metab_pd_raw),
    as.integer(baseline_table$value[baseline_table$characteristic == "Metabolites retained after filtering"][1]),
    as.integer(baseline_table$value[baseline_table$characteristic == "Zero-variance metabolites removed"][1]),
    file.exists(metabolomics_processed_file),
    sum(!is.finite(as.matrix(metab_pd_raw)))
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 8. Save outputs
# -----------------------------
csv_file <- file.path(out_dir, "baseline_characteristics.csv")
pub_csv_file <- file.path(out_dir, "baseline_characteristics_publication.csv")
rds_file <- file.path(out_dir, "baseline_characteristics.rds")
txt_file <- file.path(out_dir, "baseline_characteristics_summary.txt")

write.csv(baseline_table, csv_file, row.names = FALSE)
write.csv(publication_table, pub_csv_file, row.names = FALSE)

saveRDS(
  list(
    baseline_table = baseline_table,
    publication_table = publication_table,
    qc_summary = qc_summary,
    cohort_n = nrow(pd_model),
    sample_ids = sample_ids
  ),
  rds_file
)

summary_lines <- c(
  "Baseline characteristics completed successfully.",
  paste0("Cohort size: ", nrow(pd_model)),
  paste0("Rows in baseline table: ", nrow(baseline_table)),
  paste0("Clinical variables summarised: ", paste(continuous_candidates, collapse = ", ")),
  paste0("Sex variable used: ", gender_var),
  paste0("Metabolomics processed file found: ", file.exists(metabolomics_processed_file)),
  paste0("Raw metabolomics features: ", ncol(metab_pd_raw)),
  paste0(
    "Metabolites retained after filtering: ",
    as.integer(baseline_table$value[baseline_table$characteristic == "Metabolites retained after filtering"][1])
  ),
  paste0(
    "Zero-variance metabolites removed: ",
    as.integer(baseline_table$value[baseline_table$characteristic == "Zero-variance metabolites removed"][1])
  )
)
writeLines(summary_lines, txt_file)

# -----------------------------
# 9. Final messages
# -----------------------------
message("Baseline characteristics complete.")
message("Outputs saved to: ", out_dir)
message("Cohort size: ", nrow(pd_model))
message("Rows in baseline table: ", nrow(baseline_table))
message("Publication CSV created: ", pub_csv_file)

# 03_network_construction/01_build_clinical_networks.R
# This script builds three clinical affinity matrices from the processed
# clinical data: UPDRS, MoCA, and LEDD. The clinical variables are residualised
# against age and gender, standardised, and converted to affinity matrices using
# the SNFtool distance and affinity functions.
#
# Fixed settings:
# - SNF K = `snf_params$K`
# - SNF sigma = `snf_params$sigma`
#
# Before running:
# - Set `project_root` to the local project directory.
# - `files$clinical_processed` must point to an existing processed clinical file.
#
# Outputs are saved under `paths$clinical_networks/clinical_network/`,
# including the three affinity matrices, QC histograms, and a summary file.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

if (!file.exists(files$clinical_processed)) {
  stop("Missing clinical processed file: ", files$clinical_processed)
}

clinical_obj <- readRDS(files$clinical_processed)

if (!all(c("pd_model", "sample_ids") %in% names(clinical_obj))) {
  stop("clinical_processed.rds must contain 'pd_model' and 'sample_ids'.")
}

pd_model <- as.data.frame(clinical_obj$pd_model)
pd_ids <- as.character(clinical_obj$sample_ids)

required_cols <- c(
  "Anonymised_sampleID",
  "AGE_num",
  "UPDRS_III_num",
  "MOCA_total_num",
  "LEDD_total_num",
  "GENDER",
  "UPDRS_adj",
  "MOCA_adj",
  "LEDD_adj"
)

missing_cols <- setdiff(required_cols, names(pd_model))
if (length(missing_cols) > 0) {
  stop("Missing required columns in clinical processed data: ", paste(missing_cols, collapse = ", "))
}

if (anyDuplicated(pd_model$Anonymised_sampleID)) {
  stop("Duplicated Anonymised_sampleID values in clinical processed data.")
}

if (!identical(as.character(pd_model$Anonymised_sampleID), pd_ids)) {
  stop("Clinical processed row order does not match sample_ids.")
}

out_dir <- file.path(paths$clinical_networks, "clinical_network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(filename, plot_fun, width = 2000, height = 1500, res = 300) {
  png(filename, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  plot_fun()
}

make_raw_histogram <- function(x, out_file, main_title, xlab_text, breaks = 30) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    stop("No finite values available for histogram: ", main_title)
  }
  save_png(out_file, plot_fun = function() {
    hist(
      x,
      breaks = breaks,
      main = main_title,
      xlab = xlab_text,
      col = "grey80",
      border = "white"
    )
  })
}

save_affinity_histograms <- function(values, out_dir, file_prefix, title, xlab = "Affinity") {
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    stop("No finite values available for: ", file_prefix)
  }

  cat("\n", title, "\n", sep = "")
  print(summary(values))
  cat("Proportion exactly zero: ", mean(values == 0), "\n", sep = "")
  cat("Quantiles:\n")
  print(quantile(values, probs = c(0, 0.5, 0.9, 0.95, 0.99, 1), na.rm = TRUE))

  save_png(file.path(out_dir, paste0(file_prefix, "_histogram_linear.png")), plot_fun = function() {
    hist(
      values,
      breaks = 80,
      main = paste0(title, " (linear scale)"),
      xlab = xlab,
      col = "grey80",
      border = "white"
    )
  })

  positive_values <- values[values > 0]
  if (length(positive_values) > 0) {
    save_png(file.path(out_dir, paste0(file_prefix, "_histogram_log10.png")), plot_fun = function() {
      hist(
        log10(positive_values),
        breaks = 80,
        main = paste0(title, " (log10 scale)"),
        xlab = paste0("log10(", xlab, ")"),
        col = "grey80",
        border = "white"
      )
    })
  }

  zoom_upper <- as.numeric(quantile(values, 0.99, na.rm = TRUE))
  if (!is.finite(zoom_upper) || zoom_upper <= 0) {
    zoom_upper <- max(values, na.rm = TRUE)
  }
  values_zoom <- values[values <= zoom_upper]

  save_png(file.path(out_dir, paste0(file_prefix, "_histogram_zoomed_linear.png")), plot_fun = function() {
    hist(
      values_zoom,
      breaks = 80,
      main = paste0(title, " (linear scale, up to 99th percentile)"),
      xlab = xlab,
      xlim = c(0, zoom_upper),
      col = "grey80",
      border = "white"
    )
  })
}

residualize_clinical_variable <- function(data, outcome_name) {
  fit <- lm(
    as.formula(paste0(outcome_name, " ~ AGE_num + GENDER")),
    data = data,
    na.action = na.exclude
  )
  residuals(fit)
}

make_single_variable_affinity <- function(data, sample_ids, var_name, base_name, out_dir) {
  one_var_df <- data.frame(value = data[[var_name]], row.names = sample_ids, check.names = FALSE)
  one_var_mat <- as.matrix(one_var_df)

  x_scaled <- scale(one_var_mat)
  x_scaled <- as.matrix(x_scaled)
  rownames(x_scaled) <- sample_ids

  if (anyNA(x_scaled)) {
    stop("Scaling produced NA values for ", base_name)
  }

  K_eff <- max(1, min(snf_params$K, nrow(x_scaled) - 1))
  x_dist <- SNFtool::dist2(x_scaled, x_scaled)
  saveRDS(x_dist, file.path(out_dir, paste0(base_name, "_dist_snf.rds")))
  x_affinity <- SNFtool::affinityMatrix(x_dist, K = K_eff, sigma = snf_params$sigma)

  rownames(x_affinity) <- sample_ids
  colnames(x_affinity) <- sample_ids

  if (!isTRUE(all.equal(x_affinity, t(x_affinity)))) {
    stop("Affinity matrix is not symmetric for ", base_name)
  }
  if (anyNA(x_affinity)) {
    stop("Affinity matrix contains NA values for ", base_name)
  }

  save_affinity_histograms(
    values = x_affinity[upper.tri(x_affinity)],
    out_dir = out_dir,
    file_prefix = base_name,
    title = paste0(toupper(base_name), " Affinity Distribution"),
    xlab = "Affinity"
  )

  saveRDS(x_affinity, file.path(out_dir, paste0(base_name, "_affinity_mat.rds")))
  invisible(x_affinity)
}

pd_model$UPDRS_adj <- residualize_clinical_variable(pd_model, "UPDRS_III_num")
pd_model$MOCA_adj <- residualize_clinical_variable(pd_model, "MOCA_total_num")
pd_model$LEDD_adj <- residualize_clinical_variable(pd_model, "LEDD_total_num")

make_raw_histogram(pd_model$UPDRS_III_num, file.path(out_dir, "raw_updrs_histogram.png"), "Raw UPDRS Histogram", "UPDRS III")
make_raw_histogram(pd_model$MOCA_total_num, file.path(out_dir, "raw_moca_histogram.png"), "Raw MoCA Histogram", "MoCA total")
make_raw_histogram(pd_model$LEDD_total_num, file.path(out_dir, "raw_ledd_histogram.png"), "Raw LEDD Histogram", "LEDD total")

make_raw_histogram(pd_model$UPDRS_adj, file.path(out_dir, "adjusted_updrs_histogram.png"), "Adjusted UPDRS Histogram", "UPDRS residuals")
make_raw_histogram(pd_model$MOCA_adj, file.path(out_dir, "adjusted_moca_histogram.png"), "Adjusted MoCA Histogram", "MoCA residuals")
make_raw_histogram(pd_model$LEDD_adj, file.path(out_dir, "adjusted_ledd_histogram.png"), "Adjusted LEDD Histogram", "LEDD residuals")

clinical_network_specs <- list(
  updrs = "UPDRS_adj",
  moca = "MOCA_adj",
  ledd = "LEDD_adj"
)

clinical_affinity_mats <- list()

for (nm in names(clinical_network_specs)) {
  message("Building clinical network: ", nm)
  clinical_affinity_mats[[nm]] <- make_single_variable_affinity(
    data = pd_model,
    sample_ids = pd_ids,
    var_name = clinical_network_specs[[nm]],
    base_name = nm,
    out_dir = out_dir
  )
}

saveRDS(pd_model, file.path(out_dir, "pd_model.rds"))
saveRDS(clinical_network_specs, file.path(out_dir, "clinical_network_specs.rds"))
saveRDS(clinical_affinity_mats, file.path(out_dir, "clinical_affinity_mats.rds"))

write.csv(
  data.frame(
    metric = c("n_patients", "n_clinical_networks"),
    value = c(nrow(pd_model), length(clinical_network_specs))
  ),
  file.path(out_dir, "clinical_network_summary.csv"),
  row.names = FALSE
)

message("Clinical network construction complete.")
message("Outputs saved to: ", out_dir)

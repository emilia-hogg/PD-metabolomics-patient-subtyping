# 03_network_construction/02_build_metabolomics_network.R

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

if (!file.exists(files$analysis_cohort)) {
  stop("Missing analysis cohort file: ", files$analysis_cohort)
}

analysis_cohort <- readRDS(files$analysis_cohort)

if (!all(c("metab_pd", "sample_ids") %in% names(analysis_cohort))) {
  stop("analysis_cohort.rds must contain 'metab_pd' and 'sample_ids'.")
}

metab_pd <- as.data.frame(analysis_cohort$metab_pd)
pd_ids <- as.character(analysis_cohort$sample_ids)

if (anyDuplicated(pd_ids)) {
  stop("Duplicated sample IDs in analysis cohort.")
}
if (!identical(rownames(metab_pd), pd_ids)) {
  stop("Metabolomics row order does not match sample_ids.")
}
if (anyNA(metab_pd)) {
  stop("Metabolomics matrix contains missing values before filtering.")
}

out_dir <- file.path(paths$clinical_networks, "metabolomics_network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(filename, plot_fun, width = 2000, height = 1500, res = 300) {
  png(filename, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  plot_fun()
}

save_affinity_histograms <- function(values, out_dir, file_prefix, title, xlab = "Affinity") {
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    stop("No finite values available for: ", file_prefix)
  }

  cat("\n", title, "\n", sep = "")
  print(summary(values))
  cat("Proportion exactly zero: ", mean(values == 0), "\n", sep = "")
  cat("Number of non-zero values: ", length(values), "\n", sep = "")
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

# Zero-variance filter only; no patient removal after cohort definition.
zero_var_keep <- apply(metab_pd, 2, var, na.rm = TRUE) > 0
metab_pd <- metab_pd[, zero_var_keep, drop = FALSE]
if (ncol(metab_pd) == 0) {
  stop("No metabolomics features remained after zero-variance filtering.")
}

metab_scaled <- scale(metab_pd)
metab_scaled <- as.matrix(metab_scaled)
rownames(metab_scaled) <- pd_ids
colnames(metab_scaled) <- colnames(metab_pd)

if (anyNA(metab_scaled)) {
  stop("Scaled metabolomics matrix contains NA values.")
}
if (!identical(rownames(metab_scaled), pd_ids)) {
  stop("Scaled metabolomics row order does not match sample_ids.")
}

saveRDS(metab_pd, file.path(out_dir, "metab_pd.rds"))
saveRDS(metab_scaled, file.path(out_dir, "metab_scaled.rds"))

# SNFtool inputs
metab_dist_snf <- SNFtool::dist2(metab_scaled, metab_scaled)
rownames(metab_dist_snf) <- pd_ids
colnames(metab_dist_snf) <- pd_ids
saveRDS(metab_dist_snf, file.path(out_dir, "metab_dist_snf.rds"))

K_eff <- max(1, min(snf_params$K, nrow(metab_scaled) - 1))
metab_affinity_mat <- SNFtool::affinityMatrix(
  metab_dist_snf,
  K = K_eff,
  sigma = snf_params$sigma
)
rownames(metab_affinity_mat) <- pd_ids
colnames(metab_affinity_mat) <- pd_ids
saveRDS(metab_affinity_mat, file.path(out_dir, "metab_affinity_mat.rds"))

if (!isTRUE(all.equal(metab_affinity_mat, t(metab_affinity_mat)))) {
  stop("Metabolomics affinity matrix is not symmetric.")
}
if (anyNA(metab_affinity_mat)) {
  stop("Metabolomics affinity matrix contains NA values.")
}

save_affinity_histograms(
  values = metab_affinity_mat[upper.tri(metab_affinity_mat)],
  out_dir = out_dir,
  file_prefix = "metabolomics_affinity",
  title = "Metabolomics Affinity Distribution",
  xlab = "Affinity"
)

write.csv(
  data.frame(
    metric = c("n_patients", "n_features_after_zero_var", "K_effective"),
    value = c(nrow(metab_scaled), ncol(metab_scaled), K_eff)
  ),
  file.path(out_dir, "metabolomics_network_summary.csv"),
  row.names = FALSE
)

message("Metabolomics network construction complete.")
message("Outputs saved to: ", out_dir)
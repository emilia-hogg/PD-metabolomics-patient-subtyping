# ============================================================
# 03b_affinity_summary_all.R
#
# Reuses the upper-bound calculation from 03_network_qc.R 
# Extended to the three clinical networks

# This script summarises the affinity matrices for the three clinical networks,
# the metabolomics network, and any available combined/fused networks. It reports
# affinity distributions and, for raw affinity matrices, their ratio to a
# theoretical upper bound.
#
# Fixed settings:
# - The upper-bound calculation uses `snf_params$K` and `snf_params$sigma`.
#
# Before running:
# - Set `project_root` to the local project directory.
# - The clinical and metabolomics network outputs must already have been
#   generated.
#
# Output:
# - `affinity_summary_all_networks.csv`
# - `affinity_summary_all_networks.rds`
# ============================================================

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

clinical_dir     <- file.path(paths$clinical_networks, "clinical_network")
metabolomics_dir <- file.path(paths$clinical_networks, "metabolomics_network")
qc_dir           <- file.path(paths$clinical_networks, "qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

pd_model  <- readRDS(file.path(clinical_dir, "pd_model.rds"))
pd_ids    <- as.character(pd_model$Anonymised_sampleID)
N         <- length(pd_ids)

# ------------------------------------------------------------
# 1. Upper-bound helper
# ------------------------------------------------------------
compute_affinity_upper_bound <- function(diff, K = snf_params$K, sigma = snf_params$sigma) {
  diff <- as.matrix(diff)
  diff <- (diff + t(diff)) / 2
  diag(diff) <- 0

  K_eff <- max(1, min(K, nrow(diff) - 2))
  sortedColumns <- as.matrix(t(apply(diff, 2, sort)))
  finiteMean <- function(x) mean(x[is.finite(x)])

  means <- apply(sortedColumns[, (1:K_eff) + 1, drop = FALSE], 1, finiteMean) + .Machine$double.eps
  avg <- function(x, y) (x + y) / 2
  Sig <- outer(means, means, avg) / 3 * 2 + diff / 3 + .Machine$double.eps
  Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps

  upper_bound <- dnorm(0, mean = 0, sd = sigma * Sig, log = FALSE)
  dimnames(Sig) <- dimnames(upper_bound) <- dimnames(diff)
  list(Sig = Sig, upper_bound = upper_bound, K_eff = K_eff)
}

# ------------------------------------------------------------
# 2. Rebuild the distance matrix for a single clinical variable
# ------------------------------------------------------------
single_var_dist <- function(data, var_name, sample_ids) {
  m <- as.matrix(data.frame(value = data[[var_name]], row.names = sample_ids,
                            check.names = FALSE))
  xs <- as.matrix(scale(m))
  rownames(xs) <- sample_ids
  d <- SNFtool::dist2(xs, xs)
  dimnames(d) <- list(sample_ids, sample_ids)
  d
}

# ------------------------------------------------------------
# 3. Summarise one affinity matrix
#    row_normalised = TRUE for the fused SNF output, where the
#    diagonal is fixed at 0.5 and the off-diagonal mean is
#    determined by construction rather than estimated.
# ------------------------------------------------------------
summarise_affinity <- function(W, dist_mat = NULL, label,
                               row_normalised = FALSE, n_features = NA) {
  W <- as.matrix(W)
  ut <- upper.tri(W)
  off <- W[ut]

  out <- c(
    n_patients                = nrow(W),
    n_features                = n_features,
    n_edges_upper_triangle    = sum(ut),
    affinity_min              = min(W),
    affinity_max              = max(W),
    offdiag_min               = min(off),
    offdiag_median            = median(off),
    offdiag_mean              = mean(off),
    offdiag_q95               = as.numeric(quantile(off, 0.95, names = FALSE)),
    offdiag_q99               = as.numeric(quantile(off, 0.99, names = FALSE)),
    offdiag_max               = max(off),
    n_above_0.001             = sum(off > 0.001),
    prop_above_0.001_pct      = 100 * mean(off > 0.001),
    structural_mean_if_rownorm = if (row_normalised) 0.5 / (nrow(W) - 1) else NA_real_
  )

  # Ratio to the theoretical upper bound. Only meaningful for a raw
  # affinityMatrix() output; the fused matrix has been row-normalised
  # and no longer sits on the kernel's original scale.
  if (!is.null(dist_mat) && !row_normalised) {
    ub  <- compute_affinity_upper_bound(dist_mat)$upper_bound
    rat <- W / ub
    ubo <- ub[ut]
    rto <- rat[ut]
    out <- c(out,
      upper_bound_offdiag_min    = min(ubo),
      upper_bound_offdiag_median = median(ubo),
      upper_bound_offdiag_max    = max(ubo),
      ratio_offdiag_median       = median(rto),
      ratio_offdiag_q99          = as.numeric(quantile(rto, 0.99, names = FALSE)),
      ratio_offdiag_max          = max(rto)
    )
  }
  as.list(out)
}

# ------------------------------------------------------------
# 4. Run over every network
# ------------------------------------------------------------
res <- list()

clin_specs <- list(UPDRS_III = "UPDRS_adj", MoCA = "MOCA_adj", LEDD = "LEDD_adj")
for (nm in names(clin_specs)) {
  fn <- file.path(clinical_dir, paste0(tolower(sub("_III", "", nm)), "_affinity_mat.rds"))
  if (!file.exists(fn)) { message("skipping ", nm, ": ", fn, " not found"); next }
  message("summarising ", nm)
  W <- readRDS(fn)
  d <- single_var_dist(pd_model, clin_specs[[nm]], pd_ids)
  res[[nm]] <- summarise_affinity(W, d, nm, n_features = 1)
}

# metabolomics
mfn <- file.path(metabolomics_dir, "metab_affinity_mat.rds")
if (file.exists(mfn)) {
  message("summarising metabolomics")
  res[["Metabolomics"]] <- summarise_affinity(
    readRDS(mfn),
    readRDS(file.path(metabolomics_dir, "metab_dist_snf.rds")),
    "Metabolomics",
    n_features = ncol(readRDS(file.path(metabolomics_dir, "metab_scaled.rds")))
  )
}

# combined clinical network used by the representation-comparison branch
combined_candidates <- list.files(
  paths$clinical_networks,
  pattern = "combined.*affinity.*\\.rds$|clinical_combined.*\\.rds$",
  recursive = TRUE, full.names = TRUE)
if (length(combined_candidates)) {
  message("summarising combined clinical: ", basename(combined_candidates[1]))
  res[["Clinical combined"]] <- summarise_affinity(
    readRDS(combined_candidates[1]), NULL, "Clinical combined", n_features = 3)
} else {
  message("NOTE: no combined clinical network found; Table A8c will lack that column.")
}

# fused four-network matrix
fused_candidates <- list.files(
  dirname(paths$clinical_networks),
  pattern = "fused.*\\.rds$|W_fused.*\\.rds$|snf_.*fused.*\\.rds$",
  recursive = TRUE, full.names = TRUE)
if (length(fused_candidates)) {
  message("summarising fused: ", basename(fused_candidates[1]))
  obj <- readRDS(fused_candidates[1])
  Wf  <- if (is.matrix(obj)) obj else obj$W_fused
  if (!is.null(Wf)) {
    res[["Fused (4-network)"]] <- summarise_affinity(
      Wf, NULL, "Fused", row_normalised = TRUE, n_features = NA)
  }
} else {
  message("NOTE: no fused matrix found; set the path manually if Table A8b needs regenerating.")
}

# ------------------------------------------------------------
# 5. Assemble metrics x networks and write out
# ------------------------------------------------------------
all_metrics <- unique(unlist(lapply(res, names)))
tab <- data.frame(metric = all_metrics, stringsAsFactors = FALSE)
for (nm in names(res)) {
  tab[[nm]] <- sapply(all_metrics, function(m) {
    v <- res[[nm]][[m]]
    if (is.null(v)) NA_real_ else as.numeric(v)
  })
}

write.csv(tab, file.path(qc_dir, "affinity_summary_all_networks.csv"), row.names = FALSE)
saveRDS(tab, file.path(qc_dir, "affinity_summary_all_networks.rds"))

fmt <- function(x) {
  if (is.na(x)) return("")
  if (abs(x) >= 1000 || (abs(x) < 0.001 && x != 0)) formatC(x, format = "e", digits = 2)
  else formatC(x, format = "g", digits = 4)
}
cat("\n")
cat(sprintf("%-28s", "metric"))
for (nm in names(res)) cat(sprintf("%16s", substr(nm, 1, 15)))
cat("\n", strrep("-", 28 + 16 * length(res)), "\n", sep = "")
for (i in seq_len(nrow(tab))) {
  cat(sprintf("%-28s", tab$metric[i]))
  for (nm in names(res)) cat(sprintf("%16s", fmt(tab[[nm]][i])))
  cat("\n")
}
cat("\nWritten to: ", file.path(qc_dir, "affinity_summary_all_networks.csv"), "\n")

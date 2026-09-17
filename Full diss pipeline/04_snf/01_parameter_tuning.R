# 04_snf/01_parameter_tuning.R
#
# Tune SNF hyperparameters K and sigma for the PRIMARY FOUR-NETWORK model:
#   1) UPDRS
#   2) MoCA
#   3) LEDD
#   4) Metabolomics
#
# The tuning score is the eigengap of the fused network's normalized
# graph Laplacian, evaluated separately on train and test splits.
#
# This version is optimized to:
#   - split once
#   - scale once
#   - compute distance matrices once
#   - reuse those distances for every (K, sigma) combination
#   - use SNFtool::SNF directly (no custom fusion reimplementation)

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

if (!requireNamespace("SNFtool", quietly = TRUE)) {
  stop("SNFtool is not installed.")
}

# ------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------
clinical_file <- files$clinical_processed
metabolomics_file <- files$metabolomics_processed

if (!file.exists(clinical_file)) stop("Missing clinical processed file: ", clinical_file)
if (!file.exists(metabolomics_file)) stop("Missing metabolomics processed file: ", metabolomics_file)

clinical_obj <- readRDS(clinical_file)
metabolomics_obj <- readRDS(metabolomics_file)

if (!all(c("pd_model", "sample_ids") %in% names(clinical_obj))) {
  stop("clinical_processed.rds must contain 'pd_model' and 'sample_ids'.")
}
if (!all(c("metab_pd", "sample_ids") %in% names(metabolomics_obj))) {
  stop("metabolomics_processed.rds must contain 'metab_pd' and 'sample_ids'.")
}

pd_model <- as.data.frame(clinical_obj$pd_model)
pd_ids <- as.character(clinical_obj$sample_ids)

metab_pd <- as.matrix(metabolomics_obj$metab_pd)
metab_ids <- as.character(metabolomics_obj$sample_ids)

if (!identical(pd_ids, metab_ids)) {
  stop("Clinical and metabolomics sample_ids do not match.")
}
if (anyDuplicated(pd_ids)) {
  stop("Duplicated sample IDs found in the primary cohort.")
}
if (!identical(as.character(pd_model$Anonymised_sampleID), pd_ids)) {
  stop("Clinical row order does not match sample_ids.")
}
if (anyNA(pd_model[, c("UPDRS_adj", "MOCA_adj", "LEDD_adj")])) {
  stop("Clinical tuning variables contain NA values.")
}
if (anyNA(metab_pd)) {
  stop("metab_pd contains NA values.")
}

required_cols <- c("Anonymised_sampleID", "UPDRS_adj", "MOCA_adj", "LEDD_adj")
missing_cols <- setdiff(required_cols, names(pd_model))
if (length(missing_cols) > 0) {
  stop("Missing required clinical columns: ", paste(missing_cols, collapse = ", "))
}

# ------------------------------------------------------------
# 2. Output folder
# ------------------------------------------------------------
tuning_dir <- file.path(paths$snf, "parameter_tuning", "primary_four_network")
dir.create(tuning_dir, recursive = TRUE, showWarnings = FALSE)

results_rds <- file.path(tuning_dir, "snf_tuning_results_primary_four_network.rds")
results_csv <- file.path(tuning_dir, "snf_tuning_results_primary_four_network.csv")
best_train_csv <- file.path(tuning_dir, "snf_best_by_train_gap.csv")
best_test_csv <- file.path(tuning_dir, "snf_best_by_test_gap.csv")
recommended_csv <- file.path(tuning_dir, "recommended_parameters.csv")
split_ids_rds <- file.path(tuning_dir, "train_test_split_ids.rds")

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------
fit_scaler <- function(x) {
  x <- as.matrix(x)
  center <- colMeans(x, na.rm = TRUE)
  scale_vec <- apply(x, 2, sd, na.rm = TRUE)
  scale_vec[!is.finite(scale_vec) | scale_vec == 0] <- 1
  list(center = center, scale = scale_vec)
}

apply_scaler <- function(x, scaler) {
  x <- as.matrix(x)
  x_centered <- sweep(x, 2, scaler$center, "-")
  x_scaled <- sweep(x_centered, 2, scaler$scale, "/")
  as.matrix(x_scaled)
}

make_affinity_from_dist <- function(dist_mat, ids, K, sigma) {
  dist_mat <- as.matrix(dist_mat)
  dimnames(dist_mat) <- NULL

  K_eff <- max(1, min(K, nrow(dist_mat) - 1))
  aff <- SNFtool::affinityMatrix(dist_mat, K = K_eff, sigma = sigma)
  aff <- as.matrix(aff)

  rownames(aff) <- ids
  colnames(aff) <- ids
  aff
}

summarise_network <- function(W) {
  W <- as.matrix(W)

  if (nrow(W) != ncol(W)) stop("Fused network is not square.")
  if (anyNA(W)) stop("Fused network contains NA values.")
  if (!isTRUE(all.equal(W, t(W)))) stop("Fused network is not symmetric.")

  offdiag <- W[upper.tri(W)]

  deg <- rowSums(W)
  deg[deg <= 0] <- .Machine$double.eps
  d_inv_sqrt <- diag(1 / sqrt(deg))
  laplacian <- diag(nrow(W)) - d_inv_sqrt %*% W %*% d_inv_sqrt
  eigvals <- sort(eigen(laplacian, symmetric = TRUE, only.values = TRUE)$values)

  gap_limit <- min(10, length(eigvals) - 1)
  if (gap_limit < 1) stop("Not enough eigenvalues to compute an eigengap.")
  gaps <- diff(eigvals[1:(gap_limit + 1)])

  list(
    range = range(W),
    offdiag_summary = summary(as.vector(offdiag)),
    offdiag_quantiles = quantile(
      offdiag,
      probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1),
      na.rm = TRUE
    ),
    n_above_001 = sum(offdiag > 0.01, na.rm = TRUE),
    prop_above_001 = mean(offdiag > 0.01, na.rm = TRUE),
    best_k = which.max(gaps) + 1L,
    best_gap = max(gaps),
    eigvals = eigvals[1:min(10, length(eigvals))]
  )
}

flatten_result <- function(x) {
  data.frame(
    K = x$K,
    sigma = x$sigma,
    K_eff = x$K_eff,
    T = x$T,
    error = ifelse(is.na(x$error), "", x$error),

    train_range_min = if (is.null(x$train)) NA_real_ else x$train$range[1],
    train_range_max = if (is.null(x$train)) NA_real_ else x$train$range[2],
    train_offdiag_median = if (is.null(x$train)) NA_real_ else as.numeric(x$train$offdiag_summary["Median"]),
    train_offdiag_mean = if (is.null(x$train)) NA_real_ else as.numeric(x$train$offdiag_summary["Mean"]),
    train_offdiag_q95 = if (is.null(x$train)) NA_real_ else as.numeric(x$train$offdiag_quantiles["95%"]),
    train_n_above_001 = if (is.null(x$train)) NA_real_ else x$train$n_above_001,
    train_prop_above_001 = if (is.null(x$train)) NA_real_ else x$train$prop_above_001,
    train_best_k = if (is.null(x$train)) NA_real_ else x$train$best_k,
    train_best_gap = if (is.null(x$train)) NA_real_ else x$train$best_gap,

    test_range_min = if (is.null(x$test)) NA_real_ else x$test$range[1],
    test_range_max = if (is.null(x$test)) NA_real_ else x$test$range[2],
    test_offdiag_median = if (is.null(x$test)) NA_real_ else as.numeric(x$test$offdiag_summary["Median"]),
    test_offdiag_mean = if (is.null(x$test)) NA_real_ else as.numeric(x$test$offdiag_summary["Mean"]),
    test_offdiag_q95 = if (is.null(x$test)) NA_real_ else as.numeric(x$test$offdiag_quantiles["95%"]),
    test_n_above_001 = if (is.null(x$test)) NA_real_ else x$test$n_above_001,
    test_prop_above_001 = if (is.null(x$test)) NA_real_ else x$test$prop_above_001,
    test_best_k = if (is.null(x$test)) NA_real_ else x$test$best_k,
    test_best_gap = if (is.null(x$test)) NA_real_ else x$test$best_gap,

    k_abs_diff = if (is.null(x$train) || is.null(x$test)) NA_real_ else abs(x$train$best_k - x$test$best_k),
    gap_abs_diff = if (is.null(x$train) || is.null(x$test)) NA_real_ else abs(x$train$best_gap - x$test$best_gap),
    consistency_score = if (is.null(x$train) || is.null(x$test)) NA_real_ else pmin(x$train$best_gap, x$test$best_gap),
    stringsAsFactors = FALSE
  )
}

# SNFtool::SNF usually uses these data to identify clusters. It expects
# the matrices to be similarity/affinity matrices already.
run_fusion <- function(aff_list, K_eff, T) {
  fused <- SNFtool::SNF(aff_list, K = K_eff, t = T)
  fused <- as.matrix(fused)
  fused
}

# ------------------------------------------------------------
# 4. Fixed 80/20 split
# ------------------------------------------------------------
set.seed(20240622)
train_fraction <- 0.80
n_total <- length(pd_ids)
n_train <- floor(train_fraction * n_total)

train_ids <- sample(pd_ids, size = n_train, replace = FALSE)
train_ids <- pd_ids[pd_ids %in% train_ids]
test_ids <- pd_ids[!pd_ids %in% train_ids]

if (length(train_ids) == 0 || length(test_ids) == 0) stop("Train/test split failed.")
if (length(intersect(train_ids, test_ids)) != 0) stop("Train and test IDs overlap.")
if (length(train_ids) + length(test_ids) != length(pd_ids)) stop("Train/test split does not cover all patients.")

saveRDS(
  list(
    train_ids = train_ids,
    test_ids = test_ids,
    seed = 20240622,
    train_fraction = train_fraction
  ),
  split_ids_rds
)

message("Total samples: ", n_total)
message("Training samples: ", length(train_ids))
message("Test samples: ", length(test_ids))

# ------------------------------------------------------------
# 5. Raw matrices and scaling
# ------------------------------------------------------------
clinical_raw <- as.matrix(pd_model[, c("UPDRS_adj", "MOCA_adj", "LEDD_adj"), drop = FALSE])
rownames(clinical_raw) <- pd_ids
metab_raw <- as.matrix(metab_pd)
rownames(metab_raw) <- pd_ids

clinical_train_raw <- clinical_raw[train_ids, , drop = FALSE]
clinical_test_raw <- clinical_raw[test_ids, , drop = FALSE]
metab_train_raw <- metab_raw[train_ids, , drop = FALSE]
metab_test_raw <- metab_raw[test_ids, , drop = FALSE]

saveRDS(clinical_train_raw, file.path(tuning_dir, "clinical_train_raw.rds"))
saveRDS(clinical_test_raw, file.path(tuning_dir, "clinical_test_raw.rds"))
saveRDS(metab_train_raw, file.path(tuning_dir, "metab_train_raw.rds"))
saveRDS(metab_test_raw, file.path(tuning_dir, "metab_test_raw.rds"))

clinical_scaler <- fit_scaler(clinical_train_raw)
metab_scaler <- fit_scaler(metab_train_raw)

saveRDS(clinical_scaler, file.path(tuning_dir, "clinical_scaler_train.rds"))
saveRDS(metab_scaler, file.path(tuning_dir, "metab_scaler_train.rds"))

clinical_train_scaled <- apply_scaler(clinical_train_raw, clinical_scaler)
clinical_test_scaled <- apply_scaler(clinical_test_raw, clinical_scaler)
metab_train_scaled <- apply_scaler(metab_train_raw, metab_scaler)
metab_test_scaled <- apply_scaler(metab_test_raw, metab_scaler)

saveRDS(clinical_train_scaled, file.path(tuning_dir, "clinical_train_scaled.rds"))
saveRDS(clinical_test_scaled, file.path(tuning_dir, "clinical_test_scaled.rds"))
saveRDS(metab_train_scaled, file.path(tuning_dir, "metab_train_scaled.rds"))
saveRDS(metab_test_scaled, file.path(tuning_dir, "metab_test_scaled.rds"))

# ------------------------------------------------------------
# 6. Distances computed once
# ------------------------------------------------------------
clinical_updrs_train_dist <- SNFtool::dist2(clinical_train_scaled[, "UPDRS_adj", drop = FALSE], clinical_train_scaled[, "UPDRS_adj", drop = FALSE])
clinical_moca_train_dist  <- SNFtool::dist2(clinical_train_scaled[, "MOCA_adj", drop = FALSE], clinical_train_scaled[, "MOCA_adj", drop = FALSE])
clinical_ledd_train_dist  <- SNFtool::dist2(clinical_train_scaled[, "LEDD_adj", drop = FALSE], clinical_train_scaled[, "LEDD_adj", drop = FALSE])

clinical_updrs_test_dist <- SNFtool::dist2(clinical_test_scaled[, "UPDRS_adj", drop = FALSE], clinical_test_scaled[, "UPDRS_adj", drop = FALSE])
clinical_moca_test_dist  <- SNFtool::dist2(clinical_test_scaled[, "MOCA_adj", drop = FALSE], clinical_test_scaled[, "MOCA_adj", drop = FALSE])
clinical_ledd_test_dist  <- SNFtool::dist2(clinical_test_scaled[, "LEDD_adj", drop = FALSE], clinical_test_scaled[, "LEDD_adj", drop = FALSE])

metab_train_dist <- SNFtool::dist2(metab_train_scaled, metab_train_scaled)
metab_test_dist  <- SNFtool::dist2(metab_test_scaled, metab_test_scaled)

saveRDS(clinical_updrs_train_dist, file.path(tuning_dir, "clinical_updrs_train_dist.rds"))
saveRDS(clinical_moca_train_dist, file.path(tuning_dir, "clinical_moca_train_dist.rds"))
saveRDS(clinical_ledd_train_dist, file.path(tuning_dir, "clinical_ledd_train_dist.rds"))
saveRDS(clinical_updrs_test_dist, file.path(tuning_dir, "clinical_updrs_test_dist.rds"))
saveRDS(clinical_moca_test_dist, file.path(tuning_dir, "clinical_moca_test_dist.rds"))
saveRDS(clinical_ledd_test_dist, file.path(tuning_dir, "clinical_ledd_test_dist.rds"))
saveRDS(metab_train_dist, file.path(tuning_dir, "metab_train_dist.rds"))
saveRDS(metab_test_dist, file.path(tuning_dir, "metab_test_dist.rds"))

# ------------------------------------------------------------
# 7. Grid search
# ------------------------------------------------------------
K_grid <- c(10, 15, 20, 25, 30)
sigma_grid <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

all_combos <- expand.grid(
  K = K_grid,
  sigma = sigma_grid,
  stringsAsFactors = FALSE
)

expected_cols <- names(flatten_result(list(
  K = 1, sigma = 0.3, K_eff = 1, T = snf_params$t,
  error = NA_character_, train = NULL, test = NULL
)))

if (file.exists(results_rds)) {
  old_results <- readRDS(results_rds)
  if (identical(names(old_results), expected_cols)) {
    results_df <- old_results
  } else {
    message("Existing results file has a different schema; starting a fresh results table.")
    results_df <- data.frame()
  }
} else {
  results_df <- data.frame()
}

completed <- character(0)
if (nrow(results_df) > 0) {
  completed <- paste(results_df$K, results_df$sigma, sep = "_")
}

make_split_affinities <- function(updrs_dist, moca_dist, ledd_dist, metab_dist, ids, K, sigma) {
  K_eff <- max(
    1L,
    min(
      K,
      nrow(updrs_dist) - 1L,
      nrow(moca_dist) - 1L,
      nrow(ledd_dist) - 1L,
      nrow(metab_dist) - 1L
    )
  )

  list(
    K_eff = K_eff,
    updrs = make_affinity_from_dist(updrs_dist, ids, K_eff, sigma),
    moca  = make_affinity_from_dist(moca_dist, ids, K_eff, sigma),
    ledd  = make_affinity_from_dist(ledd_dist, ids, K_eff, sigma),
    metab = make_affinity_from_dist(metab_dist, ids, K_eff, sigma)
  )
}

run_one_setting <- function(K, sigma) {
  out <- list(
    K = K,
    sigma = sigma,
    K_eff = NA_integer_,
    T = snf_params$t,
    train = NULL,
    test = NULL,
    error = NA_character_
  )

  out <- tryCatch({
    train_aff <- make_split_affinities(
      clinical_updrs_train_dist,
      clinical_moca_train_dist,
      clinical_ledd_train_dist,
      metab_train_dist,
      train_ids,
      K,
      sigma
    )

    test_aff <- make_split_affinities(
      clinical_updrs_test_dist,
      clinical_moca_test_dist,
      clinical_ledd_test_dist,
      metab_test_dist,
      test_ids,
      K,
      sigma
    )

    out$K_eff <- min(train_aff$K_eff, test_aff$K_eff)

    fused_train <- run_fusion(
      list(train_aff$updrs, train_aff$moca, train_aff$ledd, train_aff$metab),
      K_eff = out$K_eff,
      T = snf_params$t
    )
    fused_test <- run_fusion(
      list(test_aff$updrs, test_aff$moca, test_aff$ledd, test_aff$metab),
      K_eff = out$K_eff,
      T = snf_params$t
    )

    rownames(fused_train) <- train_ids
    colnames(fused_train) <- train_ids
    rownames(fused_test) <- test_ids
    colnames(fused_test) <- test_ids

    out$train <- summarise_network(fused_train)
    out$test <- summarise_network(fused_test)

    out
  }, error = function(e) {
    out$error <- conditionMessage(e)
    out
  })

  out
}

for (i in seq_len(nrow(all_combos))) {
  K <- all_combos$K[i]
  sigma <- all_combos$sigma[i]
  combo_id <- paste(K, sigma, sep = "_")

  if (combo_id %in% completed) {
    message("Skipping already completed K = ", K, ", sigma = ", sigma)
    next
  }

  message("Running K = ", K, ", sigma = ", sigma)
  one_result <- run_one_setting(K, sigma)
  one_row <- flatten_result(one_result)


  if (nrow(results_df) == 0) {
    results_df <- one_row
  } else {
    one_row <- one_row[, names(results_df), drop = FALSE]
    results_df <- rbind(results_df, one_row)
  }

  saveRDS(results_df, results_rds)
  write.csv(results_df, results_csv, row.names = FALSE)
  saveRDS(one_result, file.path(tuning_dir, paste0("checkpoint_K", K, "_sigma", sigma, ".rds")))
}

# ------------------------------------------------------------
# 8. Best settings
# ------------------------------------------------------------
if (nrow(results_df) == 0) {
  stop("No tuning results were produced.")
}

valid_train <- results_df[is.finite(results_df$train_best_gap), , drop = FALSE]
valid_test  <- results_df[is.finite(results_df$test_best_gap), , drop = FALSE]
valid_both  <- results_df[
  is.finite(results_df$train_best_gap) & is.finite(results_df$test_best_gap),
  ,
  drop = FALSE
]

if (nrow(valid_train) == 0) stop("No valid training eigengap results found.")
if (nrow(valid_test) == 0) stop("No valid test eigengap results found.")
if (nrow(valid_both) == 0) stop("No rows have valid train and test eigengaps.")

best_train_row <- valid_train[which.max(valid_train$train_best_gap), , drop = FALSE]
best_test_row <- valid_test[which.max(valid_test$test_best_gap), , drop = FALSE]
best_recommended_row <- valid_both[which.max(valid_both$consistency_score), , drop = FALSE]

write.csv(best_train_row, best_train_csv, row.names = FALSE)
write.csv(best_test_row, best_test_csv, row.names = FALSE)
write.csv(best_recommended_row, recommended_csv, row.names = FALSE)

cat("\nFinished SNF tuning grid.\n")
cat("Results saved to: ", tuning_dir, "\n", sep = "")

cat("\nBest training-gap setting:\n")
print(best_train_row)

cat("\nBest test-gap setting:\n")
print(best_test_row)

cat("\nRecommended setting by consistency_score = min(train_gap, test_gap):\n")
print(best_recommended_row)
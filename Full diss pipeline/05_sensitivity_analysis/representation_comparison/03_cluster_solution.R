# =========================================================
# 03_cluster_solution.R
# Representation comparison branch clustering
# Re-runnable from a fresh R session.
#
# Methodological notes:
# - Reads the saved 2-network fused network produced by 02_run_2network_snf.R.
# - Computes the eigengap summary here.
# - Stores the eigengap summary, including the second-best k for reference.
# - Produces cluster labels for the best k only.
# - Does not change the fused network itself.
# =========================================================

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

if (!requireNamespace("SNFtool", quietly = TRUE)) {
  stop("SNFtool is required but is not available.")
}

# -------------------------------------------------------------------
# Input / output locations
# -------------------------------------------------------------------
branch_dir <- file.path(paths$sensitivity, "representation_comparison")
snf_file_candidates <- c(
  file.path(branch_dir, "representation_2network_snf_results.rds"),
  file.path(branch_dir, "representation_2network_fused_network.rds")
)

snf_file <- snf_file_candidates[file.exists(snf_file_candidates)][1]
if (is.na(snf_file) || !nzchar(snf_file)) {
  stop(
    "Missing 2-network SNF results file. Expected one of:\n",
    paste(snf_file_candidates, collapse = "\n")
  )
}

out_dir <- file.path(branch_dir, "cluster_solution")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------
read_square_matrix <- function(path, object_name) {
  mat <- readRDS(path)
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop(object_name, " is empty: ", path)
  }
  if (nrow(mat) != ncol(mat)) {
    stop(object_name, " must be square: ", path)
  }
  if (anyNA(mat)) {
    stop(object_name, " contains NA values: ", path)
  }
  if (!isTRUE(all.equal(mat, t(mat)))) {
    stop(object_name, " is not symmetric: ", path)
  }
  mat
}

compute_laplacian_eigendecomposition <- function(W) {
  deg <- rowSums(W)
  deg[deg <= 0] <- .Machine$double.eps
  d_inv_sqrt <- diag(1 / sqrt(deg))
  laplacian <- diag(nrow(W)) - d_inv_sqrt %*% W %*% d_inv_sqrt

  eig <- eigen(laplacian, symmetric = TRUE)
  ord <- order(eig$values)

  list(
    laplacian = laplacian,
    eigvals = eig$values[ord],
    eigvecs = eig$vectors[, ord, drop = FALSE]
  )
}

compute_eigengap_table <- function(eigvals, max_gaps = 10L) {
  if (length(eigvals) < 2) {
    stop("Not enough eigenvalues to compute eigengaps.")
  }

  gap_limit <- min(max_gaps, length(eigvals) - 1L)
  gaps <- diff(eigvals[seq_len(gap_limit + 1L)])

  gap_table <- data.frame(
    gap_rank = seq_along(gaps),
    cluster_k = seq_along(gaps) + 1L,
    gap_value = gaps,
    stringsAsFactors = FALSE
  )

  gap_table <- gap_table[order(gap_table$gap_value, decreasing = TRUE), , drop = FALSE]
  rownames(gap_table) <- NULL

  if (nrow(gap_table) == 0) {
    stop("No eigengaps were computed.")
  }

  best_row <- gap_table[1, , drop = FALSE]
  second_row <- if (nrow(gap_table) >= 2) gap_table[2, , drop = FALSE] else NULL

  list(
    gap_table = gap_table,
    best_k = as.integer(best_row$cluster_k[1]),
    best_gap = as.numeric(best_row$gap_value[1]),
    second_best_k = if (!is.null(second_row)) as.integer(second_row$cluster_k[1]) else NA_integer_,
    second_best_gap = if (!is.null(second_row)) as.numeric(second_row$gap_value[1]) else NA_real_
  )
}

run_spectral_clustering <- function(W, k) {
  if (is.na(k) || k < 2L) {
    stop("k must be at least 2.")
  }
  labels <- SNFtool::spectralClustering(as.matrix(W), K = k)
  as.integer(labels)
}

cluster_size_table <- function(labels_df) {
  tab <- table(labels_df$cluster)
  data.frame(
    cluster = as.integer(names(tab)),
    n_samples = as.integer(tab),
    proportion = as.numeric(tab) / sum(tab),
    stringsAsFactors = FALSE
  )
}

write_cluster_outputs <- function(labels_df, prefix, out_dir) {
  write.csv(
    labels_df,
    file.path(out_dir, paste0(prefix, "_cluster_labels.csv")),
    row.names = FALSE
  )

  size_df <- cluster_size_table(labels_df)
  write.csv(
    size_df,
    file.path(out_dir, paste0(prefix, "_cluster_sizes.csv")),
    row.names = FALSE
  )
}

# -------------------------------------------------------------------
# Load fused network and sample IDs
# -------------------------------------------------------------------
snf_obj <- readRDS(snf_file)

if (!all(c("fused_network", "sample_ids") %in% names(snf_obj))) {
  stop("2-network SNF results must contain 'fused_network' and 'sample_ids'.")
}

fused_network <- as.matrix(snf_obj$fused_network)
sample_ids <- as.character(snf_obj$sample_ids)

if (is.null(rownames(fused_network)) || is.null(colnames(fused_network))) {
  rownames(fused_network) <- sample_ids
  colnames(fused_network) <- sample_ids
}

if (nrow(fused_network) != length(sample_ids)) {
  stop("Fused network dimensions do not match sample_ids length.")
}
rownames(fused_network) <- sample_ids
colnames(fused_network) <- sample_ids

if (!identical(rownames(fused_network), sample_ids)) {
  stop("Row names of fused network do not match sample IDs.")
}
if (!identical(colnames(fused_network), sample_ids)) {
  stop("Column names of fused network do not match sample IDs.")
}
if (anyNA(fused_network)) {
  stop("Fused network contains NA values.")
}
if (!isTRUE(all.equal(fused_network, t(fused_network)))) {
  stop("Fused network is not symmetric.")
}

message("Loaded fused network with ", nrow(fused_network), " samples.")

# -------------------------------------------------------------------
# Eigengap analysis
# -------------------------------------------------------------------
message("Computing Laplacian eigendecomposition...")

eig_obj <- compute_laplacian_eigendecomposition(fused_network)
eigvals <- eig_obj$eigvals
eigvecs <- eig_obj$eigvecs

eigengap_obj <- compute_eigengap_table(eigvals, max_gaps = 10L)
gap_table <- eigengap_obj$gap_table

selected_summary <- data.frame(
  metric = c(
    "best_k",
    "best_gap",
    "second_best_k",
    "second_best_gap",
    "gap_difference",
    "gap_ratio_second_to_best"
  ),
  value = c(
    eigengap_obj$best_k,
    eigengap_obj$best_gap,
    eigengap_obj$second_best_k,
    eigengap_obj$second_best_gap,
    eigengap_obj$best_gap - eigengap_obj$second_best_gap,
    eigengap_obj$second_best_gap / eigengap_obj$best_gap
  ),
  stringsAsFactors = FALSE
)

write.csv(
  gap_table,
  file.path(out_dir, "eigengap_candidates.csv"),
  row.names = FALSE
)

write.csv(
  selected_summary,
  file.path(out_dir, "eigengap_selected_summary.csv"),
  row.names = FALSE
)

saveRDS(gap_table, file.path(out_dir, "eigengap_candidates.rds"))
saveRDS(selected_summary, file.path(out_dir, "eigengap_selected_summary.rds"))

png(file.path(out_dir, "eigengap_plot.png"), width = 2000, height = 1500, res = 300)
plot(
  gap_table$cluster_k,
  gap_table$gap_value,
  type = "b",
  pch = 16,
  xlab = "Candidate number of clusters (k)",
  ylab = "Eigengap",
  main = "Eigengap candidates for the 2-network fused network"
)
abline(v = eigengap_obj$best_k, lty = 2)
dev.off()

# -------------------------------------------------------------------
# Spectral clustering for best and second-best k
# -------------------------------------------------------------------
message("Running spectral clustering for best k = ", eigengap_obj$best_k, "...")

best_labels <- run_spectral_clustering(fused_network, eigengap_obj$best_k)

if (length(best_labels) != length(sample_ids)) {
  stop("Spectral clustering returned the wrong number of labels.")
}

best_labels_df <- data.frame(
  sample_id = sample_ids,
  cluster   = best_labels,
  stringsAsFactors = FALSE
)
best_labels_df <- best_labels_df[order(best_labels_df$sample_id), , drop = FALSE]
rownames(best_labels_df) <- NULL

write_cluster_outputs(best_labels_df, "best_k", out_dir)

# -------------------------------------------------------------------
# Save a compact run record
# -------------------------------------------------------------------
run_record <- data.frame(
  metric = c(
    "n_samples",
    "best_k",
    "best_gap",
    "second_best_k",
    "second_best_gap"
  ),
  value = c(
    nrow(fused_network),
    eigengap_obj$best_k,
    eigengap_obj$best_gap,
    eigengap_obj$second_best_k,
    eigengap_obj$second_best_gap
  ),
  stringsAsFactors = FALSE
)

write.csv(run_record, file.path(out_dir, "cluster_run_record.csv"), row.names = FALSE)
saveRDS(run_record, file.path(out_dir, "cluster_run_record.rds"))

# Optional compact object mirroring the main pipeline style
cluster_results <- list(
  model_name = "representation_2network",
  sample_ids = sample_ids,
  fused_network = fused_network,
  eigengap_candidates = gap_table,
  eigengap_selected_summary = selected_summary,
  best_k = eigengap_obj$best_k,
  second_best_k = eigengap_obj$second_best_k,
  best_k_cluster_labels = best_labels_df
)


saveRDS(cluster_results, file.path(out_dir, "representation_2network_cluster_results.rds"))

message("Cluster solution complete.")
message("Outputs saved to: ", out_dir)
message("Selected best k: ", eigengap_obj$best_k)
message("Selected second-best k: ", eigengap_obj$second_best_k)
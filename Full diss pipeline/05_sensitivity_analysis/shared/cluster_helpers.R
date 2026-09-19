# Shared helper functions for the sensitivity-analysis clustering scripts.
# These functions perform eigengap scans, spectral clustering, cluster summaries,
# and saving of cluster-analysis outputs.
#
# The eigengap scan evaluates candidate cluster numbers up to `max_k` and
# records the best and second-best eigengaps.
#
# No project-specific file paths or parameters need to be set in this file.

compute_eigengap_scan <- function(W, max_k = 10L) {
  W <- as.matrix(W)

  if (is.null(rownames(W)) || is.null(colnames(W))) {
    stop("Fused network must have row and column names.")
  }
  if (!isTRUE(all.equal(W, t(W)))) {
    stop("Fused network is not symmetric.")
  }
  if (anyNA(W)) {
    stop("Fused network contains NA values.")
  }

  deg <- rowSums(W)
  deg[deg <= 0] <- .Machine$double.eps
  D_inv_sqrt <- diag(1 / sqrt(deg))
  L <- diag(nrow(W)) - D_inv_sqrt %*% W %*% D_inv_sqrt

  eigvals <- sort(eigen(L, symmetric = TRUE, only.values = TRUE)$values)
  max_k <- min(max_k, length(eigvals) - 1L)

  if (max_k < 2L) {
    stop("Not enough eigenvalues available for an eigengap scan.")
  }

  gaps <- diff(eigvals[seq_len(max_k + 1L)])
  gap_table <- data.frame(
    k = seq_len(max_k) + 1L,
    eigengap = gaps,
    stringsAsFactors = FALSE
  )

  ranked_idx <- order(gaps, decreasing = TRUE)
  best_pos <- ranked_idx[1L]
  second_pos <- if (length(ranked_idx) >= 2L) ranked_idx[2L] else NA_integer_

  best_k <- best_pos + 1L
  second_best_k <- if (is.na(second_pos)) NA_integer_ else second_pos + 1L

  list(
    laplacian = L,
    eigvals = eigvals,
    gap_table = gap_table,
    best_k = best_k,
    second_best_k = second_best_k,
    best_gap = gaps[best_pos],
    second_best_gap = if (is.na(second_pos)) NA_real_ else gaps[second_pos]
  )
}

run_spectral_clustering <- function(W, k) {
  if (is.na(k) || k < 2L) {
    stop("k must be at least 2.")
  }

  labels <- SNFtool::spectralClustering(as.matrix(W), K = k)
  as.integer(labels)
}

summarise_cluster_solution <- function(
    sample_ids,
    cluster_labels,
    scan,
    chosen_k,
    model_name
) {
  cluster_df <- data.frame(
    sample_id = sample_ids,
    cluster = as.integer(cluster_labels),
    stringsAsFactors = FALSE
  )

  cluster_sizes <- as.data.frame(table(cluster_df$cluster), stringsAsFactors = FALSE)
  names(cluster_sizes) <- c("cluster", "n")

  data.frame(
    model_name = model_name,
    n_samples = length(sample_ids),
    chosen_k = chosen_k,
    best_k = scan$best_k,
    second_best_k = scan$second_best_k,
    best_gap = scan$best_gap,
    second_best_gap = scan$second_best_gap,
    n_clusters_observed = nrow(cluster_sizes),
    stringsAsFactors = FALSE
  )
}

save_cluster_outputs <- function(
    out_dir,
    prefix,
    sample_ids,
    fused_network,
    scan,
    chosen_k,
    cluster_labels,
    model_name
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cluster_df <- data.frame(
    sample_id = sample_ids,
    cluster = as.integer(cluster_labels),
    stringsAsFactors = FALSE
  )
  cluster_sizes <- as.data.frame(table(cluster_df$cluster), stringsAsFactors = FALSE)
  names(cluster_sizes) <- c("cluster", "n")

  summary_table <- summarise_cluster_solution(
    sample_ids = sample_ids,
    cluster_labels = cluster_labels,
    scan = scan,
    chosen_k = chosen_k,
    model_name = model_name
  )

  saveRDS(
    list(
      model_name = model_name,
      sample_ids = sample_ids,
      fused_network = fused_network,
      eigengap_scan = scan$gap_table,
      best_k = scan$best_k,
      second_best_k = scan$second_best_k,
      best_gap = scan$best_gap,
      second_best_gap = scan$second_best_gap,
      chosen_k = chosen_k,
      cluster_labels = cluster_df,
      cluster_sizes = cluster_sizes,
      summary_table = summary_table
    ),
    file.path(out_dir, paste0(prefix, "_cluster_results.rds"))
  )

  write.csv(cluster_df, file.path(out_dir, paste0(prefix, "_cluster_labels.csv")), row.names = FALSE)
  write.csv(cluster_sizes, file.path(out_dir, paste0(prefix, "_cluster_sizes.csv")), row.names = FALSE)
  write.csv(summary_table, file.path(out_dir, paste0(prefix, "_cluster_summary.csv")), row.names = FALSE)
  write.csv(scan$gap_table, file.path(out_dir, paste0(prefix, "_eigengap_scan.csv")), row.names = FALSE)

  png(file.path(out_dir, paste0(prefix, "_eigengap_plot.png")), width = 2000, height = 1500, res = 300)
  plot(
    scan$gap_table$k,
    scan$gap_table$eigengap,
    type = "b",
    pch = 16,
    xlab = "k",
    ylab = "Eigengap",
    main = paste0(prefix, " eigengap scan")
  )
  abline(v = scan$best_k, lty = 2)
  if (!is.na(scan$second_best_k)) abline(v = scan$second_best_k, lty = 3)
  dev.off()

  png(file.path(out_dir, paste0(prefix, "_cluster_sizes.png")), width = 2000, height = 1500, res = 300)
  barplot(
    cluster_sizes$n,
    names.arg = cluster_sizes$cluster,
    xlab = "Cluster",
    ylab = "Patients",
    main = paste0(prefix, " cluster sizes")
  )
  dev.off()

  invisible(
    list(
      cluster_df = cluster_df,
      cluster_sizes = cluster_sizes,
      summary_table = summary_table
    )
  )
}

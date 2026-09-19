# Shared helper functions for running SNF in the sensitivity analyses.
#
# These functions validate affinity matrices, check that networks have matching
# sample IDs, run SNF on a list of affinity matrices, summarise the fused
# network, and save the fused network outputs and QC summary.
#
# No project-specific file paths or parameters need to be set in this file.

check_affinity_matrix <- function(W, name = "matrix") {
  W <- as.matrix(W)

  if (is.null(W) || length(dim(W)) != 2L) {
    stop(name, " must be a 2D matrix.")
  }
  if (nrow(W) != ncol(W)) {
    stop(name, " must be square.")
  }
  if (anyNA(W)) {
    stop(name, " contains NA values.")
  }
  if (!isTRUE(all.equal(W, t(W)))) {
    stop(name, " is not symmetric.")
  }
  if (is.null(rownames(W)) || is.null(colnames(W))) {
    stop(name, " must have row and column names.")
  }
  if (!identical(rownames(W), colnames(W))) {
    stop(name, " row and column names do not match.")
  }

  invisible(TRUE)
}

check_network_list <- function(network_list) {
  if (!is.list(network_list) || length(network_list) < 2L) {
    stop("network_list must be a list containing at least two affinity matrices.")
  }

  for (i in seq_along(network_list)) {
    check_affinity_matrix(network_list[[i]], paste0("network_", i))
  }

  sample_ids <- rownames(network_list[[1]])
  for (i in seq_along(network_list)[-1]) {
    if (!identical(sample_ids, rownames(network_list[[i]]))) {
      stop("Network sample order mismatch at position ", i, ".")
    }
  }

  invisible(sample_ids)
}

run_snf_from_affinities <- function(network_list, K = 30, T = 20) {
  sample_ids <- check_network_list(network_list)

  fused_network <- SNFtool::SNF(network_list, K = K, t = T)
  fused_network <- as.matrix(fused_network)

  rownames(fused_network) <- sample_ids
  colnames(fused_network) <- sample_ids

  if (!isTRUE(all.equal(fused_network, t(fused_network)))) {
    stop("Fused network is not symmetric.")
  }
  if (anyNA(fused_network)) {
    stop("Fused network contains NA values.")
  }

  fused_network
}

summarise_fused_network <- function(fused_network, K, T, network_names) {
  offdiag <- fused_network[upper.tri(fused_network)]

  data.frame(
    metric = c(
      "n_samples",
      "n_networks",
      "K",
      "T",
      "min_similarity",
      "median_similarity",
      "mean_similarity",
      "max_similarity"
    ),
    value = c(
      nrow(fused_network),
      length(network_names),
      K,
      T,
      min(offdiag, na.rm = TRUE),
      median(offdiag, na.rm = TRUE),
      mean(offdiag, na.rm = TRUE),
      max(offdiag, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

save_fused_network_outputs <- function(
    fused_network,
    sample_ids,
    out_dir,
    prefix,
    qc_table
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(fused_network, file.path(out_dir, paste0(prefix, "_fused_network.rds")))
  saveRDS(sample_ids, file.path(out_dir, paste0(prefix, "_sample_ids.rds")))
  saveRDS(qc_table, file.path(out_dir, paste0(prefix, "_qc.rds")))
  write.csv(qc_table, file.path(out_dir, paste0(prefix, "_qc.csv")), row.names = FALSE)

  offdiag <- fused_network[upper.tri(fused_network)]
  png(
    file.path(out_dir, paste0(prefix, "_fused_network_histogram.png")),
    width = 2000,
    height = 1500,
    res = 300
  )
  hist(
    offdiag,
    breaks = 50,
    main = paste0(prefix, " fused network similarity distribution"),
    xlab = "Similarity"
  )
  dev.off()

  invisible(TRUE)
}

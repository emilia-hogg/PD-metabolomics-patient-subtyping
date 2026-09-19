# Shared helper functions for the top-20 metabolite sensitivity analyses.
#
# These functions support affinity-matrix construction, SNF, eigengap and
# spectral clustering, comparison of clustering solutions, and saving of
# analysis outputs. The file also defines the top-20 metabolite list and
# functions for matching these metabolites to the metabolomics data.
#
# No project-specific file paths need to be set in this file.

if (!requireNamespace('SNFtool', quietly = TRUE)) {
  stop('SNFtool is required but is not available.')
}

normalize_name <- function(x) {
  x <- as.character(x)
  x <- trimws(tolower(x))
  x <- gsub('[*]+', '', x)
  x <- gsub('[[:punct:]]+', ' ', x)
  x <- gsub('[[:space:]]+', ' ', x)
  trimws(x)
}

make_slug <- function(x) {
  x <- as.character(x)
  x <- tolower(x)
  x <- gsub('[*]+', '', x)
  x <- gsub('\\s+', '_', x)
  x <- gsub('[^a-z0-9_]+', '', x)
  x <- gsub('_+', '_', x)
  x <- gsub('^_|_$', '', x)
  x
}

check_affinity_matrix <- function(mat, name = 'affinity_matrix', ids = NULL) {
  mat <- as.matrix(mat)
  if (nrow(mat) != ncol(mat)) stop(name, ' is not square.')
  if (anyNA(mat)) stop(name, ' contains NA values.')
  if (!isTRUE(all.equal(mat, t(mat)))) stop(name, ' is not symmetric.')
  if (!is.null(ids)) {
    ids <- as.character(ids)
    if (!identical(rownames(mat), ids)) stop(name, ' row order does not match sample IDs.')
    if (!identical(colnames(mat), ids)) stop(name, ' column order does not match sample IDs.')
  }
  invisible(TRUE)
}

check_network_list <- function(network_list) {
  if (!is.list(network_list) || length(network_list) < 2L) {
    stop('network_list must be a list with at least two matrices.')
  }

  ref_ids <- NULL
  for (nm in names(network_list)) {
    mat <- as.matrix(network_list[[nm]])
    if (nrow(mat) != ncol(mat)) stop('Network ', nm, ' is not square.')
    if (anyNA(mat)) stop('Network ', nm, ' contains NA values.')
    if (!isTRUE(all.equal(mat, t(mat)))) stop('Network ', nm, ' is not symmetric.')
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
      stop('Network ', nm, ' must have row and column names.')
    }
    if (!identical(rownames(mat), colnames(mat))) {
      stop('Network ', nm, ' row and column names do not match.')
    }
    if (is.null(ref_ids)) {
      ref_ids <- rownames(mat)
    } else if (!identical(rownames(mat), ref_ids)) {
      stop('Network ', nm, ' does not match the sample order of the first network.')
    }
  }

  ref_ids
}

fit_scaler <- function(x) {
  x <- as.matrix(x)
  center <- colMeans(x, na.rm = TRUE)
  scale_vec <- apply(x, 2, stats::sd, na.rm = TRUE)
  scale_vec[!is.finite(scale_vec) | scale_vec == 0] <- 1
  list(center = center, scale = scale_vec)
}

apply_scaler <- function(x, scaler) {
  x <- as.matrix(x)
  x_centered <- sweep(x, 2, scaler$center, '-')
  x_scaled <- sweep(x_centered, 2, scaler$scale, '/')
  as.matrix(x_scaled)
}

build_affinity_from_data <- function(x, sample_ids, K, sigma) {
  x <- as.matrix(x)
  if (nrow(x) != length(sample_ids)) {
    stop('Row count does not match sample IDs.')
  }
  rownames(x) <- as.character(sample_ids)

  if (anyNA(x)) stop('Input matrix contains NA values.')
  if (ncol(x) == 0) stop('Input matrix has no columns.')

  if (ncol(x) == 1L) {
    if (stats::var(x[, 1], na.rm = TRUE) <= 0) {
      stop('Single-variable input has zero variance.')
    }
  } else {
    zero_var <- apply(x, 2, function(z) stats::var(z, na.rm = TRUE) <= 0)
    if (any(zero_var)) {
      stop('Input matrix contains zero-variance columns: ', paste(colnames(x)[zero_var], collapse = ', '))
    }
  }

  x_scaled <- scale(x)
  x_scaled <- as.matrix(x_scaled)
  rownames(x_scaled) <- sample_ids
  if (anyNA(x_scaled)) {
    stop('Scaling produced NA values.')
  }

  K_eff <- max(1L, min(as.integer(K), nrow(x_scaled) - 1L))
  dist_mat <- SNFtool::dist2(x_scaled, x_scaled)
  aff <- SNFtool::affinityMatrix(dist_mat, K = K_eff, sigma = sigma)
  rownames(aff) <- sample_ids
  colnames(aff) <- sample_ids
  check_affinity_matrix(aff, 'affinity_matrix', sample_ids)
  aff
}

run_snf_from_affinities <- function(network_list, K, T) {
  sample_ids <- check_network_list(network_list)
  K_eff <- max(1L, min(as.integer(K), length(sample_ids) - 1L))
  fused <- SNFtool::SNF(network_list, K = K_eff, t = as.integer(T))
  fused <- as.matrix(fused)
  rownames(fused) <- sample_ids
  colnames(fused) <- sample_ids
  check_affinity_matrix(fused, 'fused_network', sample_ids)
  fused
}

summarise_fused_network <- function(fused_network, K, T, network_names = NULL) {
  fused_network <- as.matrix(fused_network)
  check_affinity_matrix(fused_network, 'fused_network')

  offdiag <- fused_network[upper.tri(fused_network)]
  deg <- rowSums(fused_network)
  deg[deg <= 0] <- .Machine$double.eps
  d_inv_sqrt <- diag(1 / sqrt(deg))
  laplacian <- diag(nrow(fused_network)) - d_inv_sqrt %*% fused_network %*% d_inv_sqrt
  eigvals <- sort(eigen(laplacian, symmetric = TRUE, only.values = TRUE)$values)
  gap_limit <- min(10L, length(eigvals) - 1L)
  if (gap_limit < 1L) stop('Not enough eigenvalues to compute an eigengap.')
  gaps <- diff(eigvals[1:(gap_limit + 1L)])
  ordered_gaps <- order(gaps, decreasing = TRUE)
  best_k <- ordered_gaps[1L] + 1L
  second_best_k <- if (length(ordered_gaps) >= 2L) ordered_gaps[2L] + 1L else NA_integer_

  data.frame(
    metric = c(
      'n_samples',
      'n_networks',
      'K',
      'T',
      'min_affinity',
      'median_affinity',
      'mean_affinity',
      'max_affinity',
      'best_k_by_eigengap',
      'second_best_k_by_eigengap'
    ),
    value = c(
      nrow(fused_network),
      if (is.null(network_names)) NA_integer_ else length(network_names),
      as.integer(K),
      as.integer(T),
      min(offdiag, na.rm = TRUE),
      stats::median(offdiag, na.rm = TRUE),
      mean(offdiag, na.rm = TRUE),
      max(offdiag, na.rm = TRUE),
      best_k,
      second_best_k
    ),
    stringsAsFactors = FALSE
  )
}

run_eigengap_scan <- function(W, max_k = 10L) {
  W <- as.matrix(W)
  check_affinity_matrix(W, 'fused_network')

  deg <- rowSums(W)
  deg[deg <= 0] <- .Machine$double.eps
  d_inv_sqrt <- diag(1 / sqrt(deg))
  laplacian <- diag(nrow(W)) - d_inv_sqrt %*% W %*% d_inv_sqrt
  eigvals <- sort(eigen(laplacian, symmetric = TRUE, only.values = TRUE)$values)
  gap_limit <- min(as.integer(max_k), length(eigvals) - 1L)
  if (gap_limit < 1L) stop('Not enough eigenvalues to compute an eigengap.')

  gaps <- diff(eigvals[1:(gap_limit + 1L)])
  gap_table <- data.frame(
    k = 2:(gap_limit + 1L),
    eigval_left = eigvals[1:gap_limit],
    eigval_right = eigvals[2:(gap_limit + 1L)],
    gap = gaps,
    stringsAsFactors = FALSE
  )
  gap_table <- gap_table[order(-gap_table$gap, gap_table$k), , drop = FALSE]
  rownames(gap_table) <- NULL

  list(
    gap_table = gap_table,
    best_k = gap_table$k[1],
    second_best_k = if (nrow(gap_table) > 1L) gap_table$k[2] else NA_integer_,
    eigvals = eigvals[seq_len(min(length(eigvals), max_k + 1L))]
  )
}

run_spectral_clustering <- function(W, k) {
  W <- as.matrix(W)
  check_affinity_matrix(W, 'fused_network')
  k <- as.integer(k)
  if (k < 2L) stop('k must be at least 2.')
  if (k >= nrow(W)) stop('k must be smaller than the number of samples.')
  as.integer(SNFtool::spectralClustering(W, K = k))
}

save_cluster_outputs <- function(out_dir, prefix, sample_ids, fused_network, scan, chosen_k, cluster_labels, model_name) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  labels_df <- data.frame(
    sample_id = as.character(sample_ids),
    cluster = as.integer(cluster_labels),
    stringsAsFactors = FALSE
  )

  run_record <- data.frame(
    model_name = model_name,
    chosen_k = as.integer(chosen_k),
    second_best_k = if (!is.null(scan$second_best_k)) as.integer(scan$second_best_k) else NA_integer_,
    stringsAsFactors = FALSE
  )

  write.csv(labels_df, file.path(out_dir, 'best_k_cluster_labels.csv'), row.names = FALSE)
  write.csv(run_record, file.path(out_dir, 'cluster_run_record.csv'), row.names = FALSE)
  saveRDS(labels_df, file.path(out_dir, 'best_k_cluster_labels.rds'))
  saveRDS(run_record, file.path(out_dir, 'cluster_run_record.rds'))
  saveRDS(scan$gap_table, file.path(out_dir, 'eigengap_gap_table.rds'))
  write.csv(scan$gap_table, file.path(out_dir, 'eigengap_gap_table.csv'), row.names = FALSE)
  saveRDS(scan, file.path(out_dir, 'eigengap_scan.rds'))
  saveRDS(as.matrix(fused_network), file.path(out_dir, paste0(prefix, '_fused_network.rds')))
  saveRDS(list(
    model_name = model_name,
    sample_ids = as.character(sample_ids),
    fused_network = as.matrix(fused_network),
    eigengap_candidates = scan$gap_table,
    best_k = scan$best_k,
    second_best_k = scan$second_best_k,
    best_k_cluster_labels = labels_df
  ), file.path(out_dir, paste0(prefix, '_cluster_results.rds')))

  invisible(list(labels = labels_df, run_record = run_record))
}

compare_networks_upper_triangle <- function(W_main, W_branch) {
  ids <- intersect(rownames(W_main), rownames(W_branch))
  if (length(ids) < 2L) stop('Not enough shared samples to compare fused networks.')

  W_main <- as.matrix(W_main[ids, ids, drop = FALSE])
  W_branch <- as.matrix(W_branch[ids, ids, drop = FALSE])

  v1 <- W_main[upper.tri(W_main)]
  v2 <- W_branch[upper.tri(W_branch)]

  data.frame(
    metric = c('pearson_cor', 'spearman_cor', 'mean_abs_diff', 'median_abs_diff'),
    value = c(
      suppressWarnings(stats::cor(v1, v2, method = 'pearson')),
      suppressWarnings(stats::cor(v1, v2, method = 'spearman')),
      mean(abs(v1 - v2), na.rm = TRUE),
      stats::median(abs(v1 - v2), na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

pairwise_co_clustering_agreement <- function(labels_a, labels_b) {
  labels_a <- as.integer(labels_a)
  labels_b <- as.integer(labels_b)
  if (length(labels_a) != length(labels_b)) stop('Label vectors must have the same length.')
  n <- length(labels_a)
  if (n < 2L) stop('Need at least two samples for pairwise agreement.')

  same_a <- outer(labels_a, labels_a, '==')
  same_b <- outer(labels_b, labels_b, '==')
  upper_idx <- upper.tri(same_a)
  mean(same_a[upper_idx] == same_b[upper_idx], na.rm = TRUE)
}

save_representation_comparison_outputs <- function(out_dir, prefix, merged_labels, summary_obj) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(merged_labels, file.path(out_dir, paste0(prefix, '_merged_cluster_labels.csv')), row.names = FALSE)
  write.csv(as.data.frame.matrix(summary_obj$contingency), file.path(out_dir, paste0(prefix, '_contingency_table.csv')), row.names = TRUE)

  summary_table <- data.frame(
    metric = c('pairwise_co_clustering_agreement'),
    value = c(summary_obj$pairwise_agreement),
    stringsAsFactors = FALSE
  )

  write.csv(summary_table, file.path(out_dir, paste0(prefix, '_comparison_summary.csv')), row.names = FALSE)
  saveRDS(summary_table, file.path(out_dir, paste0(prefix, '_comparison_summary.rds')))
  saveRDS(summary_obj$contingency, file.path(out_dir, paste0(prefix, '_contingency_table.rds')))

  png(file.path(out_dir, paste0(prefix, '_contingency_heatmap.png')), width = 2000, height = 1500, res = 300)
  contingency_mat <- as.matrix(summary_obj$contingency)
  image(
    x = seq_len(ncol(contingency_mat)),
    y = seq_len(nrow(contingency_mat)),
    z = t(contingency_mat[nrow(contingency_mat):1, , drop = FALSE]),
    axes = FALSE,
    xlab = 'Comparison cluster 2',
    ylab = 'Comparison cluster 1',
    main = 'Cluster comparison contingency'
  )
  axis(1, at = seq_len(ncol(contingency_mat)), labels = colnames(contingency_mat))
  axis(2, at = seq_len(nrow(contingency_mat)), labels = rev(rownames(contingency_mat)))
  box()
  dev.off()

  invisible(TRUE)
}

# Top-20 metabolites from the provided table.
top20_metabolites <- c(
  'palmitoleamide (16:1)*',
  'myristoleamide (14:1)*',
  '3-methoxytyrosine',
  '3-methoxytyramine sulfate',
  'myristamide (14:0)*',
  'cyclo(leu-pro)',
  'X-21733',
  'm-tyramine sulfate',
  'margaramide (17:0)*',
  'perfluorooctanesulfonate (PFOS)',
  'heptadecenamide (17:1)*',
  'linolenamide (18:3)*',
  'linoleamide (18:2n6)',
  'fibrinopeptide B (1-13)**',
  'dopamine 3-O-sulfate',
  'X-12410',
  'threonate',
  'biliverdin',
  'p-cresol glucuronide*',
  'acetoacetate'
)

resolve_metabolite_columns <- function(df, requested_names = top20_metabolites) {
  available <- names(df)
  available_norm <- normalize_name(available)
  names(available_norm) <- available

  resolved <- character(0)
  missing <- character(0)

  for (target in requested_names) {
    if (target %in% available) {
      resolved <- c(resolved, target)
      next
    }
    target_norm <- normalize_name(target)
    hit <- names(available_norm)[available_norm == target_norm]
    if (length(hit) >= 1L) {
      resolved <- c(resolved, hit[1])
    } else {
      missing <- c(missing, target)
    }
  }

  list(resolved = resolved, missing = missing)
}

make_top20_manifest <- function(network_names, network_types, source_columns, output_files) {
  data.frame(
    network_name = network_names,
    network_type = network_types,
    source_column = source_columns,
    output_file = output_files,
    stringsAsFactors = FALSE
  )
}

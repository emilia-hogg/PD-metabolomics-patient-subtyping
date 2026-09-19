# Shared helper functions for comparing clustering solutions in the
# representation-comparison sensitivity analysis.
#
# These functions load and validate cluster results, merge cluster labels
# between the two- and four-network solutions, calculate pairwise
# co-clustering agreement, and save comparison tables and a heatmap.
#
# No project-specific file paths or parameters need to be set in this file.

load_cluster_results <- function(path, expected_name = NULL) {
  if (!file.exists(path)) {
    stop("Missing cluster results file: ", path)
  }

  obj <- readRDS(path)

  required <- c("sample_ids", "cluster_labels", "eigengap_scan", "best_k", "second_best_k")
  missing <- setdiff(required, names(obj))
  if (length(missing) > 0) {
    stop("Cluster results file is missing: ", paste(missing, collapse = ", "))
  }

  if (!is.null(expected_name) && !is.null(obj$model_name) && obj$model_name != expected_name) {
    warning("Model name mismatch for ", path, ": expected ", expected_name, ", found ", obj$model_name)
  }

  obj
}

merge_cluster_labels <- function(res_2network, res_4network) {
  df2 <- as.data.frame(res_2network$cluster_labels, stringsAsFactors = FALSE)
  df4 <- as.data.frame(res_4network$cluster_labels, stringsAsFactors = FALSE)

  names(df2) <- c("sample_id", "cluster_2network")
  names(df4) <- c("sample_id", "cluster_4network")

  merged <- merge(df2, df4, by = "sample_id", all = FALSE, sort = FALSE)

  if (nrow(merged) == 0) {
    stop("No overlapping sample IDs between the two cluster solutions.")
  }

  merged
}

pairwise_co_clustering_agreement <- function(labels_a, labels_b) {
  labels_a <- as.integer(labels_a)
  labels_b <- as.integer(labels_b)

  if (length(labels_a) != length(labels_b)) {
    stop("Label vectors must have the same length.")
  }

  n <- length(labels_a)
  if (n < 2L) {
    stop("Need at least two samples for pairwise agreement.")
  }

  same_a <- outer(labels_a, labels_a, "==")
  same_b <- outer(labels_b, labels_b, "==")

  upper_idx <- upper.tri(same_a)
  mean(same_a[upper_idx] == same_b[upper_idx], na.rm = TRUE)
}

summarise_representation_comparison <- function(merged_labels) {
  if (!all(c("cluster_2network", "cluster_4network") %in% names(merged_labels))) {
    stop("merged_labels must contain cluster_2network and cluster_4network.")
  }

  pairwise_agreement <- pairwise_co_clustering_agreement(
    merged_labels$cluster_2network,
    merged_labels$cluster_4network
  )

  contingency <- table(merged_labels$cluster_2network, merged_labels$cluster_4network)

  list(
    pairwise_agreement = pairwise_agreement,
    contingency = contingency
  )
}

save_representation_comparison_outputs <- function(out_dir, prefix, merged_labels, summary_obj) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    merged_labels,
    file.path(out_dir, paste0(prefix, "_merged_cluster_labels.csv")),
    row.names = FALSE
  )

  write.csv(
    as.data.frame.matrix(summary_obj$contingency),
    file.path(out_dir, paste0(prefix, "_contingency_table.csv")),
    row.names = TRUE
  )

  summary_table <- data.frame(
    metric = c("pairwise_co_clustering_agreement"),
    value = c(summary_obj$pairwise_agreement),
    stringsAsFactors = FALSE
  )

  write.csv(
    summary_table,
    file.path(out_dir, paste0(prefix, "_comparison_summary.csv")),
    row.names = FALSE
  )

  saveRDS(summary_table, file.path(out_dir, paste0(prefix, "_comparison_summary.rds")))
  saveRDS(summary_obj$contingency, file.path(out_dir, paste0(prefix, "_contingency_table.rds")))

  png(file.path(out_dir, paste0(prefix, "_contingency_heatmap.png")), width = 2000, height = 1500, res = 300)
  contingency_mat <- as.matrix(summary_obj$contingency)
  image(
    x = seq_len(ncol(contingency_mat)),
    y = seq_len(nrow(contingency_mat)),
    z = t(contingency_mat[nrow(contingency_mat):1, , drop = FALSE]),
    axes = FALSE,
    xlab = "4-network cluster",
    ylab = "2-network cluster",
    main = "Representation comparison contingency"
  )
  axis(1, at = seq_len(ncol(contingency_mat)), labels = colnames(contingency_mat))
  axis(2, at = seq_len(nrow(contingency_mat)), labels = rev(rownames(contingency_mat)))
  box()
  dev.off()

  invisible(TRUE)
}

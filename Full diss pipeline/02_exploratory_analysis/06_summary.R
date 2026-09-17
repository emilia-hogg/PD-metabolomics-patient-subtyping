# 02_exploratory_analysis/06_summary.R

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

summary_dir <- file.path(paths$exploratory, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

safe_read_csv <- function(path) {
  if (file.exists(path)) {
    read.csv(path, stringsAsFactors = FALSE)
  } else {
    NULL
  }
}

safe_read_rds <- function(path) {
  if (file.exists(path)) {
    readRDS(path)
  } else {
    NULL
  }
}

# -------------------------------------------------------------------
# Locate the main exploratory outputs
# -------------------------------------------------------------------
pca_var_file <- file.path(paths$exploratory, "pca", "pca_variance_summary.csv")
pca_scores_file <- file.path(paths$exploratory, "pca", "pca_scores.rds")

kmeans_pca_file <- file.path(paths$exploratory, "kmeans", "pca", "pca_kmeans_all_metrics.csv")
kmeans_umap_file <- file.path(paths$exploratory, "kmeans", "umap", "UMAP_metrics.csv")

dbscan_summary_file <- file.path(paths$exploratory, "dbscan", "dbscan_parameter_summary.csv")
hdbscan_summary_file <- file.path(paths$exploratory, "hdbscan", "hdbscan_all_minPts_summary.csv")

pca_var <- safe_read_csv(pca_var_file)
pca_scores <- safe_read_rds(pca_scores_file)
kmeans_pca <- safe_read_csv(kmeans_pca_file)
kmeans_umap <- safe_read_csv(kmeans_umap_file)
dbscan_sum <- safe_read_csv(dbscan_summary_file)
hdbscan_sum <- safe_read_csv(hdbscan_summary_file)

# -------------------------------------------------------------------
# Save a simple availability overview
# -------------------------------------------------------------------
overview <- data.frame(
  analysis = c("PCA", "K-means on PCA", "K-means on UMAP", "DBSCAN", "HDBSCAN"),
  output_found = c(
    !is.null(pca_var),
    !is.null(kmeans_pca),
    !is.null(kmeans_umap),
    !is.null(dbscan_sum),
    !is.null(hdbscan_sum)
  ),
  stringsAsFactors = FALSE
)

write.csv(overview, file.path(summary_dir, "exploratory_output_overview.csv"), row.names = FALSE)
saveRDS(overview, file.path(summary_dir, "exploratory_output_overview.rds"))

# -------------------------------------------------------------------
# PCA summary
# -------------------------------------------------------------------
if (!is.null(pca_var)) {
  write.csv(pca_var, file.path(summary_dir, "pca_variance_summary.csv"), row.names = FALSE)

  pca_top10 <- head(pca_var, 10)
  write.csv(pca_top10, file.path(summary_dir, "pca_variance_top10.csv"), row.names = FALSE)

  best_10pc <- pca_var[which.min(abs(as.numeric(sub("PC", "", pca_var$PC)) - 10)), , drop = FALSE]
  saveRDS(best_10pc, file.path(summary_dir, "pca_variance_pc10.rds"))
}

# -------------------------------------------------------------------
# K-means summaries
# -------------------------------------------------------------------
if (!is.null(kmeans_pca)) {
  write.csv(kmeans_pca, file.path(summary_dir, "kmeans_pca_all_metrics.csv"), row.names = FALSE)

  # Best k by average silhouette within each PC set
  if (all(c("pc_set", "k", "avg_silhouette", "tot_withinss") %in% names(kmeans_pca))) {
    best_by_pc <- do.call(
      rbind,
      lapply(split(kmeans_pca, kmeans_pca$pc_set), function(df) {
        df <- df[order(-df$avg_silhouette, df$tot_withinss), , drop = FALSE]
        df[1, , drop = FALSE]
      })
    )
    rownames(best_by_pc) <- NULL
    write.csv(best_by_pc, file.path(summary_dir, "kmeans_pca_best_by_pcset.csv"), row.names = FALSE)
    saveRDS(best_by_pc, file.path(summary_dir, "kmeans_pca_best_by_pcset.rds"))
  }
}

if (!is.null(kmeans_umap)) {
  write.csv(kmeans_umap, file.path(summary_dir, "kmeans_umap_metrics.csv"), row.names = FALSE)

  if (all(c("k", "avg_silhouette", "tot_withinss") %in% names(kmeans_umap))) {
    best_umap <- kmeans_umap[order(-kmeans_umap$avg_silhouette, kmeans_umap$tot_withinss), , drop = FALSE][1, , drop = FALSE]
    write.csv(best_umap, file.path(summary_dir, "kmeans_umap_best_k.csv"), row.names = FALSE)
    saveRDS(best_umap, file.path(summary_dir, "kmeans_umap_best_k.rds"))
  }
}

# -------------------------------------------------------------------
# DBSCAN summary
# -------------------------------------------------------------------
if (!is.null(dbscan_sum)) {
  write.csv(dbscan_sum, file.path(summary_dir, "dbscan_parameter_summary.csv"), row.names = FALSE)

  if (all(c("eps", "minPts", "n_clusters", "n_noise") %in% names(dbscan_sum))) {
    # Simple heuristic summary: favour more clusters and fewer noise points
    dbscan_ranked <- dbscan_sum[order(-dbscan_sum$n_clusters, dbscan_sum$n_noise, dbscan_sum$eps), , drop = FALSE]
    best_dbscan <- dbscan_ranked[1, , drop = FALSE]
    write.csv(best_dbscan, file.path(summary_dir, "dbscan_best_setting.csv"), row.names = FALSE)
    saveRDS(best_dbscan, file.path(summary_dir, "dbscan_best_setting.rds"))
  }
}

# -------------------------------------------------------------------
# HDBSCAN summary
# -------------------------------------------------------------------
if (!is.null(hdbscan_sum)) {
  write.csv(hdbscan_sum, file.path(summary_dir, "hdbscan_all_minPts_summary.csv"), row.names = FALSE)

  # Keep only the clean metrics rows if they exist
  metric_cols <- c("minPts", "n_clusters", "n_noise", "mean_membership_prob", "mean_outlier_score")
  if (all(metric_cols %in% names(hdbscan_sum))) {
    hdbscan_metrics <- unique(hdbscan_sum[, metric_cols, drop = FALSE])
    write.csv(hdbscan_metrics, file.path(summary_dir, "hdbscan_metrics_clean.csv"), row.names = FALSE)

    # Simple ranking: more clusters, less noise, higher membership probability
    hdbscan_ranked <- hdbscan_metrics[order(-hdbscan_metrics$n_clusters, hdbscan_metrics$n_noise, -hdbscan_metrics$mean_membership_prob), , drop = FALSE]
    best_hdbscan <- hdbscan_ranked[1, , drop = FALSE]
    write.csv(best_hdbscan, file.path(summary_dir, "hdbscan_best_setting.csv"), row.names = FALSE)
    saveRDS(best_hdbscan, file.path(summary_dir, "hdbscan_best_setting.rds"))
  }
}

# -------------------------------------------------------------------
# Combined compact report table
# -------------------------------------------------------------------
report_rows <- list()

if (!is.null(pca_var)) {
  report_rows[["PCA"]] <- data.frame(
    analysis = "PCA",
    summary = paste0(
      "PC1 explained ", round(100 * pca_var$variance_explained[1], 1), "% of variance; ",
      "PC10 cumulative variance = ",
      round(100 * pca_var$cumulative_variance[min(10, nrow(pca_var))], 1), "%"
    ),
    stringsAsFactors = FALSE
  )
}

if (!is.null(kmeans_pca) && all(c("pc_set", "k", "avg_silhouette") %in% names(kmeans_pca))) {
  best_by_pc <- do.call(
    rbind,
    lapply(split(kmeans_pca, kmeans_pca$pc_set), function(df) {
      df <- df[order(-df$avg_silhouette, df$tot_withinss), , drop = FALSE]
      df[1, , drop = FALSE]
    })
  )
  report_rows[["KMEANS_PCA"]] <- data.frame(
    analysis = "K-means PCA",
    summary = paste0(
      "Best silhouette across PC sets ranged from ",
      paste(best_by_pc$k, collapse = ", ")
    ),
    stringsAsFactors = FALSE
  )
}

if (!is.null(kmeans_umap) && all(c("k", "avg_silhouette") %in% names(kmeans_umap))) {
  best_umap <- kmeans_umap[order(-kmeans_umap$avg_silhouette, kmeans_umap$tot_withinss), , drop = FALSE][1, , drop = FALSE]
  report_rows[["KMEANS_UMAP"]] <- data.frame(
    analysis = "K-means UMAP",
    summary = paste0(
      "Best UMAP k = ", best_umap$k,
      " with average silhouette = ", round(best_umap$avg_silhouette, 3)
    ),
    stringsAsFactors = FALSE
  )
}

if (!is.null(dbscan_sum) && all(c("eps", "n_clusters", "n_noise") %in% names(dbscan_sum))) {
  best_dbscan <- dbscan_sum[order(-dbscan_sum$n_clusters, dbscan_sum$n_noise, dbscan_sum$eps), , drop = FALSE][1, , drop = FALSE]
  report_rows[["DBSCAN"]] <- data.frame(
    analysis = "DBSCAN",
    summary = paste0(
      "Best eps = ", best_dbscan$eps,
      " with ", best_dbscan$n_clusters, " clusters and ", best_dbscan$n_noise, " noise points"
    ),
    stringsAsFactors = FALSE
  )
}

if (!is.null(hdbscan_sum) && all(c("minPts", "n_clusters", "n_noise", "mean_membership_prob") %in% names(hdbscan_sum))) {
  hdbscan_metrics <- unique(hdbscan_sum[, c("minPts", "n_clusters", "n_noise", "mean_membership_prob", "mean_outlier_score"), drop = FALSE])
  best_hdbscan <- hdbscan_metrics[order(-hdbscan_metrics$n_clusters, hdbscan_metrics$n_noise, -hdbscan_metrics$mean_membership_prob), , drop = FALSE][1, , drop = FALSE]
  report_rows[["HDBSCAN"]] <- data.frame(
    analysis = "HDBSCAN",
    summary = paste0(
      "Best minPts = ", best_hdbscan$minPts,
      " with ", best_hdbscan$n_clusters, " clusters and mean membership probability ",
      round(best_hdbscan$mean_membership_prob, 3)
    ),
    stringsAsFactors = FALSE
  )
}

if (length(report_rows) > 0) {
  report_df <- do.call(rbind, report_rows)
  rownames(report_df) <- NULL
  write.csv(report_df, file.path(summary_dir, "exploratory_summary_report.csv"), row.names = FALSE)
  saveRDS(report_df, file.path(summary_dir, "exploratory_summary_report.rds"))
}

# -------------------------------------------------------------------
# Final text note
# -------------------------------------------------------------------
note_file <- file.path(summary_dir, "exploratory_summary_note.txt")
note_lines <- c(
  "Exploratory analysis summary generated successfully.",
  paste0("PCA file found: ", !is.null(pca_var)),
  paste0("K-means PCA file found: ", !is.null(kmeans_pca)),
  paste0("K-means UMAP file found: ", !is.null(kmeans_umap)),
  paste0("DBSCAN file found: ", !is.null(dbscan_sum)),
  paste0("HDBSCAN file found: ", !is.null(hdbscan_sum))
)
writeLines(note_lines, note_file)

message("Exploratory summary complete.")
message("Saved outputs to: ", summary_dir)
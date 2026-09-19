# 02_exploratory_analysis/05_hdbscan.R
# This script performs HDBSCAN clustering on the first 20 principal component
# scores from the metabolomics data. HDBSCAN is run across several `minPts`
# values, and the resulting cluster assignments, membership probabilities,
# outlier scores, cluster summaries, and comparison metrics are saved.
#
# The resulting HDBSCAN clusters are also visualised using the 2D UMAP
# coordinates, with one UMAP plot produced for each `minPts` value.
#
# Fixed settings:
# - PCA input: PC1-PC20
# - minPts values: 3, 5, 10, 20, 40
#
# Before running:
# - Set `project_root` to the local project directory.
# - `pca_scores.rds` must exist at `paths$exploratory/pca/pca_scores.rds`.
# - `umap_coords_clean.csv` must exist at
#   `paths$exploratory/umap/umap_coords_clean.csv`.
#
# Note:
# - HDBSCAN clustering itself is performed on the first 20 PCs.
# - UMAP is used only to visualise the resulting HDBSCAN cluster assignments.
#
# Outputs are saved under `paths$exploratory/hdbscan/`, including cluster
# assignments for each `minPts`, summary and metric tables, UMAP cluster plots,
# and combined RDS summaries of all HDBSCAN results.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

pca_file <- file.path(paths$exploratory, "pca", "pca_scores.rds")
umap_file <- file.path(paths$exploratory, "umap", "umap_coords_clean.csv")

if (!file.exists(pca_file)) stop("Missing PCA scores file: ", pca_file)
if (!file.exists(umap_file)) stop("Missing UMAP coordinates file: ", umap_file)

pca_scores <- readRDS(pca_file)
umap_coords <- read.csv(umap_file, stringsAsFactors = FALSE)

if (!"sample_id" %in% names(pca_scores)) stop("pca_scores must contain sample_id.")
if (!"sample_id" %in% names(umap_coords)) stop("UMAP file must contain sample_id.")
if (!all(c("UMAP1", "UMAP2") %in% names(umap_coords))) {
  stop("UMAP file must contain UMAP1 and UMAP2.")
}

pc_cols <- paste0("PC", 1:20)
missing_pc_cols <- setdiff(pc_cols, names(pca_scores))
if (length(missing_pc_cols) > 0) {
  stop("PCA scores missing required columns: ", paste(missing_pc_cols, collapse = ", "))
}

X <- as.matrix(pca_scores[, pc_cols, drop = FALSE])
rownames(X) <- pca_scores$sample_id

if (anyNA(X)) stop("PC score matrix contains NA values.")
if (anyDuplicated(rownames(X))) stop("Duplicated sample IDs in PCA matrix.")

out_dir <- file.path(paths$exploratory, "hdbscan")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Save input matrix for reproducibility
saveRDS(X, file.path(out_dir, "hdbscan_input_pc1_20.rds"))

minPts_values <- c(3, 5, 10, 20, 40)
all_results <- list()
summary_list <- list()

for (m in minPts_values) {
  message("Running HDBSCAN with minPts = ", m)

  fit <- dbscan::hdbscan(X, minPts = m)

  if (is.null(fit$cluster) || length(fit$cluster) != nrow(X)) {
    stop("HDBSCAN cluster vector is invalid for minPts = ", m)
  }

  membership_prob <- if (!is.null(fit$membership_prob) && length(fit$membership_prob) == nrow(X)) {
    fit$membership_prob
  } else {
    rep(NA_real_, nrow(X))
  }

  outlier_score <- if (!is.null(fit$outlier_scores) && length(fit$outlier_scores) == nrow(X)) {
    fit$outlier_scores
  } else {
    rep(NA_real_, nrow(X))
  }

  cluster_df <- data.frame(
    sample_id = rownames(X),
    cluster = factor(fit$cluster),
    membership_prob = membership_prob,
    outlier_score = outlier_score,
    stringsAsFactors = FALSE
  )

  write.csv(
    cluster_df,
    file.path(out_dir, paste0("hdbscan_minPts_", m, "_clusters.csv")),
    row.names = FALSE
  )

  tab <- table(cluster_df$cluster)
  summary_df <- data.frame(
    minPts = m,
    cluster = names(tab),
    n = as.integer(tab),
    stringsAsFactors = FALSE
  )

  # Extra compact metrics for quick comparison
  n_clusters <- length(setdiff(unique(fit$cluster), 0))
  n_noise <- sum(fit$cluster == 0)
  mean_membership <- mean(membership_prob, na.rm = TRUE)
  mean_outlier <- mean(outlier_score, na.rm = TRUE)

  metrics_df <- data.frame(
    minPts = m,
    n_clusters = n_clusters,
    n_noise = n_noise,
    mean_membership_prob = mean_membership,
    mean_outlier_score = mean_outlier,
    stringsAsFactors = FALSE
  )

  summary_list[[as.character(m)]] <- cbind(metrics_df, summary_df)
  all_results[[as.character(m)]] <- list(
    fit = fit,
    cluster_df = cluster_df,
    summary_df = summary_df,
    metrics_df = metrics_df
  )

  write.csv(
    summary_df,
    file.path(out_dir, paste0("hdbscan_minPts_", m, "_summary.csv")),
    row.names = FALSE
  )

  write.csv(
    metrics_df,
    file.path(out_dir, paste0("hdbscan_minPts_", m, "_metrics.csv")),
    row.names = FALSE
  )

  plot_df <- merge(umap_coords, cluster_df, by = "sample_id", all.x = TRUE)

  png(
    file.path(out_dir, paste0("hdbscan_umap_minPts_", m, ".png")),
    width = 2000, height = 1500, res = 300
  )
  plot(
    plot_df$UMAP1,
    plot_df$UMAP2,
    col = as.integer(plot_df$cluster),
    pch = 16,
    cex = 0.7,
    xlab = "UMAP1",
    ylab = "UMAP2",
    main = paste0("HDBSCAN on first 20 PCs (minPts = ", m, ")")
  )
  dev.off()
}

all_summary <- do.call(rbind, summary_list)
write.csv(all_summary, file.path(out_dir, "hdbscan_all_minPts_summary.csv"), row.names = FALSE)
saveRDS(all_results, file.path(out_dir, "hdbscan_all_results.rds"))
saveRDS(all_summary, file.path(out_dir, "hdbscan_all_minPts_summary.rds"))

# A quick overview plot
png(file.path(out_dir, "hdbscan_minPts_overview.png"), width = 2200, height = 1500, res = 300)
par(mfrow = c(2, 1), mar = c(5, 5, 3, 1))
plot(
  unique(all_summary$minPts),
  tapply(all_summary$n_clusters, all_summary$minPts, function(x) x[1]),
  type = "b",
  pch = 16,
  xlab = "minPts",
  ylab = "Number of clusters",
  main = "HDBSCAN clusters by minPts"
)
plot(
  unique(all_summary$minPts),
  tapply(all_summary$n, all_summary$minPts, sum),
  type = "b",
  pch = 16,
  xlab = "minPts",
  ylab = "Total points assigned across clusters",
  main = "Clustered points by minPts"
)
par(mfrow = c(1, 1))
dev.off()

message("HDBSCAN complete.")
message("Saved outputs to: ", out_dir)

# 04_dbscan.R  (CORRECTED)
# This script performs DBSCAN clustering on the first 20 principal component
# scores from the metabolomics data. It first examines k-nearest-neighbour
# distances, identifies an elbow-based epsilon value, constructs an epsilon
# grid from the observed distances, and scans DBSCAN solutions across that grid.
#
# Fixed settings:
# - PCA input: PC1-PC20
# - minPts = 40
# - 40 epsilon values are evaluated across the data-derived epsilon range
# - The elbow epsilon is identified numerically as the point furthest from
#   the straight line joining the first and last sorted kNN distances.
# - kNN diagnostic plots are produced for k = 20, 40, and 60.
#
# Before running:
# - Set `project_root` to the local project directory.
# - `pca_scores.rds` must exist at `paths$exploratory/pca/pca_scores.rds`.
# - `umap_coords_clean.csv` must exist at
#   `paths$exploratory/umap/umap_coords_clean.csv`.
#   Note: the current script checks for and loads this UMAP file, but does
#   not otherwise use the UMAP coordinates in the DBSCAN analysis (redundancy
#   related to note in '03_kmeans.R').
#
# Outputs are saved under `paths$exploratory/dbscan/`, including the PC1-PC20
# DBSCAN input matrix, kNN distance plots, parameter-scan summary and plot,
# and the DBSCAN results across all epsilon values.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

pca_file  <- file.path(paths$exploratory, "pca", "pca_scores.rds")
umap_file <- file.path(paths$exploratory, "umap", "umap_coords_clean.csv")
if (!file.exists(pca_file))  stop("Missing PCA scores file: ", pca_file)
if (!file.exists(umap_file)) stop("Missing UMAP coordinates file: ", umap_file)

pca_scores  <- readRDS(pca_file)
umap_coords <- read.csv(umap_file, stringsAsFactors = FALSE)

pc_cols <- paste0("PC", 1:20)
if (length(setdiff(pc_cols, names(pca_scores)))) {
  stop("PCA scores missing required columns.")
}

X <- as.matrix(pca_scores[, pc_cols, drop = FALSE])
rownames(X) <- pca_scores$sample_id
if (anyNA(X)) stop("PC score matrix contains NA values.")
if (anyDuplicated(rownames(X))) stop("Duplicated sample IDs in PCA matrix.")

out_dir <- file.path(paths$exploratory, "dbscan")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(X, file.path(out_dir, "dbscan_input_pc1_20.rds"))

minPts <- 40

# ------------------------------------------------------------
# 1. Observed kNN distances, and the elbow located numerically
# ------------------------------------------------------------
knn_d <- sort(dbscan::kNNdist(X, k = minPts))
n     <- length(knn_d)

# Kneedle-style elbow: the point furthest from the straight line
# joining the first and last points of the sorted curve.
x  <- seq_len(n); y <- knn_d
x1 <- x[1]; y1 <- y[1]; x2 <- x[n]; y2 <- y[n]
perp <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) /
        sqrt((y2 - y1)^2 + (x2 - x1)^2)
elbow_idx <- which.max(perp)
elbow_eps <- knn_d[elbow_idx]

cat("\n--- kNN distance summary (k = minPts = ", minPts, ") ---\n", sep = "")
print(summary(knn_d))
cat("elbow at index ", elbow_idx, " of ", n,
    "  ->  epsilon = ", round(elbow_eps, 2), "\n", sep = "")

# ------------------------------------------------------------
# 2. Epsilon grid derived from the data, not hardcoded
# ------------------------------------------------------------
eps_lo   <- floor(min(knn_d))
eps_hi   <- ceiling(quantile(knn_d, 0.999, names = FALSE))
eps_grid <- seq(eps_lo, eps_hi, length.out = 40)

cat("epsilon grid: ", round(eps_lo, 1), " to ", round(eps_hi, 1),
    " in ", length(eps_grid), " steps\n\n", sep = "")

if (elbow_eps < eps_lo || elbow_eps > eps_hi) {
  warning("Elbow falls outside the epsilon grid; widen the grid.")
}

# ------------------------------------------------------------
# 3. kNN distance plot with the elbow and grid marked
# ------------------------------------------------------------
png(file.path(out_dir, "dbscan_knn_distance_plots.png"),
    width = 2200, height = 1800, res = 300)
par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))
for (m in c(20, 40, 60)) {
  dbscan::kNNdistplot(X, k = m)
  title(main = paste0("kNN distance plot (k = ", m, ")"))
  if (m == minPts) {
    abline(h = range(eps_grid), col = "grey50", lty = 3)
    legend("topleft", bty = "n", cex = 0.8,
           legend = c(sprintf("grid: %.0f to %.0f", eps_lo, eps_hi)),
           lty = c(2, 3), col = c("red", "grey50"))
  }
}
par(mfrow = c(1, 1))
dev.off()

# ------------------------------------------------------------
# 4. Scan
# ------------------------------------------------------------
summary_rows <- list(); results <- list()

for (eps in eps_grid) {
  fit <- dbscan::dbscan(X, eps = eps, minPts = minPts)
  cl  <- fit$cluster
  tab <- table(cl[cl != 0])

  summary_rows[[as.character(eps)]] <- data.frame(
    eps                = eps,
    minPts             = minPts,
    n_clusters         = length(tab),
    n_noise            = sum(cl == 0),
    prop_noise         = mean(cl == 0),
    largest_cluster    = if (length(tab)) max(tab) else 0L,
    prop_in_largest    = if (length(tab)) max(tab) / length(cl) else 0,
    second_largest     = if (length(tab) > 1) sort(tab, decreasing = TRUE)[2] else 0L,
    stringsAsFactors = FALSE
  )
  results[[as.character(eps)]] <- data.frame(
    sample_id = rownames(X), cluster = cl, stringsAsFactors = FALSE)
}

summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[order(summary_df$eps), ]
write.csv(summary_df, file.path(out_dir, "dbscan_parameter_summary.csv"),
          row.names = FALSE)
saveRDS(results, file.path(out_dir, "dbscan_all_results.rds"))

cat("--- scan ---\n")
cat(sprintf("%8s %10s %9s %16s %16s\n",
            "eps", "clusters", "noise %", "largest cluster", "2nd largest"))
for (i in seq_len(nrow(summary_df))) {
  cat(sprintf("%8.1f %10d %8.1f%% %15s%% %16d\n",
              summary_df$eps[i], summary_df$n_clusters[i],
              100 * summary_df$prop_noise[i],
              formatC(100 * summary_df$prop_in_largest[i], format = "f", digits = 1),
              summary_df$second_largest[i]))
}

# ------------------------------------------------------------
# 5. The plot that carries the result
# ------------------------------------------------------------
png(file.path(out_dir, "dbscan_parameter_summary_plot.png"),
    width = 2000, height = 1600, res = 300)
par(mar = c(4.5, 4.5, 3, 4.5))
plot(summary_df$eps, 100 * summary_df$prop_noise, type = "l", lwd = 2,
     col = "grey30", ylim = c(0, 100),
     xlab = expression(epsilon), ylab = "% of patients",
     main = "DBSCAN across the observed distance range")
lines(summary_df$eps, 100 * summary_df$prop_in_largest, lwd = 2, col = "#1f78b4")
abline(v = elbow_eps, col = "red", lty = 2)
legend("right", bty = "n", cex = 0.85,
       legend = c("classified as noise", "in the largest cluster", "elbow"),
       lty = c(1, 1, 2), lwd = c(2, 2, 1), col = c("grey30", "#1f78b4", "red"))
par(new = TRUE)
plot(summary_df$eps, summary_df$n_clusters, type = "l", lty = 3,
     axes = FALSE, xlab = "", ylab = "")
axis(4); mtext("number of clusters", side = 4, line = 3, cex = 0.9)
dev.off()

message("DBSCAN complete. Outputs: ", out_dir)

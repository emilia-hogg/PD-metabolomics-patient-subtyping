# 02_exploratory_analysis/01_pca.R

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

if (!file.exists(files$metabolomics_processed)) {
  stop("Missing metabolomics processed file: ", files$metabolomics_processed)
}

metab_obj <- readRDS(files$metabolomics_processed)

if (!all(c("metab_scaled", "sample_ids") %in% names(metab_obj))) {
  stop("metabolomics_processed.rds must contain 'metab_scaled' and 'sample_ids'.")
}

metab_scaled <- as.matrix(metab_obj$metab_scaled)
sample_ids <- as.character(metab_obj$sample_ids)

if (anyDuplicated(sample_ids)) {
  stop("Duplicated sample IDs found in metabolomics sample_ids.")
}
if (!identical(rownames(metab_scaled), sample_ids)) {
  stop("Row order of metab_scaled does not match sample_ids.")
}
if (anyNA(metab_scaled)) {
  stop("metab_scaled contains NA values.")
}

pca_dir <- file.path(paths$exploratory, "pca")
dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)

pca <- prcomp(metab_scaled, center = FALSE, scale. = FALSE)

pca_scores <- data.frame(
  sample_id = sample_ids,
  pca$x,
  check.names = FALSE,
  row.names = NULL
)

variance_explained <- (pca$sdev^2) / sum(pca$sdev^2)
variance_summary <- data.frame(
  PC = paste0("PC", seq_along(variance_explained)),
  variance_explained = variance_explained,
  cumulative_variance = cumsum(variance_explained),
  stringsAsFactors = FALSE
)

saveRDS(pca, file.path(pca_dir, "pca_model.rds"))
saveRDS(pca_scores, file.path(pca_dir, "pca_scores.rds"))
saveRDS(variance_summary, file.path(pca_dir, "pca_variance_summary.rds"))
write.csv(variance_summary, file.path(pca_dir, "pca_variance_summary.csv"), row.names = FALSE)
write.csv(pca_scores, file.path(pca_dir, "pca_scores.csv"), row.names = FALSE)

## Full scree plot
png(file.path(pca_dir, "pca_scree_plot.png"), width = 2000, height = 1500, res = 300)
plot(
  variance_explained,
  type = "b",
  pch = 16,
  xlab = "Principal component",
  ylab = "Proportion of variance explained",
  main = "PCA scree plot"
)
dev.off()

## First 100 PCs
png(file.path(pca_dir, "pca_scree_plot_100PCs.png"), width = 2000, height = 1500, res = 300)
plot(
  variance_explained[1:100],
  type = "b",
  pch = 16,
  xlab = "Principal component",
  ylab = "Proportion of variance explained",
  main = "PCA scree plot (first 100 PCs)"
)
dev.off()

## First 50 PCs
png(file.path(pca_dir, "pca_scree_plot_50PCs.png"), width = 2000, height = 1500, res = 300)
plot(
  variance_explained[1:50],
  type = "b",
  pch = 16,
  xlab = "Principal component",
  ylab = "Proportion of variance explained",
  main = "PCA scree plot (first 50 PCs)"
)
dev.off()

## Find the first PC at which cumulative variance reaches 75%
cum_var <- cumsum(variance_explained)
pc_75 <- which(cum_var >= 0.75)[1]

## Cumulative variance plot, truncated to 0-300 PCs, with 75% threshold marked
png(file.path(pca_dir, "pca_cumulative_variance_plot.png"), width = 2000, height = 1500, res = 300)
plot(
  cum_var,
  type = "l",
  lwd = 2,
  xlab = "Principal component",
  ylab = "Cumulative variance explained",
  main = "PCA cumulative variance",
  xlim = c(0, 300)
)
abline(h = 0.75, lty = 2)
abline(v = pc_75, lty = 3)
points(pc_75, cum_var[pc_75], pch = 16, cex = 1.2)
text(
  pc_75, cum_var[pc_75],
  labels = paste0("PC", pc_75),
  pos = 3, offset = 0.7
)
dev.off()

message("75% of variance reached at PC", pc_75)

message("PCA complete.")
message("Saved to: ", pca_dir)
message("Samples: ", nrow(pca_scores))
message("Features: ", ncol(metab_scaled))
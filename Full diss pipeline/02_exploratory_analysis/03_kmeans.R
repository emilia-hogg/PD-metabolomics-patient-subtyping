# 02_exploratory_analysis/03_kmeans.R

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
if (anyDuplicated(pca_scores$sample_id)) stop("Duplicated sample IDs in PCA scores.")
if (anyDuplicated(umap_coords$sample_id)) stop("Duplicated sample IDs in UMAP coordinates.")

umap_cols <- intersect(c("UMAP1", "UMAP2"), names(umap_coords))
if (length(umap_cols) != 2) stop("UMAP coordinates must contain UMAP1 and UMAP2.")

# Load cohort metadata for site/sex/age/diagnosis overlays if available
if (!file.exists(files$clinical_processed)) {
  stop("Missing clinical processed file: ", files$clinical_processed)
}
clinical_obj <- readRDS(files$clinical_processed)
pd_model <- as.data.frame(clinical_obj$pd_model)

if (!"Anonymised_sampleID" %in% names(pd_model)) {
  stop("pd_model must contain Anonymised_sampleID.")
}

if ("site" %in% names(pd_model)) {
  site_col <- "site"
} else if ("Recruitment_site" %in% names(pd_model)) {
  site_col <- "Recruitment_site"
} else if ("SITE" %in% names(pd_model)) {
  site_col <- "SITE"
} else {
  site_col <- NA_character_
}

if (!is.na(site_col)) {
  site_lookup <- unique(pd_model[, c("Anonymised_sampleID", site_col)])
  names(site_lookup)[2] <- "site"
} else {
  site_lookup <- data.frame(Anonymised_sampleID = pd_model$Anonymised_sampleID, site = "Unknown")
}

# Helper functions
safe_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

run_kmeans_grid <- function(mat, sample_ids, output_dir, label_prefix, k_range = 2:10, nstart = 50) {
  safe_dir(output_dir)

  metrics <- list()
  cluster_results <- list()

  for (k in k_range) {
    set.seed(1)
    km <- kmeans(mat, centers = k, nstart = nstart, iter.max = 100)
    sil <- cluster::silhouette(km$cluster, dist(mat))
    avg_sil <- mean(sil[, "sil_width"])

    metrics[[as.character(k)]] <- data.frame(
      k = k,
      tot_withinss = km$tot.withinss,
      avg_silhouette = avg_sil,
      stringsAsFactors = FALSE
    )

    cluster_df <- data.frame(
      sample_id = sample_ids,
      cluster = factor(km$cluster),
      stringsAsFactors = FALSE
    )

    cluster_results[[as.character(k)]] <- cluster_df
    write.csv(
      cluster_df,
      file.path(output_dir, paste0(label_prefix, "_cluster_labels_k", k, ".csv")),
      row.names = FALSE
    )

    # Site-by-cluster diagnostics if site data are available
    if (!all(is.na(site_lookup$site))) {
      site_cluster_df <- merge(site_lookup, cluster_df, by.x = "Anonymised_sampleID", by.y = "sample_id", all.y = TRUE)
      write.csv(
        site_cluster_df,
        file.path(output_dir, paste0(label_prefix, "_site_cluster_k", k, ".csv")),
        row.names = FALSE
      )

      tab_counts <- table(site_cluster_df$site, site_cluster_df$cluster)
      tab_percent <- prop.table(tab_counts, margin = 1) * 100

      write.csv(
        as.data.frame.matrix(tab_counts),
        file.path(output_dir, paste0(label_prefix, "_site_cluster_counts_k", k, ".csv"))
      )
      write.csv(
        as.data.frame.matrix(round(tab_percent, 1)),
        file.path(output_dir, paste0(label_prefix, "_site_cluster_percent_k", k, ".csv"))
      )

      chisq_result <- suppressWarnings(chisq.test(tab_counts))
      capture.output(
        chisq_result,
        file = file.path(output_dir, paste0(label_prefix, "_site_cluster_chisq_k", k, ".txt"))
      )
    }
  }

  metrics_df <- do.call(rbind, metrics)
  write.csv(metrics_df, file.path(output_dir, paste0(label_prefix, "_metrics.csv")), row.names = FALSE)

  # Plots
  png(file.path(output_dir, paste0(label_prefix, "_elbow_plot.png")), width = 2000, height = 1500, res = 300)
  plot(metrics_df$k, metrics_df$tot_withinss, type = "b", pch = 16,
       xlab = "k", ylab = "Total within-cluster sum of squares",
       main = paste0("Elbow plot: ", label_prefix))
  dev.off()

  png(file.path(output_dir, paste0(label_prefix, "_silhouette_plot.png")), width = 2000, height = 1500, res = 300)
  plot(metrics_df$k, metrics_df$avg_silhouette, type = "b", pch = 16,
       xlab = "k", ylab = "Average silhouette width",
       main = paste0("Silhouette plot: ", label_prefix))
  dev.off()

  invisible(list(metrics = metrics_df, clusters = cluster_results))
}

# PCA clustering
pca_mat <- as.matrix(pca_scores[, grep("^PC[0-9]+$", names(pca_scores)), drop = FALSE])
rownames(pca_mat) <- pca_scores$sample_id

pca_sets <- c(5, 10, 20, 30, 40, 50)
pca_out_root <- file.path(paths$exploratory, "kmeans", "pca")
dir.create(pca_out_root, recursive = TRUE, showWarnings = FALSE)

pca_all_results <- list()
for (n_pc in pca_sets) {
  if (ncol(pca_mat) < n_pc) {
    warning("Skipping PCA set ", n_pc, " because there are only ", ncol(pca_mat), " PCs available.")
    next
  }
  message("Running PCA k-means for first ", n_pc, " PCs")
  out_dir <- file.path(pca_out_root, paste0("PC", sprintf("%02d", n_pc)))
  res <- run_kmeans_grid(
    mat = pca_mat[, seq_len(n_pc), drop = FALSE],
    sample_ids = rownames(pca_mat),
    output_dir = out_dir,
    label_prefix = paste0("PCA_PC", sprintf("%02d", n_pc)),
    k_range = 2:10,
    nstart = 50
  )
  pca_all_results[[paste0("PC", n_pc)]] <- res$metrics
}
write.csv(
  do.call(rbind, lapply(names(pca_all_results), function(nm) cbind(pc_set = nm, pca_all_results[[nm]]))),
  file.path(pca_out_root, "pca_kmeans_all_metrics.csv"),
  row.names = FALSE
)

# UMAP clustering
umap_mat <- as.matrix(umap_coords[, c("UMAP1", "UMAP2")])
rownames(umap_mat) <- umap_coords$sample_id

umap_out_root <- file.path(paths$exploratory, "kmeans", "umap")
dir.create(umap_out_root, recursive = TRUE, showWarnings = FALSE)

message("Running UMAP k-means")
umap_res <- run_kmeans_grid(
  mat = umap_mat,
  sample_ids = rownames(umap_mat),
  output_dir = umap_out_root,
  label_prefix = "UMAP",
  k_range = 2:10,
  nstart = 50
)

saveRDS(pca_all_results, file.path(pca_out_root, "pca_kmeans_all_results.rds"))
saveRDS(umap_res, file.path(umap_out_root, "umap_kmeans_results.rds"))

message("K-means complete.")
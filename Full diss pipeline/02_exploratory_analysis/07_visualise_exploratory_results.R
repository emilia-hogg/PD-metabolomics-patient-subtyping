# 02_exploratory_analysis/07_visualise_exploratory_results.R

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

viz_dir <- file.path(paths$exploratory, "visualisations")
report_dir <- file.path(paths$dissertation_figures, "report_sheet")
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

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

save_plot <- function(plot_obj, path, width = 8, height = 6) {
  ggsave(path, plot = plot_obj, width = width, height = height, dpi = 300)
}

copy_if_exists <- function(src, dst_name = basename(src)) {
  if (!file.exists(src)) {
    warning("Missing file, not copied: ", src)
    return(FALSE)
  }
  dst <- file.path(report_dir, dst_name)
  ok <- file.copy(src, dst, overwrite = TRUE)
  if (!ok) warning("Could not copy: ", src)
  ok
}

detect_site_column <- function(df) {
  candidates <- c("site", "Site", "SITE", "Recruitment_site", "RecruitmentSite", "centre", "Center", "Centre")
  found <- intersect(candidates, names(df))[1]
  if (is.na(found)) NA_character_ else found
}

detect_diag_column <- function(df) {
  candidates <- c("PATIENT_CASE", "diagnostic_group", "Diagnosis", "diagnosis", "group", "Group")
  found <- intersect(candidates, names(df))[1]
  if (is.na(found)) NA_character_ else found
}

make_age_group <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (sum(is.finite(x)) < 10) return(NULL)

  qs <- unique(as.numeric(quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)))
  if (length(qs) < 3) return(NULL)

  cut(
    x,
    breaks = qs,
    include.lowest = TRUE,
    ordered_result = TRUE
  )
}

plot_umap_by_group <- function(df, group_col, title, out_file) {
  if (!group_col %in% names(df)) return(invisible(FALSE))
  vals <- df[[group_col]]
  if (length(unique(na.omit(vals))) < 2) return(invisible(FALSE))

  p <- ggplot(df, aes(x = UMAP1, y = UMAP2, colour = .data[[group_col]])) +
    geom_point(size = 1.1, alpha = 0.75) +
    theme_minimal() +
    labs(title = title, colour = group_col)

  save_plot(p, out_file)
  invisible(TRUE)
}

plot_cluster_umap <- function(umap_df, cluster_df, title, out_file) {
  if (is.null(cluster_df) || !"sample_id" %in% names(cluster_df) || !"cluster" %in% names(cluster_df)) {
    return(invisible(FALSE))
  }

  plot_df <- merge(umap_df, cluster_df, by = "sample_id", all.x = TRUE)
  if (!all(c("UMAP1", "UMAP2") %in% names(plot_df))) {
    return(invisible(FALSE))
  }

  plot_df$cluster <- factor(plot_df$cluster)

  p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
    geom_point(size = 1.1, alpha = 0.75) +
    theme_minimal() +
    labs(title = title, colour = "Cluster")

  save_plot(p, out_file)
  invisible(TRUE)
}

make_site_cluster_outputs <- function(site_cluster_df, out_prefix, out_dir) {
  tab_counts <- table(site_cluster_df$site, site_cluster_df$cluster)
  tab_counts_df <- as.data.frame.matrix(tab_counts)
  tab_counts_df <- tibble::rownames_to_column(tab_counts_df, var = "site")

  tab_percent <- prop.table(tab_counts, margin = 1) * 100
  tab_percent_df <- as.data.frame.matrix(round(tab_percent, 1))
  tab_percent_df <- tibble::rownames_to_column(tab_percent_df, var = "site")

  write.csv(tab_counts_df, file.path(out_dir, paste0(out_prefix, "_site_cluster_counts.csv")), row.names = FALSE)
  write.csv(tab_percent_df, file.path(out_dir, paste0(out_prefix, "_site_cluster_percent.csv")), row.names = FALSE)

  long_df <- as.data.frame(tab_percent)
  names(long_df) <- c("site", "cluster", "percent")

  p_heat <- ggplot(long_df, aes(x = cluster, y = site, fill = percent)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(percent, 1)), size = 3) +
    theme_minimal() +
    labs(title = paste0("Site vs cluster heatmap: ", out_prefix), x = "Cluster", y = "Site", fill = "%")

  save_plot(p_heat, file.path(out_dir, paste0(out_prefix, "_site_cluster_heatmap.png")), width = 9, height = 6)

  if (nrow(tab_counts_df) > 1 && ncol(tab_counts_df) > 1) {
    chisq_res <- tryCatch(chisq.test(tab_counts), error = function(e) e)
    capture.output(chisq_res, file = file.path(out_dir, paste0(out_prefix, "_site_cluster_chisq.txt")))
  }

  invisible(TRUE)
}

# -------------------------------------------------------------------
# Load processed clinical and UMAP data
# -------------------------------------------------------------------
clinical_obj <- safe_read_rds(files$clinical_processed)
if (is.null(clinical_obj) || is.null(clinical_obj$pd_model)) {
  stop("Missing clinical processed object or pd_model.")
}
pd_model <- as.data.frame(clinical_obj$pd_model)

umap_file <- file.path(paths$exploratory, "umap", "umap_coords_clean.csv")
if (!file.exists(umap_file)) {
  stop("Missing UMAP coordinates file: ", umap_file)
}
umap_coords <- read.csv(umap_file, stringsAsFactors = FALSE)

if (!"sample_id" %in% names(umap_coords)) stop("UMAP coordinates must contain sample_id.")
if (!all(c("UMAP1", "UMAP2") %in% names(umap_coords))) stop("UMAP coordinates must contain UMAP1 and UMAP2.")
if (!"Anonymised_sampleID" %in% names(pd_model)) stop("pd_model must contain Anonymised_sampleID.")

pd_model$Anonymised_sampleID <- as.character(pd_model$Anonymised_sampleID)

umap_df <- merge(
  umap_coords,
  pd_model,
  by.x = "sample_id",
  by.y = "Anonymised_sampleID",
  all.x = TRUE
)

# -------------------------------------------------------------------
# 1) UMAP coloured by covariates
# -------------------------------------------------------------------
cov_dir <- file.path(viz_dir, "umap_covariates")
dir.create(cov_dir, recursive = TRUE, showWarnings = FALSE)

sex_col <- intersect(c("GENDER", "Sex", "sex"), names(umap_df))[1]
if (!is.na(sex_col)) {
  umap_df[[sex_col]] <- as.factor(umap_df[[sex_col]])
  plot_umap_by_group(
    umap_df,
    sex_col,
    "UMAP coloured by sex",
    file.path(cov_dir, "umap_by_sex.png")
  )
}

age_col <- intersect(c("AGE_num", "AGE", "age"), names(umap_df))[1]
if (!is.na(age_col)) {
  umap_df$age_group <- make_age_group(umap_df[[age_col]])
  if (!is.null(umap_df$age_group)) {
    plot_umap_by_group(
      umap_df,
      "age_group",
      "UMAP coloured by age group",
      file.path(cov_dir, "umap_by_age_group.png")
    )
  }
}

site_col <- detect_site_column(umap_df)
if (!is.na(site_col)) {
  umap_df[[site_col]] <- as.factor(umap_df[[site_col]])
  plot_umap_by_group(
    umap_df,
    site_col,
    "UMAP coloured by recruitment site",
    file.path(cov_dir, "umap_by_site.png")
  )
}

diag_col <- detect_diag_column(umap_df)
if (!is.na(diag_col)) {
  umap_df[[diag_col]] <- as.factor(umap_df[[diag_col]])
  if (length(unique(na.omit(umap_df[[diag_col]]))) > 1) {
    plot_umap_by_group(
      umap_df,
      diag_col,
      "UMAP coloured by diagnostic group",
      file.path(cov_dir, "umap_by_diagnostic_group.png")
    )
  }
}

# -------------------------------------------------------------------
# 2) K-means visualisations on PCA-based clustering
# -------------------------------------------------------------------
kmeans_pca_root <- file.path(paths$exploratory, "kmeans", "pca")
kmeans_pca_viz <- file.path(viz_dir, "kmeans_pca")
dir.create(kmeans_pca_viz, recursive = TRUE, showWarnings = FALSE)

pc_dirs <- list.dirs(kmeans_pca_root, full.names = TRUE, recursive = FALSE)
if (length(pc_dirs) == 0) {
  warning("No PCA k-means directories found in: ", kmeans_pca_root)
}

for (pc_dir in pc_dirs) {
  metrics_file <- list.files(pc_dir, pattern = "_metrics\\.csv$", full.names = TRUE)
  if (length(metrics_file) == 0) next
  metrics_file <- metrics_file[1]

  metrics_df <- read.csv(metrics_file, stringsAsFactors = FALSE)
  if (!all(c("k", "avg_silhouette") %in% names(metrics_df))) next

  metrics_df <- metrics_df[order(-metrics_df$avg_silhouette, metrics_df$tot_withinss), , drop = FALSE]
  best_k <- metrics_df$k[1]
  selected_ks <- sort(unique(intersect(c(best_k - 1, best_k, best_k + 1), 2:10)))

  pc_name <- basename(pc_dir)
  out_subdir <- file.path(kmeans_pca_viz, pc_name)
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)

  for (k in selected_ks) {
    cluster_file <- file.path(
      pc_dir,
      paste0(sub("_metrics\\.csv$", "", basename(metrics_file)), "_cluster_labels_k", k, ".csv")
    )
    if (!file.exists(cluster_file)) next

    cluster_df <- read.csv(cluster_file, stringsAsFactors = FALSE)

    plot_cluster_umap(
      umap_df,
      cluster_df,
      paste0("UMAP coloured by PCA k-means clusters (", pc_name, ", k = ", k, ")"),
      file.path(out_subdir, paste0(pc_name, "_k", k, "_umap_clusters.png"))
    )

    if (!is.na(site_col)) {
      site_cluster_df <- merge(
        umap_df[, c("sample_id", site_col), drop = FALSE],
        cluster_df,
        by = "sample_id",
        all.y = TRUE
      )
      names(site_cluster_df)[2] <- "site"
      make_site_cluster_outputs(
        site_cluster_df = site_cluster_df,
        out_prefix = paste0(pc_name, "_k", k),
        out_dir = out_subdir
      )
    }
  }
}

# -------------------------------------------------------------------
# 3) K-means visualisations on UMAP-based clustering
# -------------------------------------------------------------------
kmeans_umap_root <- file.path(paths$exploratory, "kmeans", "umap")
kmeans_umap_viz <- file.path(viz_dir, "kmeans_umap")
dir.create(kmeans_umap_viz, recursive = TRUE, showWarnings = FALSE)

umap_metrics_file <- file.path(kmeans_umap_root, "UMAP_metrics.csv")
if (file.exists(umap_metrics_file)) {
  umap_metrics <- read.csv(umap_metrics_file, stringsAsFactors = FALSE)
  if (all(c("k", "avg_silhouette") %in% names(umap_metrics))) {
    umap_metrics <- umap_metrics[order(-umap_metrics$avg_silhouette, umap_metrics$tot_withinss), , drop = FALSE]
    best_k <- umap_metrics$k[1]
    selected_ks <- sort(unique(intersect(c(best_k - 1, best_k, best_k + 1), 2:10)))

    for (k in selected_ks) {
      cluster_file <- file.path(kmeans_umap_root, paste0("UMAP_cluster_labels_k", k, ".csv"))
      if (!file.exists(cluster_file)) next

      cluster_df <- read.csv(cluster_file, stringsAsFactors = FALSE)

      plot_cluster_umap(
        umap_df,
        cluster_df,
        paste0("UMAP coloured by UMAP k-means clusters (k = ", k, ")"),
        file.path(kmeans_umap_viz, paste0("UMAP_k", k, "_clusters.png"))
      )
    }
  }
}

# -------------------------------------------------------------------
# 4) DBSCAN representative plot from best setting
# -------------------------------------------------------------------
dbscan_dir <- file.path(paths$exploratory, "dbscan")
dbscan_sum_file <- file.path(dbscan_dir, "dbscan_parameter_summary.csv")
if (file.exists(dbscan_sum_file)) {
  dbscan_sum <- read.csv(dbscan_sum_file, stringsAsFactors = FALSE)

  if (all(c("eps", "minPts", "n_clusters", "n_noise") %in% names(dbscan_sum))) {
    dbscan_best <- dbscan_sum[order(-dbscan_sum$n_clusters, dbscan_sum$n_noise, dbscan_sum$eps), , drop = FALSE][1, , drop = FALSE]

    eps_val <- dbscan_best$eps
    minPts_val <- dbscan_best$minPts
    cluster_file <- file.path(dbscan_dir, paste0("dbscan_eps_", eps_val, "_minPts_", minPts_val, "_clusters.csv"))

    if (file.exists(cluster_file)) {
      cluster_df <- read.csv(cluster_file, stringsAsFactors = FALSE)
      plot_cluster_umap(
        umap_df,
        cluster_df,
        paste0("UMAP coloured by best DBSCAN clusters (eps = ", eps_val, ", minPts = ", minPts_val, ")"),
        file.path(viz_dir, "dbscan_best_umap_clusters.png")
      )
    }
  }
}

# -------------------------------------------------------------------
# 5) HDBSCAN representative plot from best setting
# -------------------------------------------------------------------
hdbscan_dir <- file.path(paths$exploratory, "hdbscan")
hdbscan_sum_file <- file.path(hdbscan_dir, "hdbscan_all_minPts_summary.csv")
if (file.exists(hdbscan_sum_file)) {
  hdbscan_sum <- read.csv(hdbscan_sum_file, stringsAsFactors = FALSE)

  if (all(c("minPts", "n_clusters", "n_noise", "mean_membership_prob") %in% names(hdbscan_sum))) {
    hdbscan_clean <- unique(hdbscan_sum[, c("minPts", "n_clusters", "n_noise", "mean_membership_prob", "mean_outlier_score"), drop = FALSE])
    hdbscan_best <- hdbscan_clean[order(-hdbscan_clean$n_clusters, hdbscan_clean$n_noise, -hdbscan_clean$mean_membership_prob), , drop = FALSE][1, , drop = FALSE]

    minPts_val <- hdbscan_best$minPts
    cluster_file <- file.path(hdbscan_dir, paste0("hdbscan_minPts_", minPts_val, "_clusters.csv"))

    if (file.exists(cluster_file)) {
      cluster_df <- read.csv(cluster_file, stringsAsFactors = FALSE)
      plot_cluster_umap(
        umap_df,
        cluster_df,
        paste0("UMAP coloured by best HDBSCAN clusters (minPts = ", minPts_val, ")"),
        file.path(viz_dir, "hdbscan_best_umap_clusters.png")
      )
    }
  }
}

# -------------------------------------------------------------------
# 6) Copy the key figures and tables into the final report sheet folder
# -------------------------------------------------------------------
report_map <- c(
  # PCA
  file.path(paths$exploratory, "pca", "pca_scree_plot.png"),
  file.path(paths$exploratory, "pca", "pca_cumulative_variance_plot.png"),
  file.path(paths$exploratory, "pca", "pca_variance_summary.csv"),

  # UMAP
  file.path(cov_dir, "umap_by_sex.png"),
  file.path(cov_dir, "umap_by_age_group.png"),
  file.path(cov_dir, "umap_by_site.png"),
  file.path(cov_dir, "umap_by_diagnostic_group.png"),
  file.path(viz_dir, "umap_scatter.png"),

  # K-means PCA and UMAP
  file.path(kmeans_pca_viz),
  file.path(kmeans_umap_viz),

  # DBSCAN / HDBSCAN
  file.path(viz_dir, "dbscan_best_umap_clusters.png"),
  file.path(viz_dir, "hdbscan_best_umap_clusters.png"),
  file.path(paths$exploratory, "dbscan", "dbscan_knn_distance_plots.png"),
  file.path(paths$exploratory, "dbscan", "dbscan_parameter_summary_plot.png"),
  file.path(paths$exploratory, "hdbscan", "hdbscan_minPts_overview.png")
)

copied_index <- data.frame(
  source = character(0),
  destination = character(0),
  copied = logical(0),
  stringsAsFactors = FALSE
)

for (src in report_map) {
  if (dir.exists(src)) {
    # Copy the whole directory tree into report_sheet
    rel_name <- basename(src)
    dst_dir <- file.path(report_dir, rel_name)
    dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)

    files_to_copy <- list.files(src, recursive = TRUE, full.names = TRUE)
    for (f in files_to_copy) {
      rel <- substring(f, nchar(src) + 2)
      dst <- file.path(dst_dir, rel)
      dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
      ok <- file.copy(f, dst, overwrite = TRUE)
      copied_index <- rbind(
        copied_index,
        data.frame(source = f, destination = dst, copied = ok, stringsAsFactors = FALSE)
      )
    }
  } else {
    dst_name <- basename(src)
    dst <- file.path(report_dir, dst_name)
    ok <- copy_if_exists(src, dst_name)
    copied_index <- rbind(
      copied_index,
      data.frame(source = src, destination = dst, copied = ok, stringsAsFactors = FALSE)
    )
  }
}

write.csv(copied_index, file.path(report_dir, "report_sheet_copy_index.csv"), row.names = FALSE)
writeLines(
  c(
    "Exploratory report sheet created successfully.",
    paste0("Files copied: ", sum(copied_index$copied, na.rm = TRUE)),
    paste0("Files attempted: ", nrow(copied_index))
  ),
  file.path(report_dir, "report_sheet_readme.txt")
)

# -------------------------------------------------------------------
# 7) Final log
# -------------------------------------------------------------------
log_file <- file.path(viz_dir, "visualisation_outputs_log.txt")
produced_files <- list.files(viz_dir, recursive = TRUE, full.names = TRUE)
writeLines(c(
  "Exploratory visualisation script completed.",
  paste0("Total files written: ", length(produced_files)),
  produced_files
), log_file)

message("Exploratory visualisations complete.")
message("Saved outputs to: ", viz_dir)
message("Report sheet copied to: ", report_dir)
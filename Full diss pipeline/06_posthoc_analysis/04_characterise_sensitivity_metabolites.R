# 04_characterise_sensitivity_metabolites.R
# This script characterises metabolite differences between clusters for each
# sensitivity analysis. It performs metabolite-level statistical tests, applies
# BH correction, adds chemical annotations, summarises metabolite classes, and
# produces diagnostic plots, volcano plots, and heatmaps.
#
# Fixed settings:
# - All sensitivity cluster solutions are expected to contain two clusters.
# - Welch and Wilcoxon tests are performed for each metabolite.
# - P-values are adjusted using the Benjamini-Hochberg (BH) method.
# - Significance threshold = 0.05.
# - The Top-20 branches are restricted to the specified 20 metabolites.
#
# Before running:
# - Set `project_root` to the local project directory.
# - The processed analysis data and sensitivity cluster solutions must exist.
# - `chemical_annotation.xlsx` must exist under `06_posthoc_analysis/inputs/`.
#
# Outputs are saved under `paths$posthoc/sensitivity_metabolites/`, with
# separate results for each sensitivity branch and combined summary outputs.

# -----------------------------------------------------------------------------
# 0. Load config, packages, and shared helpers
# -----------------------------------------------------------------------------
project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))
source(file.path(project_root, "06_posthoc_analysis", "00_shared_characterisation_helpers.R"))

# -----------------------------------------------------------------------------
# 1. Required packages
# -----------------------------------------------------------------------------
required_pkgs <- c("ggplot2", "readxl")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required package: ", pkg, call. = FALSE)
  }
}

# -----------------------------------------------------------------------------
# 2. Basic checks and output root
# -----------------------------------------------------------------------------
assert_file_exists(files$analysis_cohort, "analysis_cohort")
assert_file_exists(files$metabolomics_processed, "metabolomics_processed")
assert_file_exists(files$clinical_processed, "clinical_processed")
assert_file_exists(
  file.path(project_root, "06_posthoc_analysis", "inputs", "chemical_annotation.xlsx"),
  "chemical_annotation.xlsx"
)

out_root <- get_posthoc_output_dir("sensitivity_metabolites")
safe_dir(out_root)

# -----------------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------------
sanitize_filename <- function(x) {
  x <- as.character(x)
  x <- gsub("[\\\\/:*?\"<>|]+", "_", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) == 0, "feature", x)
}

find_first_existing <- function(paths) {
  ok <- paths[file.exists(paths)]
  if (length(ok) == 0) return(NA_character_)
  ok[1]
}

safe_shapiro_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L || length(unique(x)) < 3L) return(NA_real_)
  if (length(x) > 5000L) {
    set.seed(1)
    x <- sample(x, 5000L)
  }
  out <- tryCatch(stats::shapiro.test(x)$p.value, error = function(e) NA_real_)
  as.numeric(out)
}

safe_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(NA_real_)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  m <- mean(x)
  mean(((x - m) / s)^3)
}

safe_wilcox <- function(x, g) {
  x <- as.numeric(x)
  g <- as.factor(g)
  if (nlevels(g) != 2L) {
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  out <- tryCatch(stats::wilcox.test(x ~ g, exact = FALSE), error = function(e) NULL)
  if (is.null(out)) {
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  list(statistic = unname(out$statistic), p.value = unname(out$p.value))
}

safe_welch <- function(x, g) {
  x <- as.numeric(x)
  g <- as.factor(g)
  if (nlevels(g) != 2L) {
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  out <- tryCatch(stats::t.test(x ~ g), error = function(e) NULL)
  if (is.null(out)) {
    return(list(statistic = NA_real_, p.value = NA_real_))
  }
  list(statistic = unname(out$statistic), p.value = unname(out$p.value))
}

cohens_d <- function(x, g) {
  x <- as.numeric(x)
  g <- as.factor(g)
  lv <- levels(g)
  if (length(lv) != 2L) return(NA_real_)
  x1 <- x[g == lv[1]]
  x2 <- x[g == lv[2]]
  m1 <- mean(x1, na.rm = TRUE)
  m2 <- mean(x2, na.rm = TRUE)
  s1 <- stats::sd(x1, na.rm = TRUE)
  s2 <- stats::sd(x2, na.rm = TRUE)
  n1 <- sum(is.finite(x1))
  n2 <- sum(is.finite(x2))
  if (n1 + n2 < 3L) return(NA_real_)
  pooled <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(pooled) || pooled == 0) return(NA_real_)
  (m2 - m1) / pooled
}

safe_ci_mean_diff <- function(x, g, conf.level = 0.95) {
  x <- as.numeric(x)
  g <- as.factor(g)
  lv <- levels(g)
  if (length(lv) != 2L) return(c(NA_real_, NA_real_))
  x1 <- x[g == lv[1]]
  x2 <- x[g == lv[2]]
  n1 <- sum(is.finite(x1))
  n2 <- sum(is.finite(x2))
  if (n1 < 2L || n2 < 2L) return(c(NA_real_, NA_real_))
  m1 <- mean(x1, na.rm = TRUE)
  m2 <- mean(x2, na.rm = TRUE)
  s1 <- stats::sd(x1, na.rm = TRUE)
  s2 <- stats::sd(x2, na.rm = TRUE)
  se <- sqrt(s1^2 / n1 + s2^2 / n2)
  if (!is.finite(se) || se == 0) return(c(NA_real_, NA_real_))
  df <- (s1^2 / n1 + s2^2 / n2)^2 / ((s1^2 / n1)^2 / (n1 - 1) + (s2^2 / n2)^2 / (n2 - 1))
  alpha <- 1 - conf.level
  q <- stats::qt(1 - alpha / 2, df = df)
  diff <- m2 - m1
  c(diff - q * se, diff + q * se)
}

make_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

save_png <- function(file, width = 2000, height = 1500, res = 300, expr) {
  grDevices::png(file, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  invisible(file)
}

# Flexible cluster-label extraction:
# - data frame with sample_id + cluster
# - data frame with sample_id + label
# - 2-column data frame with any names
# - 1-column data frame with rownames as sample IDs
# - list with best_k_cluster_labels / cluster_labels / labels / cluster_df
# - named vector
extract_cluster_labels <- function(obj) {
  if (is.list(obj) && !is.data.frame(obj)) {
    candidate_fields <- c("best_k_cluster_labels", "cluster_labels", "labels", "cluster_df")
    for (fld in candidate_fields) {
      if (!is.null(obj[[fld]])) {
        return(extract_cluster_labels(obj[[fld]]))
      }
    }
    if (!is.null(obj$sample_id) && !is.null(obj$cluster)) {
      obj <- data.frame(sample_id = obj$sample_id, cluster = obj$cluster, stringsAsFactors = FALSE)
    } else {
      stop("Could not interpret the sensitivity cluster object.", call. = FALSE)
    }
  }

  if (is.data.frame(obj)) {
    nm <- names(obj)

    if (all(c("sample_id", "cluster") %in% nm)) {
      out <- obj[, c("sample_id", "cluster"), drop = FALSE]
    } else if (all(c("sample_id", "label") %in% nm)) {
      out <- obj[, c("sample_id", "label"), drop = FALSE]
      names(out) <- c("sample_id", "cluster")
    } else if (ncol(obj) == 2L) {
      out <- obj
      names(out) <- c("sample_id", "cluster")
    } else if (ncol(obj) == 1L && !is.null(rownames(obj))) {
      out <- data.frame(
        sample_id = rownames(obj),
        cluster = obj[[1]],
        stringsAsFactors = FALSE
      )
    } else {
      stop("Could not identify sample_id/cluster columns in a labels data frame.", call. = FALSE)
    }

  } else if (is.vector(obj)) {
    if (is.null(names(obj))) {
      stop(
        "Cluster vector must have names equal to sample IDs, or be stored in a data frame.",
        call. = FALSE
      )
    }
    out <- data.frame(
      sample_id = as.character(names(obj)),
      cluster = as.vector(obj),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Could not interpret the sensitivity cluster labels object.", call. = FALSE)
  }

  out$sample_id <- as.character(out$sample_id)
  out$cluster <- suppressWarnings(as.integer(as.character(out$cluster)))

  if (anyDuplicated(out$sample_id)) {
    stop("Cluster labels contain duplicated sample IDs.", call. = FALSE)
  }
  if (anyNA(out$sample_id) || anyNA(out$cluster)) {
    stop("Cluster labels contain missing values.", call. = FALSE)
  }

  out
}

align_cluster_metadata <- function(cluster_obj, sample_order) {
  sample_order <- as.character(sample_order)
  out <- extract_cluster_labels(cluster_obj)
  out <- out[match(sample_order, out$sample_id), , drop = FALSE]

  if (anyNA(out$sample_id) || anyNA(out$cluster)) {
    stop("Failed to align cluster labels to the metabolomics sample order.", call. = FALSE)
  }
  if (!identical(out$sample_id, sample_order)) {
    stop("Aligned cluster labels are not in the expected sample order.", call. = FALSE)
  }
  out
}

read_chemical_annotation <- function() {
  annot_file <- file.path(project_root, "06_posthoc_analysis", "inputs", "chemical_annotation.xlsx")
  assert_file_exists(annot_file, "chemical_annotation.xlsx")
  ann <- readxl::read_excel(annot_file)
  ann <- as.data.frame(ann)

  if (!"CHEMICAL_NAME" %in% names(ann)) {
    stop("chemical_annotation.xlsx must contain CHEMICAL_NAME.", call. = FALSE)
  }
  if (!"SUPER_PATHWAY" %in% names(ann)) {
    ann$SUPER_PATHWAY <- NA_character_
  }
  if (!"SUB_PATHWAY" %in% names(ann)) {
    ann$SUB_PATHWAY <- NA_character_
  }

  ann$CHEMICAL_NAME <- as.character(ann$CHEMICAL_NAME)
  ann$SUPER_PATHWAY <- as.character(ann$SUPER_PATHWAY)
  ann$SUB_PATHWAY <- as.character(ann$SUB_PATHWAY)
  ann <- ann[!duplicated(ann$CHEMICAL_NAME), , drop = FALSE]
  ann
}

add_annotation_columns <- function(feature_table, annotation) {
  feature_table <- as.data.frame(feature_table)
  annotation <- as.data.frame(annotation)

  if (!"feature_name" %in% names(feature_table)) {
    stop("feature_table must contain feature_name.", call. = FALSE)
  }

  idx <- match(feature_table$feature_name, annotation$CHEMICAL_NAME)
  feature_table$SUPER_PATHWAY <- annotation$SUPER_PATHWAY[idx]
  feature_table$SUB_PATHWAY <- annotation$SUB_PATHWAY[idx]

  extra_cols <- setdiff(names(annotation), c("CHEMICAL_NAME", "SUPER_PATHWAY", "SUB_PATHWAY"))
  if (length(extra_cols) > 0) {
    for (nm in extra_cols) {
      feature_table[[nm]] <- annotation[[nm]][idx]
    }
  }

  feature_table
}

plot_volcano <- function(df, out_file, title_text, p_col = "welch_p_adj", effect_col = "mean_diff") {
  df <- as.data.frame(df)
  if (!all(c(p_col, effect_col, "feature_name") %in% names(df))) return(invisible(NULL))

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Missing required package: ggrepel", call. = FALSE)
  }

  pvals <- suppressWarnings(as.numeric(df[[p_col]]))
  effects <- suppressWarnings(as.numeric(df[[effect_col]]))
  df$neg_log10_p <- -log10(pmax(pvals, .Machine$double.xmin))
  df[[effect_col]] <- effects
  df$significance <- ifelse(pvals < 0.05, "BH < 0.05", "Not significant")

sig_idx <- which(pvals < 0.05 & is.finite(pvals) & is.finite(effects))

label_df <- df[integer(0), , drop = FALSE]
if (length(sig_idx) > 0) {
  sig_df <- df[sig_idx, , drop = FALSE]

  # Rank by adjusted p-value first, then effect size
  sig_df <- sig_df[order(sig_df[[p_col]], -abs(sig_df[[effect_col]])), , drop = FALSE]

  # Fewer labels to avoid clutter
  n_label <- min(8, nrow(sig_df))
  label_df <- sig_df[seq_len(n_label), , drop = FALSE]
  label_df$label <- as.character(label_df$feature_name)
}

  save_png(out_file, expr = {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[effect_col]], y = neg_log10_p, colour = significance)) +
      ggplot2::geom_point(alpha = 0.8, size = 1.8) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2) +
      ggplot2::scale_colour_manual(values = c("BH < 0.05" = "red", "Not significant" = "grey60")) +
      ggplot2::labs(
        title = title_text,
        x = "Mean difference",
        y = "-log10(BH adjusted p)",
        colour = "Significance"
      ) +
      ggplot2::theme_bw(base_size = 12)

    if (nrow(label_df) > 0) {
      p <- p + ggrepel::geom_text_repel(
  data = label_df,
  ggplot2::aes(label = label),
  size = 3,
  box.padding = 0.25,
  point.padding = 0.15,
  max.overlaps = Inf,
  min.segment.length = 0,
  seed = 1,
  show.legend = FALSE
)
    }

    print(p)
  })

  invisible(NULL)
}

plot_top_heatmap <- function(metab_scaled, features, cluster_factor, out_file, title_text) {
  features <- intersect(features, colnames(metab_scaled))
  if (length(features) < 2L) return(invisible(NULL))

  mat <- metab_scaled[, features, drop = FALSE]
  mat <- scale(mat)
  mat <- as.matrix(mat)

  ord <- order(as.integer(cluster_factor), names(cluster_factor))
  mat <- mat[ord, , drop = FALSE]

  df <- as.data.frame(as.table(mat))
  names(df) <- c("sample_id", "feature_name", "value")
  df$sample_id <- factor(df$sample_id, levels = rownames(mat))
  df$feature_name <- factor(df$feature_name, levels = colnames(mat))

  save_png(out_file, width = 2400, height = 1800, res = 250, expr = {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = feature_name, y = sample_id, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::labs(title = title_text, x = "Metabolites", y = "Samples", fill = "Scaled") +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
    print(p)
  })

  invisible(NULL)
}

plot_qq <- function(x, out_file, title_text) {
  x <- x[is.finite(x)]
  if (length(x) < 5L) return(invisible(NULL))
  save_png(out_file, expr = {
    qqnorm(x, main = title_text)
    qqline(x)
  })
  invisible(NULL)
}

summarise_classes <- function(feature_table, class_col, metric_col, top_n = 15L) {
  if (!class_col %in% names(feature_table) || !metric_col %in% names(feature_table)) {
    return(data.frame())
  }
  df <- feature_table[!is.na(feature_table[[class_col]]) & is.finite(feature_table[[metric_col]]), , drop = FALSE]
  if (nrow(df) == 0L) return(data.frame())

  agg <- aggregate(
    df[[metric_col]],
    list(class = df[[class_col]]),
    function(z) c(n = length(z), median = stats::median(z, na.rm = TRUE), mean = mean(z, na.rm = TRUE))
  )
  agg <- do.call(data.frame, agg)
  names(agg) <- c("class", paste0(metric_col, c("_n", "_median", "_mean")))
  agg <- agg[order(-agg[[paste0(metric_col, "_median")]]), , drop = FALSE]
  if (!is.null(top_n) && nrow(agg) > top_n) agg <- agg[seq_len(top_n), , drop = FALSE]
  agg
}

# -----------------------------------------------------------------------------
# 4. Load cohort and metabolomics objects
# -----------------------------------------------------------------------------
objs <- load_analysis_objects()
pd_model <- as.data.frame(objs$pd_model)

metab_obj <- readRDS(files$metabolomics_processed)
if (!all(c("metab_pd", "metab_scaled", "sample_ids") %in% names(metab_obj))) {
  stop("metabolomics_processed.rds must contain metab_pd, metab_scaled, and sample_ids.", call. = FALSE)
}

metab_pd <- as.data.frame(metab_obj$metab_pd)
metab_scaled <- as.matrix(metab_obj$metab_scaled)
metab_sample_ids <- as.character(metab_obj$sample_ids)

if (!identical(rownames(metab_pd), metab_sample_ids)) {
  stop("Row order mismatch in metab_pd.", call. = FALSE)
}
if (!identical(rownames(metab_scaled), metab_sample_ids)) {
  stop("Row order mismatch in metab_scaled.", call. = FALSE)
}
if (!"Anonymised_sampleID" %in% names(pd_model)) {
  stop("pd_model must contain Anonymised_sampleID.", call. = FALSE)
}

pd_model$sample_id <- as.character(pd_model$Anonymised_sampleID)
pd_model <- pd_model[match(metab_sample_ids, pd_model$sample_id), , drop = FALSE]
if (!identical(pd_model$sample_id, metab_sample_ids)) {
  stop("Failed to align clinical data to metabolomics sample order.", call. = FALSE)
}

annotation <- read_chemical_annotation()

top20_metabolites <- c(
  "palmitoleamide (16:1)*",
  "myristoleamide (14:1)*",
  "3-methoxytyrosine",
  "3-methoxytyramine sulfate",
  "myristamide (14:0)*",
  "cyclo(leu-pro)",
  "X-21733",
  "m-tyramine sulfate",
  "margaramide (17:0)*",
  "perfluorooctanesulfonate (PFOS)",
  "heptadecenamide (17:1)*",
  "linolenamide (18:3)*",
  "linoleamide (18:2n6)",
  "fibrinopeptide B (1-13)**",
  "dopamine 3-O-sulfate",
  "X-12410",
  "threonate",
  "biliverdin",
  "p-cresol glucuronide*",
  "acetoacetate"
)

# -----------------------------------------------------------------------------
# 5. Branch definitions
# -----------------------------------------------------------------------------
branch_specs <- list(
  clinical_only = list(
    branch_dir = file.path(paths$sensitivity, "clinical_only"),
    result_candidates = c(
      file.path(paths$sensitivity, "clinical_only", "clinical_only_cluster_results.rds"),
      file.path(paths$sensitivity, "clinical_only", "cluster_results.rds"),
      file.path(paths$sensitivity, "clinical_only", "best_k_cluster_labels.rds")
    )
  ),
  no_LEDD = list(
    branch_dir = file.path(paths$sensitivity, "no_LEDD"),
    result_candidates = c(
      file.path(paths$sensitivity, "no_LEDD", "no_ledd_cluster_results.rds"),
      file.path(paths$sensitivity, "no_LEDD", "cluster_results.rds"),
      file.path(paths$sensitivity, "no_LEDD", "best_k_cluster_labels.rds")
    )
  ),
  representation_comparison = list(
    branch_dir = file.path(paths$sensitivity, "representation_comparison", "cluster_solution"),
    result_candidates = c(
      file.path(paths$sensitivity, "representation_comparison", "cluster_solution", "representation_2network_cluster_results.rds"),
      file.path(paths$sensitivity, "representation_comparison", "cluster_solution", "best_k_cluster_labels.rds"),
      file.path(paths$sensitivity, "representation_comparison", "cluster_solution", "cluster_run_record.rds")
    )
  ),
  top20_combined = list(
    branch_dir = file.path(paths$sensitivity, "top20_combined", "cluster_solution"),
    result_candidates = c(
      file.path(paths$sensitivity, "top20_combined", "cluster_solution", "top20_combined_cluster_results.rds"),
      file.path(paths$sensitivity, "top20_combined", "cluster_solution", "best_k_cluster_labels.rds"),
      file.path(paths$sensitivity, "top20_combined", "cluster_solution", "cluster_run_record.rds")
    )
  ),
  top20_23network = list(
    branch_dir = file.path(paths$sensitivity, "top20_23network", "cluster_solution"),
    result_candidates = c(
      file.path(paths$sensitivity, "top20_23network", "cluster_solution", "top20_23network_cluster_results.rds"),
      file.path(paths$sensitivity, "top20_23network", "cluster_solution", "best_k_cluster_labels.rds"),
      file.path(paths$sensitivity, "top20_23network", "cluster_solution", "cluster_run_record.rds")
    )
  )
)

# -----------------------------------------------------------------------------
# 6. Branch loader
# -----------------------------------------------------------------------------
load_branch_solution <- function(branch_name, result_candidates) {
  result_file <- find_first_existing(result_candidates)
  if (is.na(result_file) || !nzchar(result_file)) {
    stop(
      "Could not find sensitivity result file for branch: ", branch_name, "\nLooked in:\n",
      paste(result_candidates, collapse = "\n"),
      call. = FALSE
    )
  }

  obj <- readRDS(result_file)

  if (branch_name == "representation_comparison" && is.list(obj)) {
    if (!is.null(obj$best_k_cluster_labels)) {
      cluster_obj <- obj$best_k_cluster_labels
    } else if (!is.null(obj$cluster_labels)) {
      cluster_obj <- obj$cluster_labels
    } else if (!is.null(obj$labels)) {
      cluster_obj <- obj$labels
    } else if (!is.null(obj$cluster_df)) {
      cluster_obj <- obj$cluster_df
    } else if (!is.null(obj$sample_id) && !is.null(obj$cluster)) {
      cluster_obj <- data.frame(sample_id = obj$sample_id, cluster = obj$cluster, stringsAsFactors = FALSE)
    } else {
      stop(
        "representation_comparison result file was found, but no cluster labels field was recognisable.",
        call. = FALSE
      )
    }
    sample_ids_branch <- if (!is.null(obj$sample_ids)) as.character(obj$sample_ids) else metab_sample_ids
    fused_network <- if (!is.null(obj$fused_network)) as.matrix(obj$fused_network) else NULL

  } else if (is.list(obj) && !is.null(obj$cluster_labels) && !is.null(obj$sample_ids)) {
    cluster_obj <- obj$cluster_labels
    sample_ids_branch <- as.character(obj$sample_ids)
    fused_network <- if (!is.null(obj$fused_network)) as.matrix(obj$fused_network) else NULL

  } else if (is.list(obj) && !is.null(obj$best_k_cluster_labels)) {
    cluster_obj <- obj$best_k_cluster_labels
    sample_ids_branch <- if (!is.null(obj$sample_ids)) as.character(obj$sample_ids) else metab_sample_ids
    fused_network <- if (!is.null(obj$fused_network)) as.matrix(obj$fused_network) else NULL

  } else if (is.data.frame(obj) || (is.vector(obj) && !is.null(names(obj)))) {
    cluster_obj <- obj
    sample_ids_branch <- metab_sample_ids
    fused_network <- NULL

  } else {
    stop(
      "Could not interpret cluster labels for branch: ", branch_name, "\n",
      "Please inspect names(readRDS(result_file)) or the object structure.",
      call. = FALSE
    )
  }

  cluster_df <- align_cluster_metadata(cluster_obj, sample_ids_branch)

  if (!is.null(fused_network)) {
    if (!identical(rownames(fused_network), sample_ids_branch)) {
      fused_network <- reorder_matrix_to_ids(fused_network, sample_ids_branch)
    }
    check_symmetric_matrix(fused_network, paste0(branch_name, " fused network"))
  }

  list(
    branch_name = branch_name,
    result_file = result_file,
    cluster_df = cluster_df,
    sample_ids = sample_ids_branch,
    fused_network = fused_network
  )
}

# -----------------------------------------------------------------------------
# 7. Binary metabolite characterisation
# -----------------------------------------------------------------------------
run_binary_characterisation <- function(metab_pd, metab_scaled, cluster_df, branch_name,
                                        out_table_dir, out_plot_dir, annotation, diag_dir, class_dir) {
  make_dir(out_table_dir)
  make_dir(out_plot_dir)
  make_dir(diag_dir)
  make_dir(class_dir)

  cluster_factor <- factor(cluster_df$cluster, levels = sort(unique(cluster_df$cluster)))
  names(cluster_factor) <- cluster_df$sample_id

  if (nlevels(cluster_factor) != 2L) {
    stop("All sensitivity analyses should be binary, but found ", nlevels(cluster_factor), " clusters.", call. = FALSE)
  }

  feature_names <- intersect(names(metab_pd), colnames(metab_scaled))
  if (length(feature_names) == 0L) {
    stop("No shared metabolite features found between metab_pd and metab_scaled.", call. = FALSE)
  }

  n_feat <- length(feature_names)

  feature_name <- feature_names
  n_obs <- rep(NA_integer_, n_feat)
  cluster_1_mean <- rep(NA_real_, n_feat)
  cluster_2_mean <- rep(NA_real_, n_feat)
  mean_diff <- rep(NA_real_, n_feat)
  ci_low <- rep(NA_real_, n_feat)
  ci_high <- rep(NA_real_, n_feat)
  welch_statistic <- rep(NA_real_, n_feat)
  welch_p <- rep(NA_real_, n_feat)
  wilcox_statistic <- rep(NA_real_, n_feat)
  wilcox_p <- rep(NA_real_, n_feat)
  cohens_d_vec <- rep(NA_real_, n_feat)
  shapiro_p <- rep(NA_real_, n_feat)
  skewness <- rep(NA_real_, n_feat)

  for (i in seq_along(feature_names)) {
    feat <- feature_names[i]
    x <- suppressWarnings(as.numeric(metab_pd[[feat]]))
    g <- cluster_factor[rownames(metab_pd)]

    if (length(x) != length(g)) {
      stop("Length mismatch for feature ", feat, ".", call. = FALSE)
    }

    ok <- is.finite(x) & !is.na(g)
    x <- x[ok]
    g <- droplevels(g[ok])

    n_obs[i] <- length(x)
    shapiro_p[i] <- safe_shapiro_p(x)
    skewness[i] <- safe_skewness(x)

    if (length(x) < 3L || nlevels(g) != 2L) {
      next
    }

    levs <- levels(g)
    cluster_1_mean[i] <- mean(x[g == levs[1]], na.rm = TRUE)
    cluster_2_mean[i] <- mean(x[g == levs[2]], na.rm = TRUE)
    mean_diff[i] <- cluster_2_mean[i] - cluster_1_mean[i]

    ci <- safe_ci_mean_diff(x, g)
    ci_low[i] <- ci[1]
    ci_high[i] <- ci[2]

    wel <- safe_welch(x, g)
    wil <- safe_wilcox(x, g)

    welch_statistic[i] <- wel$statistic
    welch_p[i] <- wel$p.value
    wilcox_statistic[i] <- wil$statistic
    wilcox_p[i] <- wil$p.value
    cohens_d_vec[i] <- cohens_d(x, g)
  }

  results_df <- data.frame(
    feature_name = feature_name,
    n = n_obs,
    cluster_1_mean = cluster_1_mean,
    cluster_2_mean = cluster_2_mean,
    mean_diff = mean_diff,
    ci_low = ci_low,
    ci_high = ci_high,
    welch_statistic = welch_statistic,
    welch_p = welch_p,
    wilcox_statistic = wilcox_statistic,
    wilcox_p = wilcox_p,
    cohens_d = cohens_d_vec,
    shapiro_p = shapiro_p,
    skewness = skewness,
    stringsAsFactors = FALSE
  )

  results_df$welch_p_adj <- stats::p.adjust(results_df$welch_p, method = "BH")
  results_df$wilcox_p_adj <- stats::p.adjust(results_df$wilcox_p, method = "BH")

  numeric_cols <- c(
    "n", "cluster_1_mean", "cluster_2_mean", "mean_diff", "ci_low", "ci_high",
    "welch_statistic", "welch_p", "wilcox_statistic", "wilcox_p",
    "cohens_d", "shapiro_p", "skewness", "welch_p_adj", "wilcox_p_adj"
  )
  for (nm in numeric_cols) {
    results_df[[nm]] <- suppressWarnings(as.numeric(results_df[[nm]]))
  }

  results_annotated <- add_annotation_columns(results_df, annotation)
  results_annotated <- results_annotated[order(results_annotated$welch_p_adj, -abs(results_annotated$mean_diff)), , drop = FALSE]
  rownames(results_annotated) <- NULL

  for (nm in numeric_cols) {
    if (nm %in% names(results_annotated)) {
      results_annotated[[nm]] <- suppressWarnings(as.numeric(results_annotated[[nm]]))
    }
  }

  summary_table <- data.frame(
    metric = c(
      "n_features_total",
      "n_features_with_bh_significant_welch",
      "prop_features_with_bh_significant_welch",
      "median_welch_p_adj",
      "median_abs_mean_diff",
      "median_abs_cohens_d",
      "prop_features_shapiro_p_lt_0.05",
      "prop_features_abs_skewness_gt_1"
    ),
    value = c(
      nrow(results_annotated),
      sum(results_annotated$welch_p_adj < 0.05, na.rm = TRUE),
      mean(results_annotated$welch_p_adj < 0.05, na.rm = TRUE),
      median(results_annotated$welch_p_adj, na.rm = TRUE),
      median(abs(results_annotated$mean_diff), na.rm = TRUE),
      median(abs(results_annotated$cohens_d), na.rm = TRUE),
      mean(results_annotated$shapiro_p < 0.05, na.rm = TRUE),
      mean(abs(results_annotated$skewness) > 1, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  top_features <- head(results_annotated$feature_name, 30)
  top_features <- top_features[top_features %in% colnames(metab_scaled)]

  save_csv(results_annotated, file.path(out_table_dir, paste0(branch_name, "_binary_feature_table.csv")))
  save_csv(summary_table, file.path(out_table_dir, paste0(branch_name, "_binary_summary.csv")))
  saveRDS(results_annotated, file.path(out_table_dir, paste0(branch_name, "_binary_feature_table.rds")))
  saveRDS(summary_table, file.path(out_table_dir, paste0(branch_name, "_binary_summary.rds")))

  plot_volcano(
    results_annotated,
    file.path(out_plot_dir, paste0(branch_name, "_volcano.png")),
    paste0("Binary metabolite volcano plot (", branch_name, ")"),
    p_col = "welch_p_adj",
    effect_col = "mean_diff"
  )

  if (length(top_features) >= 2L) {
    plot_top_heatmap(
      metab_scaled = metab_scaled,
      features = top_features,
      cluster_factor = cluster_factor,
      out_file = file.path(out_plot_dir, paste0(branch_name, "_top_metabolites_heatmap.png")),
      title_text = paste0("Top metabolite differences (", branch_name, ")")
    )
  }

  qq_candidates <- unique(na.omit(c(
    results_annotated$feature_name[1],
    results_annotated$feature_name[which.max(abs(results_annotated$mean_diff))],
    results_annotated$feature_name[which.max(abs(results_annotated$cohens_d))]
  )))
  for (feat in qq_candidates) {
    plot_qq(
      as.numeric(metab_pd[[feat]]),
      file.path(diag_dir, paste0(branch_name, "_qq_", sanitize_filename(feat), ".png")),
      paste0("QQ plot: ", feat, " (", branch_name, ")")
    )
  }

class_summaries <- list(
    super_pathway = summarise_classes_annotated(
      results_annotated, class_col = "SUPER_PATHWAY",
      p_col = "welch_p_adj", diff_col = "mean_diff",
      direction = "lowest", top_n = NULL          # only ~9 classes, keep all
    ),
    sub_pathway_strongest = summarise_classes_annotated(
      results_annotated, class_col = "SUB_PATHWAY",
      p_col = "welch_p_adj", diff_col = "mean_diff",
      direction = "lowest", top_n = 15L
    ),
    sub_pathway_weakest = summarise_classes_annotated(
      results_annotated, class_col = "SUB_PATHWAY",
      p_col = "welch_p_adj", diff_col = "mean_diff",
      direction = "highest", top_n = 15L
    ),
    sub_pathway_all = summarise_classes_annotated(
      results_annotated, class_col = "SUB_PATHWAY",
      p_col = "welch_p_adj", diff_col = "mean_diff",
      direction = "lowest", top_n = NULL
    )
  )

  write.csv(class_summaries$super_pathway,
            file.path(class_dir, paste0(branch_name, "_super_pathway_summary.csv")), row.names = FALSE)
  write.csv(class_summaries$sub_pathway_strongest,
            file.path(class_dir, paste0(branch_name, "_sub_pathway_strongest15.csv")), row.names = FALSE)
  write.csv(class_summaries$sub_pathway_weakest,
            file.path(class_dir, paste0(branch_name, "_sub_pathway_weakest15.csv")), row.names = FALSE)
  write.csv(class_summaries$sub_pathway_all,
            file.path(class_dir, paste0(branch_name, "_sub_pathway_all.csv")), row.names = FALSE)
  saveRDS(class_summaries, file.path(class_dir, paste0(branch_name, "_class_summaries.rds")))

  list(
    branch_name = branch_name,
    cluster_factor = cluster_factor,
    results = results_df,
    results_annotated = results_annotated,
    summary_table = summary_table,
    class_summaries = class_summaries
  )
}

# -----------------------------------------------------------------------------
# 8. Run all sensitivity branches
# -----------------------------------------------------------------------------
branch_results <- list()
branch_overview <- list()
branch_class_summaries <- list()

for (branch_name in names(branch_specs)) {
  message("Running sensitivity metabolite characterisation for branch: ", branch_name)
  spec <- branch_specs[[branch_name]]
  branch <- load_branch_solution(branch_name = branch_name, result_candidates = spec$result_candidates)

  branch_out_dir <- file.path(out_root, branch_name)
  branch_table_dir <- file.path(branch_out_dir, "tables")
  branch_plot_dir <- file.path(branch_out_dir, "plots")
  branch_diag_dir <- file.path(branch_out_dir, "diagnostics")
  branch_class_dir <- file.path(branch_out_dir, "class_summaries")
  for (d in c(branch_out_dir, branch_table_dir, branch_plot_dir, branch_diag_dir, branch_class_dir)) {
    make_dir(d)
  }

  branch_metab_pd <- metab_pd
  branch_metab_scaled <- metab_scaled

  if (branch_name %in% c("top20_combined", "top20_23network")) {
    missing_features <- setdiff(top20_metabolites, colnames(branch_metab_pd))
    if (length(missing_features) > 0L) {
      stop(
        "The following Top-20 metabolites could not be found in the metabolomics matrix:\n",
        paste(missing_features, collapse = "\n"),
        call. = FALSE
      )
    }
    branch_metab_pd <- branch_metab_pd[, top20_metabolites, drop = FALSE]
    branch_metab_scaled <- branch_metab_scaled[, top20_metabolites, drop = FALSE]
  }

  res <- run_binary_characterisation(
    metab_pd = branch_metab_pd,
    metab_scaled = branch_metab_scaled,
    cluster_df = branch$cluster_df,
    branch_name = branch_name,
    out_table_dir = branch_table_dir,
    out_plot_dir = branch_plot_dir,
    annotation = annotation,
    diag_dir = branch_diag_dir,
    class_dir = branch_class_dir
  )

  branch_results[[branch_name]] <- list(branch = branch, results = res)

  branch_overview[[branch_name]] <- data.frame(
    branch_name = branch_name,
    n_clusters = length(unique(branch$cluster_df$cluster)),
    n_features = nrow(res$results_annotated),
    n_bh_significant = sum(res$results_annotated$welch_p_adj < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  branch_class_summaries[[branch_name]] <- res$class_summaries
}

branch_overview_df <- do.call(rbind, branch_overview)
rownames(branch_overview_df) <- NULL
save_csv(branch_overview_df, file.path(out_root, "sensitivity_metabolite_branch_overview.csv"))
saveRDS(branch_overview_df, file.path(out_root, "sensitivity_metabolite_branch_overview.rds"))
saveRDS(branch_results, file.path(out_root, "sensitivity_metabolite_characterisation_results.rds"))
saveRDS(branch_class_summaries, file.path(out_root, "sensitivity_metabolite_class_summaries.rds"))

# -----------------------------------------------------------------------------
# 9. Combined summary
# -----------------------------------------------------------------------------
combined_summary <- do.call(rbind, lapply(names(branch_results), function(nm) {
  x <- branch_results[[nm]]$results$summary_table
  data.frame(branch_name = nm, x, stringsAsFactors = FALSE)
}))
rownames(combined_summary) <- NULL
save_csv(combined_summary, file.path(out_root, "sensitivity_metabolite_combined_summary.csv"))
saveRDS(combined_summary, file.path(out_root, "sensitivity_metabolite_combined_summary.rds"))

# -----------------------------------------------------------------------------
# 10. Session info
# -----------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_root, "session_info.txt"))
message("All sensitivity metabolite characterisation complete.")
message("Outputs saved to: ", out_root)

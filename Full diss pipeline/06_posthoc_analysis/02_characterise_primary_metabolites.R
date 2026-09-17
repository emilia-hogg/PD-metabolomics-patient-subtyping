# 06_posthoc_analysis/02_characterise_primary_metabolites.R
# Re-runnable from a fresh R session.
#
# - load the primary cluster labels
# - merge them onto the metabolomics matrix
# - run per-metabolite Welch tests + Wilcoxon tests
# - calculate mean differences, confidence intervals, Cohen's d
# - apply BH correction
# - generate coloured volcano plots and heatmaps
# - merge the metabolite annotation file with SUPER_PATHWAY, SUB_PATHWAY,
#   and CHEMICAL_NAME
# - summarise metabolite classes
# - save explicit QQ-plot selection reasons
# - stop short of pathway enrichment

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

# -------------------------------------------------------------------
# Extra packages used only in this post hoc script
# -------------------------------------------------------------------
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

invisible(lapply(c("readxl", "pheatmap", "ggrepel"), install_if_missing))

# -------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------
if (!dir.exists(paths$posthoc)) {
  dir.create(paths$posthoc, recursive = TRUE, showWarnings = FALSE)
}

out_root <- file.path(paths$posthoc, "primary_metabolites")
plots_dir <- file.path(out_root, "plots")
tables_dir <- file.path(out_root, "tables")
diag_dir <- file.path(out_root, "diagnostics")
pairwise_dir <- file.path(out_root, "pairwise")
inputs_dir <- file.path(project_root, "06_posthoc_analysis", "inputs")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pairwise_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------
clean_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[nchar(x) == 0] <- NA_character_
  tolower(x)
}

safe_filename <- function(x, max_len = 120) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[/\\\\:*?\"<>|]", "_", x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (nchar(x) > max_len) {
    x <- substr(x, 1, max_len)
  }
  if (!nzchar(x)) {
    x <- "unnamed_feature"
  }
  x
}

parse_cluster_id <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.numeric(x) || is.integer(x)) {
    out <- as.integer(x)
    if (anyNA(out)) {
      stop("Cluster labels contain NA values after numeric conversion.", call. = FALSE)
    }
    return(out)
  }

  x <- as.character(x)
  digits <- stringr::str_extract(x, "-?[0-9]+")
  out <- suppressWarnings(as.integer(digits))

  if (all(!is.na(out))) {
    return(out)
  }

  levs <- unique(x)
  mapped <- setNames(seq_along(levs), levs)
  out2 <- unname(mapped[x])
  if (anyNA(out2)) {
    stop("Some cluster labels could not be converted to integers.", call. = FALSE)
  }
  as.integer(out2)
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(is.finite(x)) < 2) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}

safe_shapiro_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || length(unique(x)) < 3) return(NA_real_)
  if (length(x) > 5000) {
    set.seed(1)
    x <- sample(x, 5000)
  }
  out <- tryCatch(stats::shapiro.test(x)$p.value, error = function(e) NA_real_)
  as.numeric(out)
}

safe_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  m <- mean(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean(((x - m) / s)^3)
}

cohens_d <- function(x1, x2) {
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  n1 <- length(x1)
  n2 <- length(x2)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  s1 <- stats::sd(x1)
  s2 <- stats::sd(x2)
  pooled_var <- (((n1 - 1) * s1^2) + ((n2 - 1) * s2^2)) / (n1 + n2 - 2)
  if (!is.finite(pooled_var) || pooled_var <= 0) return(NA_real_)
  (mean(x1) - mean(x2)) / sqrt(pooled_var)
}

welch_ci <- function(x1, x2, conf_level = 0.95) {
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  n1 <- length(x1)
  n2 <- length(x2)
  if (n1 < 2 || n2 < 2) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  m1 <- mean(x1)
  m2 <- mean(x2)
  s1 <- stats::sd(x1)
  s2 <- stats::sd(x2)
  se <- sqrt((s1^2 / n1) + (s2^2 / n2))
  if (!is.finite(se) || se <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  num <- (s1^2 / n1 + s2^2 / n2)^2
  den <- ((s1^2 / n1)^2 / (n1 - 1)) + ((s2^2 / n2)^2 / (n2 - 1))
  if (!is.finite(den) || den <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  df <- num / den
  if (!is.finite(df) || df <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  alpha <- 1 - conf_level
  crit <- stats::qt(1 - alpha / 2, df = df)
  diff <- m1 - m2
  c(lower = diff - crit * se, upper = diff + crit * se)
}

find_first_existing <- function(paths_vec) {
  existing <- paths_vec[file.exists(paths_vec)]
  if (length(existing) == 0) return(NA_character_)
  existing[1]
}

extract_cluster_labels <- function(obj) {
  if (is.data.frame(obj)) {
    nm <- names(obj)
    sid_col <- nm[tolower(nm) == "sample_id"]
    cl_col <- nm[tolower(nm) %in% c("cluster", "label", "labels", "cluster_id")]

    if (length(sid_col) == 0 || length(cl_col) == 0) {
      if (ncol(obj) >= 2) {
        out <- data.frame(
          sample_id = as.character(obj[[1]]),
          cluster = obj[[2]],
          stringsAsFactors = FALSE
        )
      } else {
        stop("Could not identify sample_id and cluster columns in the labels file.", call. = FALSE)
      }
    } else {
      out <- data.frame(
        sample_id = as.character(obj[[sid_col[1]]]),
        cluster = obj[[cl_col[1]]],
        stringsAsFactors = FALSE
      )
    }
  } else if (is.list(obj) && !is.null(obj$sample_id) && !is.null(obj$cluster)) {
    out <- data.frame(
      sample_id = as.character(obj$sample_id),
      cluster = obj$cluster,
      stringsAsFactors = FALSE
    )
  } else if (is.vector(obj) && !is.null(names(obj))) {
    out <- data.frame(
      sample_id = as.character(names(obj)),
      cluster = as.vector(obj),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Could not interpret the primary cluster labels object.", call. = FALSE)
  }

  out$sample_id <- clean_key(out$sample_id)
  out <- out[!is.na(out$sample_id), , drop = FALSE]
  out$cluster <- parse_cluster_id(out$cluster)

  if (anyDuplicated(out$sample_id)) {
    stop("Primary cluster labels contain duplicated sample IDs.", call. = FALSE)
  }

  out[order(out$sample_id), , drop = FALSE]
}

save_qq_plot <- function(x, file, title_text) {
  x <- x[is.finite(x)]
  png(file, width = 2000, height = 1500, res = 300)
  on.exit(dev.off(), add = TRUE)
  qqnorm(x, main = title_text, pch = 16, cex = 0.55)
  qqline(x, col = "red", lwd = 2)
}

plot_volcano <- function(df, out_file, title_text, xlab_text) {
  df$neglog10_welch_p_adj <- -log10(pmax(df$welch_p_adj, .Machine$double.eps))
  df$sig_class <- dplyr::case_when(
    df$welch_p_adj < 0.05 & df$mean_diff > 0 ~ "Significant positive",
    df$welch_p_adj < 0.05 & df$mean_diff < 0 ~ "Significant negative",
    TRUE ~ "Not significant"
  )
  df$point_size <- pmin(abs(df$mean_diff), stats::quantile(abs(df$mean_diff), 0.95, na.rm = TRUE))

  df$label_text <- ifelse(
    !is.na(df$CHEMICAL_NAME) & nzchar(df$CHEMICAL_NAME),
    df$CHEMICAL_NAME,
    df$feature_name
  )

  top_for_labels <- df[order(df$welch_p_adj, -abs(df$mean_diff)), , drop = FALSE]
  top_for_labels <- top_for_labels[seq_len(min(12L, nrow(top_for_labels))), , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean_diff, y = neglog10_welch_p_adj)) +
    ggplot2::geom_point(ggplot2::aes(color = sig_class, size = point_size), alpha = 0.8) +
    ggplot2::scale_color_manual(
      values = c(
        "Significant positive" = "#D55E00",
        "Significant negative" = "#0072B2",
        "Not significant" = "grey70"
      )
    ) +
    ggplot2::scale_size_continuous(name = "|mean difference|", range = c(1.2, 5)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggplot2::labs(
      title = title_text,
      x = xlab_text,
      y = "-log10(BH-adjusted Welch p)",
      color = "Category"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggrepel::geom_text_repel(
      data = top_for_labels,
      ggplot2::aes(label = label_text),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    )

  ggplot2::ggsave(out_file, plot = p, width = 10.5, height = 7.5, dpi = 300)
  invisible(p)
}

plot_heatmap <- function(mat_samples_by_features, clusters_df, out_file, title_text) {
  clusters_df <- as.data.frame(clusters_df, stringsAsFactors = FALSE)

  if (!all(c("sample_id", "cluster") %in% names(clusters_df))) {
    stop("clusters_df must contain sample_id and cluster columns.", call. = FALSE)
  }

  mat_samples_by_features <- as.matrix(mat_samples_by_features)
  if (is.null(rownames(mat_samples_by_features))) {
    stop("Matrix must have row names equal to sample IDs.", call. = FALSE)
  }

  mat_ids <- rownames(mat_samples_by_features)
  idx <- match(mat_ids, clusters_df$sample_id)

  if (anyNA(idx)) {
    missing_ids <- mat_ids[is.na(idx)]
    stop(
      "Heatmap sample annotations do not match matrix row names. First missing sample_id: ",
      missing_ids[1],
      call. = FALSE
    )
  }

  cl <- clusters_df$cluster[idx]
  cl <- factor(cl, levels = sort(unique(cl)))

  ordered_samples <- mat_ids[order(cl)]
  mat_ord <- mat_samples_by_features[ordered_samples, , drop = FALSE]

  hm <- t(mat_ord)
  hm <- t(scale(t(hm)))
  hm[!is.finite(hm)] <- 0

  anno_col <- data.frame(cluster = factor(cl[match(ordered_samples, mat_ids)], levels = levels(cl)))
  rownames(anno_col) <- ordered_samples

  pheatmap::pheatmap(
    hm,
    annotation_col = anno_col,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_colnames = FALSE,
    fontsize_row = 8,
    main = title_text,
    filename = out_file,
    width = 13,
    height = 9
  )

  invisible(hm)
}

pairwise_metabolite_comparison <- function(mat, clusters, cluster_a, cluster_b) {
  idx_a <- which(clusters == cluster_a)
  idx_b <- which(clusters == cluster_b)

  if (length(idx_a) < 2 || length(idx_b) < 2) {
    stop("Not enough samples in one of the clusters for pairwise testing.", call. = FALSE)
  }

  feature_names <- colnames(mat)

  res_list <- lapply(seq_len(ncol(mat)), function(j) {
    x1 <- as.numeric(mat[idx_a, j])
    x2 <- as.numeric(mat[idx_b, j])
    x1 <- x1[is.finite(x1)]
    x2 <- x2[is.finite(x2)]

    mean_a <- safe_mean(x1)
    mean_b <- safe_mean(x2)
    sd_a <- safe_sd(x1)
    sd_b <- safe_sd(x2)
    mean_diff <- mean_a - mean_b

    tt <- tryCatch(stats::t.test(x1, x2, var.equal = FALSE), error = function(e) NULL)
    wt <- tryCatch(stats::wilcox.test(x1, x2, exact = FALSE), error = function(e) NULL)

    ci <- welch_ci(x1, x2, conf_level = 0.95)
    d <- cohens_d(x1, x2)

    data.frame(
      feature_name = feature_names[j],
      n_cluster_a = length(x1),
      n_cluster_b = length(x2),
      mean_cluster_a = mean_a,
      mean_cluster_b = mean_b,
      sd_cluster_a = sd_a,
      sd_cluster_b = sd_b,
      mean_diff = mean_diff,
      ci_lower = ci[["lower"]],
      ci_upper = ci[["upper"]],
      welch_t = if (is.null(tt)) NA_real_ else unname(tt$statistic),
      welch_df = if (is.null(tt)) NA_real_ else unname(tt$parameter),
      welch_p = if (is.null(tt)) NA_real_ else tt$p.value,
      wilcox_w = if (is.null(wt)) NA_real_ else unname(wt$statistic),
      wilcox_p = if (is.null(wt)) NA_real_ else wt$p.value,
      cohens_d = d,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, res_list)
  out$welch_p_adj <- stats::p.adjust(out$welch_p, method = "BH")
  out$wilcox_p_adj <- stats::p.adjust(out$wilcox_p, method = "BH")
  out
}

summarise_class_table <- function(df, class_col, significant_cutoff = 0.05) {
  if (!class_col %in% names(df)) return(data.frame())

  x <- df[!is.na(df[[class_col]]) & nzchar(as.character(df[[class_col]])), , drop = FALSE]
  if (nrow(x) == 0) return(data.frame())

  grp <- as.character(x[[class_col]])
  split_idx <- split(seq_len(nrow(x)), grp)

  out <- do.call(rbind, lapply(names(split_idx), function(g) {
    sub <- x[split_idx[[g]], , drop = FALSE]
    data.frame(
      class = g,
      n_features = nrow(sub),
      n_significant = sum(sub$welch_p_adj < significant_cutoff, na.rm = TRUE),
      prop_significant = mean(sub$welch_p_adj < significant_cutoff, na.rm = TRUE),
      median_welch_p_adj = median(sub$welch_p_adj, na.rm = TRUE),
      median_wilcox_p_adj = median(sub$wilcox_p_adj, na.rm = TRUE),
      median_abs_mean_diff = median(abs(sub$mean_diff), na.rm = TRUE),
      median_cohens_d = median(sub$cohens_d, na.rm = TRUE),
      min_welch_p_adj = min(sub$welch_p_adj, na.rm = TRUE),
      n_up = sum(sub$mean_diff > 0, na.rm = TRUE),
      n_down = sum(sub$mean_diff < 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  out[order(out$n_significant, decreasing = TRUE), , drop = FALSE]
}

summarise_best_feature_classes <- function(best_feature_df, class_col) {
  if (!class_col %in% names(best_feature_df)) return(data.frame())

  x <- best_feature_df[!is.na(best_feature_df[[class_col]]) & nzchar(as.character(best_feature_df[[class_col]])), , drop = FALSE]
  if (nrow(x) == 0) return(data.frame())

  grp <- as.character(x[[class_col]])
  split_idx <- split(seq_len(nrow(x)), grp)

  out <- do.call(rbind, lapply(names(split_idx), function(g) {
    sub <- x[split_idx[[g]], , drop = FALSE]
    data.frame(
      class = g,
      n_features = nrow(sub),
      n_significant_best = sum(sub$best_welch_p_adj < 0.05, na.rm = TRUE),
      prop_significant_best = mean(sub$best_welch_p_adj < 0.05, na.rm = TRUE),
      median_best_welch_p_adj = median(sub$best_welch_p_adj, na.rm = TRUE),
      median_best_abs_mean_diff = median(abs(sub$best_mean_diff), na.rm = TRUE),
      median_best_cohens_d = median(sub$best_cohens_d, na.rm = TRUE),
      min_best_welch_p_adj = min(sub$best_welch_p_adj, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  out[order(out$n_significant_best, decreasing = TRUE), , drop = FALSE]
}

plot_class_bar <- function(df, x_col, y_col, out_file, title_text, xlab_text, ylab_text, n_top = 15) {
  if (nrow(df) == 0) return(invisible(NULL))
  df <- df[order(df[[y_col]], decreasing = TRUE), , drop = FALSE]
  df <- df[seq_len(min(n_top, nrow(df))), , drop = FALSE]
  df[[x_col]] <- factor(df[[x_col]], levels = rev(df[[x_col]]))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title_text, x = xlab_text, y = ylab_text) +
    ggplot2::theme_bw(base_size = 12)

  ggplot2::ggsave(out_file, plot = p, width = 10, height = 6.5, dpi = 300)
  invisible(p)
}

# -------------------------------------------------------------------
# Load metabolomics data
# -------------------------------------------------------------------
if (!file.exists(files$metabolomics_processed)) {
  stop("Missing metabolomics processed file: ", files$metabolomics_processed, call. = FALSE)
}

metab_obj <- readRDS(files$metabolomics_processed)
if (!all(c("metab_pd", "sample_ids") %in% names(metab_obj))) {
  stop("metabolomics_processed.rds must contain 'metab_pd' and 'sample_ids'.", call. = FALSE)
}

metab_pd <- as.data.frame(metab_obj$metab_pd, check.names = FALSE)
sample_ids <- clean_key(metab_obj$sample_ids)

if (is.null(rownames(metab_pd))) {
  stop("metab_pd must have row names equal to sample IDs.", call. = FALSE)
}

rownames(metab_pd) <- clean_key(rownames(metab_pd))
if (anyDuplicated(rownames(metab_pd))) {
  stop("Metabolomics matrix contains duplicated sample IDs.", call. = FALSE)
}
if (anyNA(metab_pd)) {
  stop("Metabolomics matrix contains NA values.", call. = FALSE)
}
if (!identical(sort(rownames(metab_pd)), sort(sample_ids))) {
  stop("Sample IDs in metabolomics data do not match the expected sample_ids object.", call. = FALSE)
}

# -------------------------------------------------------------------
# Load the primary cluster labels
# -------------------------------------------------------------------
cluster_candidates <- c(
  file.path(paths$snf, "main_four_network", "cluster_solution", "best_k_cluster_labels.rds"),
  file.path(paths$snf, "main_four_network", "cluster_solution", "best_k_cluster_labels.csv"),
  file.path(paths$snf, "main_four_network", "best_k_cluster_labels.rds"),
  file.path(paths$snf, "main_four_network", "best_k_cluster_labels.csv"),
  file.path(paths$snf, "cluster_solution", "best_k_cluster_labels.rds"),
  file.path(paths$snf, "cluster_solution", "best_k_cluster_labels.csv")
)

labels_file <- find_first_existing(cluster_candidates)
if (is.na(labels_file) || !nzchar(labels_file)) {
  stop(
    "Could not find the primary cluster labels file. Looked in:\n",
    paste(cluster_candidates, collapse = "\n"),
    call. = FALSE
  )
}

if (grepl("\\.rds$", labels_file, ignore.case = TRUE)) {
  labels_obj <- readRDS(labels_file)
} else {
  labels_obj <- read.csv(labels_file, stringsAsFactors = FALSE, check.names = FALSE)
}

labels_df <- extract_cluster_labels(labels_obj)

if (anyDuplicated(labels_df$sample_id)) {
  stop("Primary cluster labels contain duplicated sample IDs after import.", call. = FALSE)
}
if (!all(labels_df$sample_id %in% rownames(metab_pd))) {
  missing_ids <- setdiff(labels_df$sample_id, rownames(metab_pd))
  stop(
    "Some primary cluster sample IDs are missing from metabolomics data. First missing ID: ",
    missing_ids[1],
    call. = FALSE
  )
}

metab_pd <- metab_pd[labels_df$sample_id, , drop = FALSE]
if (!identical(rownames(metab_pd), labels_df$sample_id)) {
  stop("Failed to align metabolomics rows to primary cluster order.", call. = FALSE)
}

cluster_membership <- data.frame(
  sample_id = labels_df$sample_id,
  cluster = labels_df$cluster,
  stringsAsFactors = FALSE
)
rownames(cluster_membership) <- cluster_membership$sample_id

write.csv(cluster_membership, file.path(out_root, "primary_cluster_membership.csv"), row.names = FALSE)
saveRDS(cluster_membership, file.path(out_root, "primary_cluster_membership.rds"))

cluster_sizes <- as.data.frame(table(cluster_membership$cluster), stringsAsFactors = FALSE)
names(cluster_sizes) <- c("cluster", "n_samples")
cluster_sizes$prop <- cluster_sizes$n_samples / sum(cluster_sizes$n_samples)
write.csv(cluster_sizes, file.path(tables_dir, "primary_cluster_sizes.csv"), row.names = FALSE)

cluster_vec <- cluster_membership$cluster
names(cluster_vec) <- cluster_membership$sample_id
cluster_levels <- sort(unique(cluster_vec))
n_clusters <- length(cluster_levels)

if (n_clusters < 2) {
  stop("Primary cluster solution contains fewer than two clusters.", call. = FALSE)
}

# -------------------------------------------------------------------
# Load chemical annotation file
# -------------------------------------------------------------------
annotation_candidates <- c(
  file.path(inputs_dir, "chemical_annotation.xlsx"),
  file.path(inputs_dir, "chemical_annotation.xls"),
  file.path(inputs_dir, "chemical_annotation.csv")
)

annotation_file <- find_first_existing(annotation_candidates)
if (is.na(annotation_file) || !nzchar(annotation_file)) {
  stop(
    "Could not find the chemical annotation file. Looked in:\n",
    paste(annotation_candidates, collapse = "\n"),
    call. = FALSE
  )
}

if (grepl("\\.csv$", annotation_file, ignore.case = TRUE)) {
  annotation <- read.csv(annotation_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  annotation <- readxl::read_excel(annotation_file, sheet = 1)
  annotation <- as.data.frame(annotation, stringsAsFactors = FALSE, check.names = FALSE)
}

required_annotation_cols <- c("SUPER_PATHWAY", "SUB_PATHWAY", "CHEMICAL_NAME")
missing_annotation_cols <- setdiff(required_annotation_cols, names(annotation))
if (length(missing_annotation_cols) > 0) {
  stop(
    "Annotation file is missing required columns: ",
    paste(missing_annotation_cols, collapse = ", "),
    call. = FALSE
  )
}

annotation$CHEMICAL_NAME <- as.character(annotation$CHEMICAL_NAME)
annotation$SUPER_PATHWAY <- as.character(annotation$SUPER_PATHWAY)
annotation$SUB_PATHWAY <- as.character(annotation$SUB_PATHWAY)
annotation$CHEMICAL_NAME_key <- clean_key(annotation$CHEMICAL_NAME)
annotation <- annotation[!is.na(annotation$CHEMICAL_NAME_key), , drop = FALSE]
annotation <- annotation[!duplicated(annotation$CHEMICAL_NAME_key), , drop = FALSE]

feature_annotation <- data.frame(
  feature_name = colnames(metab_pd),
  feature_name_key = clean_key(colnames(metab_pd)),
  stringsAsFactors = FALSE
)

feature_annotation <- merge(
  feature_annotation,
  annotation,
  by.x = "feature_name_key",
  by.y = "CHEMICAL_NAME_key",
  all.x = TRUE,
  sort = FALSE
)

feature_annotation$label <- ifelse(
  !is.na(feature_annotation$CHEMICAL_NAME) & nzchar(feature_annotation$CHEMICAL_NAME),
  feature_annotation$CHEMICAL_NAME,
  feature_annotation$feature_name
)

annotation_match_summary <- data.frame(
  metric = c(
    "n_features_total",
    "n_features_matched_to_annotation",
    "prop_features_matched_to_annotation"
  ),
  value = c(
    nrow(feature_annotation),
    sum(!is.na(feature_annotation$SUPER_PATHWAY) | !is.na(feature_annotation$SUB_PATHWAY) | !is.na(feature_annotation$CHEMICAL_NAME)),
    mean(!is.na(feature_annotation$SUPER_PATHWAY) | !is.na(feature_annotation$SUB_PATHWAY) | !is.na(feature_annotation$CHEMICAL_NAME))
  ),
  stringsAsFactors = FALSE
)
write.csv(annotation_match_summary, file.path(tables_dir, "annotation_match_summary.csv"), row.names = FALSE)
write.csv(feature_annotation, file.path(tables_dir, "feature_annotation_lookup.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# Overall metabolite distribution diagnostics
# -------------------------------------------------------------------
metab_mat <- as.matrix(metab_pd)
storage.mode(metab_mat) <- "numeric"

diagnostics_df <- data.frame(
  feature_name = colnames(metab_mat),
  shapiro_p = vapply(seq_len(ncol(metab_mat)), function(j) safe_shapiro_p(metab_mat[, j]), numeric(1)),
  skewness = vapply(seq_len(ncol(metab_mat)), function(j) safe_skewness(metab_mat[, j]), numeric(1)),
  mean = vapply(seq_len(ncol(metab_mat)), function(j) safe_mean(metab_mat[, j]), numeric(1)),
  sd = vapply(seq_len(ncol(metab_mat)), function(j) safe_sd(metab_mat[, j]), numeric(1)),
  stringsAsFactors = FALSE
)

diagnostics_df <- merge(
  diagnostics_df,
  feature_annotation[, c("feature_name", "CHEMICAL_NAME", "SUPER_PATHWAY", "SUB_PATHWAY", "label")],
  by = "feature_name",
  all.x = TRUE,
  sort = FALSE
)

write.csv(diagnostics_df, file.path(tables_dir, "metabolite_distribution_diagnostics.csv"), row.names = FALSE)

diagnostic_summary <- data.frame(
  metric = c(
    "n_metabolites",
    "prop_shapiro_p_below_0.05",
    "median_shapiro_p",
    "median_abs_skewness",
    "prop_abs_skewness_above_1"
  ),
  value = c(
    nrow(diagnostics_df),
    mean(diagnostics_df$shapiro_p < 0.05, na.rm = TRUE),
    median(diagnostics_df$shapiro_p, na.rm = TRUE),
    median(abs(diagnostics_df$skewness), na.rm = TRUE),
    mean(abs(diagnostics_df$skewness) > 1, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write.csv(diagnostic_summary, file.path(tables_dir, "metabolite_distribution_summary.csv"), row.names = FALSE)

png(file.path(diag_dir, "shapiro_p_histogram.png"), width = 2000, height = 1500, res = 300)
hist(
  diagnostics_df$shapiro_p,
  breaks = 40,
  main = "Shapiro-Wilk p-values across metabolites",
  xlab = "Shapiro-Wilk p-value",
  col = "grey80",
  border = "white"
)
abline(v = 0.05, lty = 2)
dev.off()

png(file.path(diag_dir, "skewness_histogram.png"), width = 2000, height = 1500, res = 300)
hist(
  diagnostics_df$skewness,
  breaks = 40,
  main = "Metabolite skewness across features",
  xlab = "Skewness",
  col = "grey80",
  border = "white"
)
abline(v = c(-1, 0, 1), lty = c(2, 1, 2))
dev.off()

# -------------------------------------------------------------------
# Pairwise metabolite comparisons
# -------------------------------------------------------------------
pair_list <- utils::combn(cluster_levels, 2, simplify = FALSE)
all_pair_results <- list()
pair_summaries <- list()

for (pair in pair_list) {
  cluster_a <- pair[1]
  cluster_b <- pair[2]
  pair_label <- paste0("cluster_", cluster_a, "_vs_", cluster_b)
  pair_out_dir <- file.path(pairwise_dir, pair_label)
  dir.create(pair_out_dir, recursive = TRUE, showWarnings = FALSE)

  message("Running pairwise comparison: ", pair_label)

  pair_res <- pairwise_metabolite_comparison(
    mat = metab_mat,
    clusters = cluster_vec[rownames(metab_mat)],
    cluster_a = cluster_a,
    cluster_b = cluster_b
  )

  pair_res <- merge(
    pair_res,
    feature_annotation[, c("feature_name", "CHEMICAL_NAME", "SUPER_PATHWAY", "SUB_PATHWAY", "label")],
    by = "feature_name",
    all.x = TRUE,
    sort = FALSE
  )

  pair_res$pair_label <- pair_label
  pair_res$cluster_a <- cluster_a
  pair_res$cluster_b <- cluster_b

  pair_res <- pair_res[order(pair_res$welch_p_adj, -abs(pair_res$mean_diff)), , drop = FALSE]
  rownames(pair_res) <- NULL

  write.csv(pair_res, file.path(pair_out_dir, paste0(pair_label, "_metabolite_results.csv")), row.names = FALSE)
  saveRDS(pair_res, file.path(pair_out_dir, paste0(pair_label, "_metabolite_results.rds")))

  top_hits <- pair_res[seq_len(min(30L, nrow(pair_res))), , drop = FALSE]
  write.csv(top_hits, file.path(pair_out_dir, paste0(pair_label, "_top_hits.csv")), row.names = FALSE)

  pair_summary <- data.frame(
    pair_label = pair_label,
    cluster_a = cluster_a,
    cluster_b = cluster_b,
    n_samples_a = sum(cluster_vec == cluster_a),
    n_samples_b = sum(cluster_vec == cluster_b),
    n_metabolites = nrow(pair_res),
    n_significant_welch = sum(pair_res$welch_p_adj < 0.05, na.rm = TRUE),
    n_significant_wilcox = sum(pair_res$wilcox_p_adj < 0.05, na.rm = TRUE),
    top_feature = pair_res$feature_name[1],
    top_feature_label = pair_res$label[1],
    top_welch_p_adj = pair_res$welch_p_adj[1],
    top_abs_mean_diff = abs(pair_res$mean_diff[1]),
    stringsAsFactors = FALSE
  )
  pair_summaries[[pair_label]] <- pair_summary

  plot_volcano(
    pair_res,
    out_file = file.path(pair_out_dir, paste0(pair_label, "_volcano.png")),
    title_text = paste0("Primary metabolite comparison: ", pair_label),
    xlab_text = paste0("Mean difference (cluster ", cluster_a, " - cluster ", cluster_b, ")")
  )

  top_for_heatmap <- pair_res$feature_name[seq_len(min(25L, nrow(pair_res)))]
  heatmap_mat <- metab_mat[, top_for_heatmap, drop = FALSE]

  plot_heatmap(
    mat_samples_by_features = heatmap_mat,
    clusters_df = cluster_membership,
    out_file = file.path(pair_out_dir, paste0(pair_label, "_top_metabolites_heatmap.png")),
    title_text = paste0("Top metabolites: ", pair_label)
  )

  super_summary <- summarise_class_table(pair_res, "SUPER_PATHWAY", significant_cutoff = 0.05)
  sub_summary <- summarise_class_table(pair_res, "SUB_PATHWAY", significant_cutoff = 0.05)

  write.csv(
    super_summary,
    file.path(pair_out_dir, paste0(pair_label, "_superpathway_summary.csv")),
    row.names = FALSE
  )
  write.csv(
    sub_summary,
    file.path(pair_out_dir, paste0(pair_label, "_subpathway_summary.csv")),
    row.names = FALSE
  )

  all_pair_results[[pair_label]] <- pair_res
}

pairwise_summary_df <- do.call(rbind, pair_summaries)
write.csv(pairwise_summary_df, file.path(tables_dir, "pairwise_summary.csv"), row.names = FALSE)

all_pairwise_results_df <- do.call(rbind, all_pair_results)
write.csv(all_pairwise_results_df, file.path(tables_dir, "all_pairwise_metabolite_results.csv"), row.names = FALSE)
saveRDS(all_pairwise_results_df, file.path(out_root, "all_pairwise_metabolite_results.rds"))

# -------------------------------------------------------------------
# Global feature-level summary:
# one strongest observed comparison per metabolite across all pairs
# -------------------------------------------------------------------
best_feature_df <- dplyr::group_by(all_pairwise_results_df, .data$feature_name) |>
  dplyr::slice_min(order_by = .data$welch_p_adj, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

# Make the key columns explicit and easy to interpret.
best_feature_df$best_pair_label <- best_feature_df$pair_label
best_feature_df$best_welch_p_adj <- best_feature_df$welch_p_adj
best_feature_df$best_mean_diff <- best_feature_df$mean_diff
best_feature_df$best_cohens_d <- best_feature_df$cohens_d
best_feature_df$best_abs_mean_diff <- abs(best_feature_df$mean_diff)
best_feature_df$best_abs_cohens_d <- abs(best_feature_df$cohens_d)

best_feature_df <- best_feature_df[
  order(best_feature_df$best_welch_p_adj, -best_feature_df$best_abs_mean_diff),
  ,
  drop = FALSE
]
rownames(best_feature_df) <- NULL

write.csv(best_feature_df, file.path(tables_dir, "best_feature_overview.csv"), row.names = FALSE)
saveRDS(best_feature_df, file.path(out_root, "best_feature_overview.rds"))

# -------------------------------------------------------------------
# Overall result summary
# -------------------------------------------------------------------
overall_summary <- data.frame(
  metric = c(
    "n_features_total",
    "n_features_with_bh_significant_best_hit",
    "prop_features_with_bh_significant_best_hit",
    "min_best_welch_p_adj",
    "median_best_welch_p_adj",
    "median_best_abs_mean_diff",
    "median_best_abs_cohens_d"
  ),
  value = c(
    nrow(best_feature_df),
    sum(best_feature_df$best_welch_p_adj < 0.05, na.rm = TRUE),
    mean(best_feature_df$best_welch_p_adj < 0.05, na.rm = TRUE),
    min(best_feature_df$best_welch_p_adj, na.rm = TRUE),
    median(best_feature_df$best_welch_p_adj, na.rm = TRUE),
    median(best_feature_df$best_abs_mean_diff, na.rm = TRUE),
    median(best_feature_df$best_abs_cohens_d, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write.csv(overall_summary, file.path(tables_dir, "overall_metabolite_summary.csv"), row.names = FALSE)

no_signal_note <- if (sum(best_feature_df$best_welch_p_adj < 0.05, na.rm = TRUE) == 0) {
  "No metabolite reached BH-adjusted significance in any pairwise comparison. The downstream class summaries therefore describe the strongest non-significant signals only."
} else {
  "Some metabolites reached BH-adjusted significance in pairwise comparisons. Class summaries still report the strongest signals within each class."
}
writeLines(no_signal_note, con = file.path(tables_dir, "interpretation_note.txt"))

# -------------------------------------------------------------------
# Primary column for the cross-branch SUPER_PATHWAY comparison table
# Mirrors summarise_classes() in 04_characterise_sensitivity_metabolites.R
# so the numbers are directly comparable across branches.
# -------------------------------------------------------------------
format_p_display <- function(p, digits = 3) {
  ifelse(
    is.na(p), "-",
    ifelse(p >= 0.001,
           formatC(p, format = "f", digits = digits),
           formatC(p, format = "e", digits = 1))
  )
}

comparison_class_table <- function(df, class_col, p_col, diff_col, alpha = 0.05) {
  keep <- !is.na(df[[class_col]]) & nzchar(as.character(df[[class_col]])) &
    is.finite(df[[p_col]])
  x <- df[keep, , drop = FALSE]
  if (nrow(x) == 0) return(data.frame())

  sp <- split(seq_len(nrow(x)), as.character(x[[class_col]]))

  out <- do.call(rbind, lapply(names(sp), function(g) {
    sub <- x[sp[[g]], , drop = FALSE]
    n_sig <- sum(sub[[p_col]] < alpha, na.rm = TRUE)
    data.frame(
      super_pathway        = g,
      n_metabolites        = nrow(sub),
      n_significant        = n_sig,
      metabolites_display  = if (n_sig > 0) paste0(nrow(sub), " (", n_sig, ")") else as.character(nrow(sub)),
      median_welch_p_adj   = median(sub[[p_col]], na.rm = TRUE),
      mean_welch_p_adj     = mean(sub[[p_col]], na.rm = TRUE),
      min_welch_p_adj      = min(sub[[p_col]], na.rm = TRUE),
      median_abs_mean_diff = median(abs(sub[[diff_col]]), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  out$primary_display <- format_p_display(out$median_welch_p_adj)
  out[order(out$super_pathway), , drop = FALSE]
}

# n_clusters == 2 -> there is only one pairwise comparison, so best_feature_df
# is exactly the single binary test used by the sensitivity branches.
# n_clusters > 2  -> best_feature_df is a minimum over pairs, which is NOT
# comparable to a single binary test; the per-pair tables below are.
primary_super_comparison <- comparison_class_table(
  as.data.frame(best_feature_df),
  class_col = "SUPER_PATHWAY",
  p_col     = "best_welch_p_adj",
  diff_col  = "best_mean_diff"
)

write.csv(
  primary_super_comparison,
  file.path(tables_dir, "primary_superpathway_comparison_column.csv"),
  row.names = FALSE
)

# Same table computed separately for every pairwise comparison, so you can
# report a single, like-for-like contrast if the primary solution has k > 2.
per_pair_super_comparison <- do.call(rbind, lapply(names(all_pair_results), function(nm) {
  tab <- comparison_class_table(
    as.data.frame(all_pair_results[[nm]]),
    class_col = "SUPER_PATHWAY",
    p_col     = "welch_p_adj",
    diff_col  = "mean_diff"
  )
  if (nrow(tab) == 0) return(NULL)
  cbind(pair_label = nm, tab, stringsAsFactors = FALSE)
}))

write.csv(
  per_pair_super_comparison,
  file.path(tables_dir, "primary_superpathway_comparison_by_pair.csv"),
  row.names = FALSE
)

message("Primary superpathway comparison column written to: ",
        file.path(tables_dir, "primary_superpathway_comparison_column.csv"))

# -------------------------------------------------------------------
# Primary SUB_PATHWAY summaries, matching the sensitivity branches
# -------------------------------------------------------------------
primary_sub_strongest <- summarise_classes_annotated(
  as.data.frame(best_feature_df),
  class_col = "SUB_PATHWAY", p_col = "best_welch_p_adj",
  diff_col  = "best_mean_diff", direction = "lowest", top_n = 15L
)
primary_sub_weakest <- summarise_classes_annotated(
  as.data.frame(best_feature_df),
  class_col = "SUB_PATHWAY", p_col = "best_welch_p_adj",
  diff_col  = "best_mean_diff", direction = "highest", top_n = 15L
)
primary_sub_all <- summarise_classes_annotated(
  as.data.frame(best_feature_df),
  class_col = "SUB_PATHWAY", p_col = "best_welch_p_adj",
  diff_col  = "best_mean_diff", direction = "lowest", top_n = NULL
)

write.csv(primary_sub_strongest, file.path(tables_dir, "primary_sub_pathway_strongest15.csv"), row.names = FALSE)
write.csv(primary_sub_weakest,   file.path(tables_dir, "primary_sub_pathway_weakest15.csv"),   row.names = FALSE)
write.csv(primary_sub_all,       file.path(tables_dir, "primary_sub_pathway_all.csv"),         row.names = FALSE)

message("Primary sub-pathway summaries written to: ", tables_dir)

# -------------------------------------------------------------------
# Global superpathway / subpathway summaries:
# strongest non-significant signals only
# -------------------------------------------------------------------
summarise_best_feature_classes <- function(best_feature_df, class_col) {
  if (!class_col %in% names(best_feature_df)) return(data.frame())

  x <- best_feature_df[
    !is.na(best_feature_df[[class_col]]) & nzchar(as.character(best_feature_df[[class_col]])),
    ,
    drop = FALSE
  ]
  if (nrow(x) == 0) return(data.frame())

  grp <- as.character(x[[class_col]])
  split_idx <- split(seq_len(nrow(x)), grp)

  out <- do.call(rbind, lapply(names(split_idx), function(g) {
    sub <- x[split_idx[[g]], , drop = FALSE]
    top_idx <- which.min(sub$best_welch_p_adj)

    data.frame(
      class = g,
      n_features = nrow(sub),
      n_bh_significant = sum(sub$best_welch_p_adj < 0.05, na.rm = TRUE),
      prop_bh_significant = mean(sub$best_welch_p_adj < 0.05, na.rm = TRUE),
      median_best_welch_p_adj = median(sub$best_welch_p_adj, na.rm = TRUE),
      median_best_abs_mean_diff = median(sub$best_abs_mean_diff, na.rm = TRUE),
      median_best_abs_cohens_d = median(sub$best_abs_cohens_d, na.rm = TRUE),
      min_best_welch_p_adj = min(sub$best_welch_p_adj, na.rm = TRUE),
      top_feature = sub$feature_name[top_idx],
      top_feature_label = sub$label[top_idx],
      top_feature_welch_p_adj = sub$best_welch_p_adj[top_idx],
      top_feature_abs_mean_diff = sub$best_abs_mean_diff[top_idx],
      top_feature_abs_cohens_d = sub$best_abs_cohens_d[top_idx],
      stringsAsFactors = FALSE
    )
  }))

  out[order(out$median_best_abs_mean_diff, decreasing = TRUE), , drop = FALSE]
}

best_super_summary <- summarise_best_feature_classes(best_feature_df, "SUPER_PATHWAY")
best_sub_summary <- summarise_best_feature_classes(best_feature_df, "SUB_PATHWAY")

write.csv(best_super_summary, file.path(tables_dir, "best_superpathway_summary.csv"), row.names = FALSE)
write.csv(best_sub_summary, file.path(tables_dir, "best_subpathway_summary.csv"), row.names = FALSE)

# These plots are now explicitly about the strongest non-significant signals.
plot_class_bar(
  best_super_summary,
  x_col = "class",
  y_col = "median_best_abs_mean_diff",
  out_file = file.path(plots_dir, "top_superpathways_by_median_abs_mean_diff.png"),
  title_text = "Superpathways ranked by median absolute mean difference",
  xlab_text = "Super pathway",
  ylab_text = "Median absolute mean difference",
  n_top = 15
)

plot_class_bar(
  best_sub_summary,
  x_col = "class",
  y_col = "median_best_abs_mean_diff",
  out_file = file.path(plots_dir, "top_subpathways_by_median_abs_mean_diff.png"),
  title_text = "Subpathways ranked by median absolute mean difference",
  xlab_text = "Sub pathway",
  ylab_text = "Median absolute mean difference",
  n_top = 15
)

# -------------------------------------------------------------------
# QQ plot selection:
# representative diagnostics, not drivers of inference
# -------------------------------------------------------------------
best_pair_label <- pairwise_summary_df$pair_label[which.min(pairwise_summary_df$top_welch_p_adj)]
best_pair_results <- all_pair_results[[best_pair_label]]

qq_candidates <- data.frame(
  feature_name = character(),
  reason = character(),
  stringsAsFactors = FALSE
)

if (nrow(best_pair_results) > 0) {
  # Top features from the strongest observed pairwise comparison
  top_pair_features <- unique(na.omit(
    best_pair_results$feature_name[seq_len(min(3L, nrow(best_pair_results)))]
  ))

  if (length(top_pair_features) > 0) {
    qq_candidates <- rbind(
      qq_candidates,
      data.frame(
        feature_name = top_pair_features,
        reason = paste0(
          "Top features from the strongest observed pairwise comparison: ",
          best_pair_label
        ),
        stringsAsFactors = FALSE
      )
    )
  }

  # Feature with the largest absolute mean difference in the best pair
  idx_mean <- which.max(abs(best_pair_results$mean_diff))
  if (length(idx_mean) == 1 && is.finite(idx_mean) && nrow(best_pair_results) >= idx_mean) {
    qq_candidates <- rbind(
      qq_candidates,
      data.frame(
        feature_name = best_pair_results$feature_name[idx_mean],
        reason = paste0(
          "Largest absolute mean difference in the strongest observed pairwise comparison: ",
          best_pair_label
        ),
        stringsAsFactors = FALSE
      )
    )
  }

  # Feature with the largest absolute effect size in the best pair
  idx_d <- which.max(abs(best_pair_results$cohens_d))
  if (length(idx_d) == 1 && is.finite(idx_d) && nrow(best_pair_results) >= idx_d) {
    qq_candidates <- rbind(
      qq_candidates,
      data.frame(
        feature_name = best_pair_results$feature_name[idx_d],
        reason = paste0(
          "Largest absolute effect size in the strongest observed pairwise comparison: ",
          best_pair_label
        ),
        stringsAsFactors = FALSE
      )
    )
  }
}

# Add global distribution diagnostics
qq_candidates <- rbind(
  qq_candidates,
  data.frame(
    feature_name = diagnostics_df$feature_name[which.max(abs(diagnostics_df$skewness))],
    reason = "Largest absolute skewness across all metabolites",
    stringsAsFactors = FALSE
  ),
  data.frame(
    feature_name = diagnostics_df$feature_name[which.min(diagnostics_df$shapiro_p)],
    reason = "Smallest Shapiro-Wilk p-value across all metabolites",
    stringsAsFactors = FALSE
  )
)

qq_selection_df <- qq_candidates |>
  dplyr::group_by(.data$feature_name) |>
  dplyr::summarise(reason = paste(unique(.data$reason), collapse = " | "), .groups = "drop")

qq_selection_df <- merge(
  qq_selection_df,
  feature_annotation[, c("feature_name", "CHEMICAL_NAME", "SUPER_PATHWAY", "SUB_PATHWAY", "label")],
  by = "feature_name",
  all.x = TRUE,
  sort = FALSE
)

write.csv(qq_selection_df, file.path(tables_dir, "qq_selection_summary.csv"), row.names = FALSE)

representative_features <- unique(qq_selection_df$feature_name)
representative_features <- representative_features[representative_features %in% colnames(metab_mat)]

for (feat in representative_features) {
  safe_feat <- safe_filename(feat)
  out_file <- file.path(diag_dir, paste0("qq_", safe_feat, ".png"))
  save_qq_plot(
    x = metab_mat[, feat],
    file = out_file,
    title_text = paste0("QQ plot: ", feat)
  )
}

# -------------------------------------------------------------------
# Final pairwise overview table
# -------------------------------------------------------------------
pairwise_best_hit_overview <- do.call(
  rbind,
  lapply(names(all_pair_results), function(nm) {
    df <- all_pair_results[[nm]]
    data.frame(
      pair_label = nm,
      feature_name = df$feature_name[1],
      label = df$label[1],
      welch_p_adj = df$welch_p_adj[1],
      mean_diff = df$mean_diff[1],
      cohens_d = df$cohens_d[1],
      stringsAsFactors = FALSE
    )
  })
)
pairwise_best_hit_overview <- pairwise_best_hit_overview[
  order(pairwise_best_hit_overview$welch_p_adj),
  ,
  drop = FALSE
]
write.csv(pairwise_best_hit_overview, file.path(tables_dir, "pairwise_best_hit_overview.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# Table A12: distributional diagnostics and test concordance
# -------------------------------------------------------------------
 
n_metab <- nrow(diagnostics_df)
 
# --- normality and skew, with counts as well as proportions --------
n_shapiro_sig <- sum(diagnostics_df$shapiro_p < 0.05, na.rm = TRUE)
n_skew_gt1    <- sum(abs(diagnostics_df$skewness) > 1, na.rm = TRUE)
 
idx_min_shapiro <- which.min(diagnostics_df$shapiro_p)
idx_max_skew    <- which.max(abs(diagnostics_df$skewness))
 
lab_of <- function(i) {
  if (length(i) != 1 || is.na(i)) return(NA_character_)
  lb <- diagnostics_df$label[i]
  if (is.na(lb) || !nzchar(lb)) diagnostics_df$feature_name[i] else lb
}
 
# --- Welch vs Wilcoxon concordance, computed per pairwise contrast --
concordance_for <- function(df) {
  ok <- is.finite(df$welch_p_adj) & is.finite(df$wilcox_p_adj)
  w  <- df$welch_p_adj[ok]
  u  <- df$wilcox_p_adj[ok]
 
  sig_w <- w < 0.05
  sig_u <- u < 0.05
 
  data.frame(
    n_tested                = sum(ok),
    n_sig_welch             = sum(sig_w),
    n_sig_wilcox            = sum(sig_u),
    n_sig_both              = sum(sig_w & sig_u),
    n_sig_welch_only        = sum(sig_w & !sig_u),
    n_sig_wilcox_only       = sum(!sig_w & sig_u),
    n_disagree              = sum(sig_w != sig_u),
    prop_disagree           = mean(sig_w != sig_u),
    median_abs_diff_p_adj   = median(abs(w - u)),
    max_abs_diff_p_adj      = max(abs(w - u)),
    spearman_rho            = suppressWarnings(cor(w, u, method = "spearman")),
    spearman_rho_rank_pval  = suppressWarnings(
                                cor(rank(w), rank(u), method = "pearson")),
    stringsAsFactors = FALSE
  )
}
 
conc_by_pair <- do.call(rbind, lapply(names(all_pair_results), function(nm) {
  cbind(pair_label = nm, concordance_for(all_pair_results[[nm]]),
        stringsAsFactors = FALSE)
}))
write.csv(conc_by_pair,
          file.path(tables_dir, "welch_wilcox_concordance_by_pair.csv"),
          row.names = FALSE)
 
# For a two-cluster solution there is a single contrast, which is the
# like-for-like equivalent of each sensitivity branch. For k > 2 the
# pooled row below is a summary across contrasts and the per-pair file
# above is the one to report.
conc <- concordance_for(as.data.frame(all_pairwise_results_df))
 
# --- assemble the table -------------------------------------------
fmt_np <- function(n, N) sprintf("%s (%.1f)", format(n, big.mark = ","), 100 * n / N)
 
table_a12 <- data.frame(
  diagnostic = c(
    "Metabolites tested, n",
    "Metabolites with Shapiro-Wilk p < 0.05, n (%)",
    "Metabolites with absolute skewness > 1, n (%)",
    "Median Shapiro-Wilk p-value",
    "Smallest Shapiro-Wilk p-value",
    "Median absolute skewness",
    "Largest absolute skewness",
    "Metabolites significant by Welch, n",
    "Metabolites significant by Wilcoxon, n",
    "Metabolites where Welch and Wilcoxon disagree, n",
    "Median absolute difference in adjusted p-values",
    "Maximum absolute difference in adjusted p-values",
    "Spearman correlation of adjusted p-values"
  ),
  value = c(
    format(n_metab, big.mark = ","),
    fmt_np(n_shapiro_sig, n_metab),
    fmt_np(n_skew_gt1, n_metab),
    format_p_display(median(diagnostics_df$shapiro_p, na.rm = TRUE)),
    format_p_display(min(diagnostics_df$shapiro_p, na.rm = TRUE)),
    formatC(median(abs(diagnostics_df$skewness), na.rm = TRUE), format = "f", digits = 2),
    formatC(max(abs(diagnostics_df$skewness), na.rm = TRUE), format = "f", digits = 2),
    format(conc$n_sig_welch, big.mark = ","),
    format(conc$n_sig_wilcox, big.mark = ","),
    format(conc$n_disagree, big.mark = ","),
    formatC(conc$median_abs_diff_p_adj, format = "f", digits = 3),
    formatC(conc$max_abs_diff_p_adj, format = "f", digits = 3),
    formatC(conc$spearman_rho, format = "f", digits = 3)
  ),
  notes = c(
    "after zero-variance filtering",
    "",
    "",
    "",
    lab_of(idx_min_shapiro),
    "",
    lab_of(idx_max_skew),
    "BH-adjusted, alpha = 0.05",
    "BH-adjusted, alpha = 0.05",
    sprintf("%d Welch only, %d Wilcoxon only",
            conc$n_sig_welch_only, conc$n_sig_wilcox_only),
    "",
    "",
    "across all metabolites"
  ),
  stringsAsFactors = FALSE
)
 
write.csv(table_a12,
          file.path(tables_dir, "table_A12_distributional_diagnostics.csv"),
          row.names = FALSE)
 
cat("\n--- Table A12 ---\n")
for (i in seq_len(nrow(table_a12))) {
  cat(sprintf("  %-48s %-14s %s\n",
              table_a12$diagnostic[i], table_a12$value[i], table_a12$notes[i]))
}
 
if (conc$n_sig_welch == 0 && conc$n_sig_wilcox == 0) {
  cat("\nNOTE: neither test returned a significant result, so the count of\n")
  cat("disagreements is zero by construction and does not demonstrate\n")
  cat("agreement between the two tests. Quote the Spearman correlation\n")
  cat("above, and the disagreement counts from the sensitivity branches,\n")
  cat("when supporting the claim that the two tests were concordant.\n")
}
 

# -------------------------------------------------------------------
# Final messages
# -------------------------------------------------------------------
message("Primary metabolite characterisation complete.")
message("Outputs saved to: ", out_root)
message("Best pair for representative plots: ", best_pair_label)
message("Number of primary clusters: ", n_clusters)
message(no_signal_note)
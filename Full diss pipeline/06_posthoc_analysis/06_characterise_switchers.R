# 06_posthoc_analysis/06_characterise_switchers.R
# Re-runnable from a fresh R session.
#
# - characterise patients whose cluster assignments change between the primary
#   integrated SNF solution and each sensitivity branch
# - compare switchers vs non-switchers clinically
# - compare switchers vs non-switchers metabolomically
# - assess whether switchers look like borderline / unstable cases

# -----------------------------------------------------------------------------
# 0. Load config, packages, and shared helpers
# -----------------------------------------------------------------------------
project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))
source(file.path(project_root, "06_posthoc_analysis", "00_shared_characterisation_helpers.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Missing required package: ggplot2", call. = FALSE)
}

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# -----------------------------------------------------------------------------
# 1. Utilities
# -----------------------------------------------------------------------------
safe_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

assert_file_exists <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    stop("Missing required file for ", label, ": ", path, call. = FALSE)
  }
  invisible(TRUE)
}

clean_sample_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[nchar(x) == 0] <- NA_character_
  x
}

parse_cluster_id <- function(x) {
  if (is.factor(x)) x <- as.character(x)

  if (is.numeric(x) || is.integer(x)) {
    out <- as.integer(x)
    if (anyNA(out)) stop("Cluster labels contain NA after numeric conversion.", call. = FALSE)
    return(out)
  }

  x <- as.character(x)
  if (length(x) == 0) return(integer(0))

  extracted <- suppressWarnings(as.integer(gsub("^.*?(-?[0-9]+).*$", "\\1", x)))
  if (!anyNA(extracted)) return(as.integer(extracted))

  levs <- unique(x)
  map <- setNames(seq_along(levs), levs)
  out <- unname(map[x])
  if (anyNA(out)) stop("Some cluster labels could not be converted to integers.", call. = FALSE)
  as.integer(out)
}

standardise_cluster_df <- function(obj, fallback_sample_ids = NULL, label = "cluster object") {
  if (is.data.frame(obj)) {
    nm <- names(obj)
    if (all(c("sample_id", "cluster") %in% nm)) {
      out <- obj[, c("sample_id", "cluster"), drop = FALSE]
    } else if (all(c("sample_id", "label") %in% nm)) {
      out <- obj[, c("sample_id", "label"), drop = FALSE]
      names(out)[2] <- "cluster"
    } else if (ncol(obj) == 2) {
      out <- obj
      names(out) <- c("sample_id", "cluster")
    } else {
      stop(label, " does not have obvious sample_id/cluster columns.", call. = FALSE)
    }

    out$sample_id <- clean_sample_ids(out$sample_id)
    out$cluster <- parse_cluster_id(out$cluster)
    out <- out[!is.na(out$sample_id), , drop = FALSE]

    if (anyDuplicated(out$sample_id)) {
      stop(label, " contains duplicated sample IDs.", call. = FALSE)
    }
    return(out)
  }

  if (is.list(obj) && !is.null(obj$sample_id) && !is.null(obj$cluster)) {
    out <- data.frame(
      sample_id = clean_sample_ids(obj$sample_id),
      cluster = parse_cluster_id(obj$cluster),
      stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$sample_id), , drop = FALSE]
    if (anyDuplicated(out$sample_id)) {
      stop(label, " contains duplicated sample IDs.", call. = FALSE)
    }
    return(out)
  }

  if ((is.vector(obj) || is.factor(obj)) && !is.null(names(obj))) {
    out <- data.frame(
      sample_id = clean_sample_ids(names(obj)),
      cluster = parse_cluster_id(obj),
      stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$sample_id), , drop = FALSE]
    if (anyDuplicated(out$sample_id)) {
      stop(label, " contains duplicated sample IDs.", call. = FALSE)
    }
    return(out)
  }

  if (!is.null(fallback_sample_ids) && length(fallback_sample_ids) == length(obj)) {
    out <- data.frame(
      sample_id = clean_sample_ids(fallback_sample_ids),
      cluster = parse_cluster_id(obj),
      stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$sample_id), , drop = FALSE]
    if (anyDuplicated(out$sample_id)) {
      stop(label, " contains duplicated sample IDs.", call. = FALSE)
    }
    return(out)
  }

  stop("Could not standardise ", label, ".", call. = FALSE)
}

find_first_existing <- function(paths_vec) {
  hit <- paths_vec[file.exists(paths_vec)][1]
  if (length(hit) == 0L || is.na(hit) || !nzchar(hit)) return(NA_character_)
  hit
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

cohens_d <- function(x, g) {
  x <- as.numeric(x)
  g <- factor(g)
  lv <- levels(g)
  if (length(lv) != 2) return(NA_real_)
  x1 <- x[g == lv[1]]
  x2 <- x[g == lv[2]]
  m1 <- mean(x1, na.rm = TRUE)
  m2 <- mean(x2, na.rm = TRUE)
  s1 <- stats::sd(x1, na.rm = TRUE)
  s2 <- stats::sd(x2, na.rm = TRUE)
  n1 <- sum(is.finite(x1))
  n2 <- sum(is.finite(x2))
  pooled <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(pooled) || pooled == 0) return(NA_real_)
  (m2 - m1) / pooled
}

summarise_binary_group <- function(df, group_col, vars) {
  out <- list()
  for (v in vars) {
    if (!is.numeric(df[[v]])) next
    tmp <- aggregate(
      df[[v]],
      list(group = df[[group_col]]),
      function(z) c(
        n = sum(is.finite(z)),
        mean = mean(z, na.rm = TRUE),
        sd = stats::sd(z, na.rm = TRUE),
        median = stats::median(z, na.rm = TRUE),
        q25 = stats::quantile(z, 0.25, na.rm = TRUE, names = FALSE),
        q75 = stats::quantile(z, 0.75, na.rm = TRUE, names = FALSE)
      )
    )
    tmp <- do.call(data.frame, tmp)
    tmp$variable <- v
    out[[v]] <- tmp
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

# Simple plotting helpers
plot_boxplot <- function(df, x, y, out_file, title, ylab = NULL) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[[x]])) +
    ggplot2::geom_boxplot(outlier.alpha = 0.2, width = 0.6) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "none") +
    ggplot2::labs(x = NULL, y = ylab %||% y, title = title)
  ggplot2::ggsave(out_file, p, width = 6.2, height = 4.6, dpi = 300)
  invisible(p)
}

plot_volcano <- function(df, out_file, title) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean_diff, y = neglog10_p)) +
    ggplot2::geom_point(ggplot2::aes(color = sig), alpha = 0.7, size = 1.5) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "bottom") +
    ggplot2::labs(
      x = "Mean difference (switcher - non-switcher)",
      y = "-log10 BH-adjusted p",
      color = NULL,
      title = title
    )
  ggplot2::ggsave(out_file, p, width = 7.0, height = 5.4, dpi = 300)
  invisible(p)
}

plot_top_metab_heatmap <- function(df, metabolite_cols, out_file, title) {
  long_df <- do.call(rbind, lapply(metabolite_cols, function(m) {
    data.frame(
      sample_id = df$sample_id,
      switch_group = df$switch_group,
      metabolite = m,
      value = scale(df[[m]])[, 1],
      stringsAsFactors = FALSE
    )
  }))
  long_df$metabolite <- factor(long_df$metabolite, levels = rev(metabolite_cols))
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = metabolite, y = switch_group, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(face = "bold")) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Z", title = title)
  ggplot2::ggsave(out_file, p, width = 10.0, height = 4.4, dpi = 300)
  invisible(p)
}

# -----------------------------------------------------------------------------
# 2. Load analysis objects
# -----------------------------------------------------------------------------
objs <- load_analysis_objects()
pd_model <- as.data.frame(objs$pd_model, check.names = FALSE)
master_ids <- as.character(objs$sample_ids)

if (!"Anonymised_sampleID" %in% names(pd_model)) {
  stop("pd_model must contain Anonymised_sampleID.", call. = FALSE)
}
pd_model$sample_id <- clean_sample_ids(pd_model$Anonymised_sampleID)
pd_model <- pd_model[match(master_ids, pd_model$sample_id), , drop = FALSE]
if (!identical(pd_model$sample_id, master_ids)) {
  stop("Failed to align pd_model to master sample order.", call. = FALSE)
}

# Load metabolomics scale matrix from the processed object.
metab_processed_file <- file.path(paths$analysis_dataset, "metabolomics_processed.rds")
assert_file_exists(metab_processed_file, "metabolomics_processed")
metab_obj <- readRDS(metab_processed_file)
if (!all(c("metab_scaled", "sample_ids") %in% names(metab_obj))) {
  stop("metabolomics_processed.rds must contain metab_scaled and sample_ids.", call. = FALSE)
}
metab_scaled <- as.matrix(metab_obj$metab_scaled)
metab_ids <- as.character(metab_obj$sample_ids)
metab_scaled <- metab_scaled[match(master_ids, metab_ids), , drop = FALSE]
rownames(metab_scaled) <- master_ids
if (!identical(rownames(metab_scaled), master_ids)) {
  stop("Failed to align metab_scaled to master sample order.", call. = FALSE)
}

required_clinical <- c("AGE_num", "UPDRS_III_num", "MOCA_total_num", "LEDD_total_num", "GENDER")
missing_clinical <- setdiff(required_clinical, names(pd_model))
if (length(missing_clinical) > 0) {
  stop("pd_model is missing required columns: ", paste(missing_clinical, collapse = ", "), call. = FALSE)
}

# -----------------------------------------------------------------------------
# 3. Branch definitions and output folders
# -----------------------------------------------------------------------------
branch_names <- c(
  "clinical_only",
  "no_LEDD",
  "representation_comparison",
  "top20_combined",
  "top20_23network"
)

out_root <- get_posthoc_output_dir("switcher_characterisation")
plot_dir <- file.path(out_root, "plots")
table_dir <- file.path(out_root, "tables")
diag_dir <- file.path(out_root, "diagnostics")
safe_dir(plot_dir)
safe_dir(table_dir)
safe_dir(diag_dir)

load_branch_switchers <- function(branch_name) {
  branch_dir <- file.path(paths$posthoc, "cluster_switchers")
  candidates <- c(
    file.path(branch_dir, "tables", paste0(branch_name, "_switchers.rds")),
    file.path(branch_dir, paste0(branch_name, "_switchers.rds"))
  )
  f <- find_first_existing(candidates)
  if (is.na(f)) {
    stop(
      "Could not find switcher table for branch '", branch_name, "'. Looked in:\n",
      paste(candidates, collapse = "\n"),
      call. = FALSE
    )
  }
  obj <- readRDS(f)
  needed <- c("sample_id", "main_cluster", "branch_cluster", "aligned_cluster", "switched")
  missing <- setdiff(needed, names(obj))
  if (length(missing) > 0) {
    stop(
      "Switcher table for branch '", branch_name, "' is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  obj$sample_id <- clean_sample_ids(obj$sample_id)
  obj$main_cluster <- parse_cluster_id(obj$main_cluster)
  obj$branch_cluster <- parse_cluster_id(obj$branch_cluster)
  obj$aligned_cluster <- parse_cluster_id(obj$aligned_cluster)
  obj$switched <- as.logical(obj$switched)
  obj
}

# Optional borderline proxy from the switcher script, if present.
load_margin_proxy <- function() {
  f <- file.path(paths$posthoc, "cluster_switchers", "combined_switcher_cluster_distance_proxy.rds")
  if (file.exists(f)) {
    readRDS(f)
  } else {
    NULL
  }
}

# -----------------------------------------------------------------------------
# 4. Per-branch characterisation
# -----------------------------------------------------------------------------
clinical_summary_rows <- list()
clinical_test_rows <- list()
metabolite_test_rows <- list()
branch_overviews <- list()

margin_proxy <- load_margin_proxy()
margin_proxy_rows <- list()

for (branch_name in branch_names) {
  message("Characterising switchers for branch: ", branch_name)

  switch_df <- load_branch_switchers(branch_name)
  switch_flag <- data.frame(
    sample_id = switch_df$sample_id,
    switched = switch_df$switched,
    stringsAsFactors = FALSE
  )

  branch_data <- merge(
    pd_model[, c("sample_id", required_clinical), drop = FALSE],
    switch_flag,
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )
  branch_data <- branch_data[match(master_ids, branch_data$sample_id), , drop = FALSE]
  if (!identical(branch_data$sample_id, master_ids)) {
    stop("Failed to align branch clinical data for branch: ", branch_name, call. = FALSE)
  }
  if (anyNA(branch_data$switched)) {
    stop("Missing switch flags for branch: ", branch_name, call. = FALSE)
  }

  branch_data$switch_group <- factor(ifelse(branch_data$switched, "Switcher", "Non-switcher"),
                                     levels = c("Non-switcher", "Switcher"))
  branch_data$GENDER <- factor(branch_data$GENDER)

  # Clinical summary.
  cl_sum <- data.frame(
    branch = branch_name,
    n_total = nrow(branch_data),
    n_switchers = sum(branch_data$switched),
    prop_switchers = mean(branch_data$switched),
    switcher_age_median = median(branch_data$AGE_num[branch_data$switched], na.rm = TRUE),
    nonswitcher_age_median = median(branch_data$AGE_num[!branch_data$switched], na.rm = TRUE),
    switcher_updrs_median = median(branch_data$UPDRS_III_num[branch_data$switched], na.rm = TRUE),
    nonswitcher_updrs_median = median(branch_data$UPDRS_III_num[!branch_data$switched], na.rm = TRUE),
    switcher_moca_median = median(branch_data$MOCA_total_num[branch_data$switched], na.rm = TRUE),
    nonswitcher_moca_median = median(branch_data$MOCA_total_num[!branch_data$switched], na.rm = TRUE),
    switcher_ledd_median = median(branch_data$LEDD_total_num[branch_data$switched], na.rm = TRUE),
    nonswitcher_ledd_median = median(branch_data$LEDD_total_num[!branch_data$switched], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  clinical_summary_rows[[branch_name]] <- cl_sum
  write.csv(cl_sum, file.path(table_dir, paste0(branch_name, "_clinical_summary.csv")), row.names = FALSE)
  saveRDS(cl_sum, file.path(table_dir, paste0(branch_name, "_clinical_summary.rds")))

  # Clinical tests.
  clin_vars <- c("AGE_num", "UPDRS_III_num", "MOCA_total_num", "LEDD_total_num")
  ct <- list()
  for (v in clin_vars) {
    x <- branch_data[[v]]
    g <- branch_data$switch_group
    if (length(unique(g)) < 2) next
    wt <- suppressWarnings(wilcox.test(x ~ g, exact = FALSE))
    tt <- suppressWarnings(t.test(x ~ g))
    ct[[v]] <- data.frame(
      branch = branch_name,
      variable = v,
      switcher_median = median(x[g == "Switcher"], na.rm = TRUE),
      nonswitcher_median = median(x[g == "Non-switcher"], na.rm = TRUE),
      wilcox_p = wt$p.value,
      t_p = tt$p.value,
      cohens_d = cohens_d(x, g),
      stringsAsFactors = FALSE
    )
  }
  if (length(ct) > 0) {
    ct <- do.call(rbind, ct)
    ct$wilcox_p_adj <- p.adjust(ct$wilcox_p, method = "BH")
    ct$t_p_adj <- p.adjust(ct$t_p, method = "BH")
    clinical_test_rows[[branch_name]] <- ct
    write.csv(ct, file.path(table_dir, paste0(branch_name, "_clinical_tests.csv")), row.names = FALSE)
    saveRDS(ct, file.path(table_dir, paste0(branch_name, "_clinical_tests.rds")))
  }

  # Gender table.
  gender_tab <- table(branch_data$switch_group, branch_data$GENDER)
  gender_df <- as.data.frame.matrix(gender_tab)
  gender_df$switch_group <- rownames(gender_df)
  rownames(gender_df) <- NULL
  write.csv(gender_df, file.path(table_dir, paste0(branch_name, "_gender_table.csv")), row.names = FALSE)
  saveRDS(gender_df, file.path(table_dir, paste0(branch_name, "_gender_table.rds")))

  # Switchers by cluster.
  cluster_tab <- aggregate(
    switched ~ aligned_cluster,
    data = switch_df,
    FUN = function(z) c(n = length(z), switched = sum(z), prop = mean(z))
  )
  cluster_tab <- do.call(data.frame, cluster_tab)
  names(cluster_tab) <- c("aligned_cluster", "n_samples", "n_switchers", "prop_switchers")
  cluster_tab$branch <- branch_name
  write.csv(cluster_tab, file.path(table_dir, paste0(branch_name, "_switchers_by_cluster.csv")), row.names = FALSE)
  saveRDS(cluster_tab, file.path(table_dir, paste0(branch_name, "_switchers_by_cluster.rds")))

  # Borderline proxy if available: cluster-margin summaries.
  if (!is.null(margin_proxy)) {
    mp <- margin_proxy[margin_proxy$branch == branch_name, , drop = FALSE]
    if (nrow(mp) > 0) {
      mp$branch_name <- branch_name
      margin_proxy_rows[[branch_name]] <- mp
      write.csv(mp, file.path(table_dir, paste0(branch_name, "_margin_proxy.csv")), row.names = FALSE)
      saveRDS(mp, file.path(table_dir, paste0(branch_name, "_margin_proxy.rds")))
    }
  }

  # Metabolomics comparisons: switcher vs non-switcher.
  metab_df <- data.frame(
    sample_id = master_ids,
    switch_group = branch_data$switch_group,
    stringsAsFactors = FALSE
  )
  metab_df <- cbind(metab_df, as.data.frame(metab_scaled, check.names = FALSE))

  metabolite_names <- setdiff(names(metab_df), c("sample_id", "switch_group"))
  mt <- vector("list", length(metabolite_names))
  names(mt) <- metabolite_names

  for (m in metabolite_names) {
    x <- metab_df[[m]]
    g <- metab_df$switch_group
    if (!is.numeric(x)) next
    wt <- suppressWarnings(wilcox.test(x ~ g, exact = FALSE))
    tt <- suppressWarnings(t.test(x ~ g))
    mt[[m]] <- data.frame(
      branch = branch_name,
      metabolite = m,
      switcher_mean = mean(x[g == "Switcher"], na.rm = TRUE),
      nonswitcher_mean = mean(x[g == "Non-switcher"], na.rm = TRUE),
      switcher_median = median(x[g == "Switcher"], na.rm = TRUE),
      nonswitcher_median = median(x[g == "Non-switcher"], na.rm = TRUE),
      mean_diff = mean(x[g == "Switcher"], na.rm = TRUE) - mean(x[g == "Non-switcher"], na.rm = TRUE),
      wilcox_p = wt$p.value,
      t_p = tt$p.value,
      cohens_d = cohens_d(x, g),
      stringsAsFactors = FALSE
    )
  }

  mt <- mt[!vapply(mt, is.null, logical(1))]
  if (length(mt) > 0) {
    mt <- do.call(rbind, mt)
    mt$wilcox_p_adj <- p.adjust(mt$wilcox_p, method = "BH")
    mt$t_p_adj <- p.adjust(mt$t_p, method = "BH")
    mt$neglog10_p <- -log10(pmax(mt$wilcox_p_adj, .Machine$double.xmin))
    mt$sig <- mt$wilcox_p_adj < 0.05
    metabolite_test_rows[[branch_name]] <- mt

    write.csv(mt, file.path(table_dir, paste0(branch_name, "_metabolite_tests.csv")), row.names = FALSE)
    saveRDS(mt, file.path(table_dir, paste0(branch_name, "_metabolite_tests.rds")))

    # Volcano plot.
    plot_volcano(
      mt,
      file.path(plot_dir, paste0(branch_name, "_metabolite_volcano.png")),
      title = paste0(branch_name, ": metabolite switcher volcano plot")
    )

    # Top hits.
    top_hits <- mt[order(mt$wilcox_p_adj), , drop = FALSE]
    top_hits <- top_hits[seq_len(min(20, nrow(top_hits))), , drop = FALSE]
    write.csv(top_hits, file.path(table_dir, paste0(branch_name, "_top_metabolite_hits.csv")), row.names = FALSE)
    saveRDS(top_hits, file.path(table_dir, paste0(branch_name, "_top_metabolite_hits.rds")))

    # Heatmap of top metabolites.
    top_n <- min(30, nrow(mt))
    top_metabs <- mt$metabolite[order(mt$wilcox_p_adj)][seq_len(top_n)]
    plot_top_metab_heatmap(
      metab_df[, c("sample_id", "switch_group", top_metabs), drop = FALSE],
      top_metabs,
      file.path(plot_dir, paste0(branch_name, "_top_metabolite_heatmap.png")),
      title = paste0(branch_name, ": top metabolite pattern")
    )

    if (has_ggrepel) {
      sig_hits <- mt[order(mt$wilcox_p_adj), , drop = FALSE]
      sig_hits <- sig_hits[sig_hits$sig, , drop = FALSE]
      sig_hits <- sig_hits[seq_len(min(10, nrow(sig_hits))), , drop = FALSE]
      if (nrow(sig_hits) > 0) {
        p <- ggplot2::ggplot(mt, ggplot2::aes(x = mean_diff, y = neglog10_p)) +
          ggplot2::geom_point(ggplot2::aes(color = sig), alpha = 0.7, size = 1.5) +
          ggrepel::geom_text_repel(
            data = sig_hits,
            ggplot2::aes(label = metabolite),
            size = 3,
            max.overlaps = Inf
          ) +
          ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
          ggplot2::theme_minimal(base_size = 12) +
          ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "bottom") +
          ggplot2::labs(
            x = "Mean difference (switcher - non-switcher)",
            y = "-log10 BH-adjusted p",
            color = NULL,
            title = paste0(branch_name, ": labelled metabolite volcano plot")
          )
        ggplot2::ggsave(file.path(plot_dir, paste0(branch_name, "_metabolite_volcano_labelled.png")),
                        p, width = 7.2, height = 5.8, dpi = 300)
      }
    }
  }

  # Clinical boxplots.
  for (v in clin_vars) {
    df <- branch_data[, c("switch_group", v), drop = FALSE]
    names(df)[2] <- "value"
    plot_boxplot(
      df,
      x = "switch_group",
      y = "value",
      out_file = file.path(plot_dir, paste0(branch_name, "_", v, "_switcher_boxplot.png")),
      title = paste0(branch_name, ": ", v, " in switchers vs non-switchers"),
      ylab = v
    )
  }

  branch_overviews[[branch_name]] <- data.frame(
    branch = branch_name,
    n_total = nrow(branch_data),
    n_switchers = sum(branch_data$switched),
    prop_switchers = mean(branch_data$switched),
    stringsAsFactors = FALSE
  )
  write.csv(branch_overviews[[branch_name]],
            file.path(diag_dir, paste0(branch_name, "_switcher_overview.csv")),
            row.names = FALSE)
  saveRDS(branch_overviews[[branch_name]],
          file.path(diag_dir, paste0(branch_name, "_switcher_overview.rds")))
}

# -----------------------------------------------------------------------------
# 5. Combined outputs
# -----------------------------------------------------------------------------
combined_clinical_summary <- do.call(rbind, clinical_summary_rows)
rownames(combined_clinical_summary) <- NULL
write.csv(combined_clinical_summary, file.path(out_root, "combined_switcher_clinical_summary.csv"), row.names = FALSE)
saveRDS(combined_clinical_summary, file.path(out_root, "combined_switcher_clinical_summary.rds"))

if (length(clinical_test_rows) > 0) {
  combined_clinical_tests <- do.call(rbind, clinical_test_rows)
  rownames(combined_clinical_tests) <- NULL
  write.csv(combined_clinical_tests, file.path(out_root, "combined_switcher_clinical_tests.csv"), row.names = FALSE)
  saveRDS(combined_clinical_tests, file.path(out_root, "combined_switcher_clinical_tests.rds"))
}

if (length(metabolite_test_rows) > 0) {
  combined_metab_tests <- do.call(rbind, metabolite_test_rows)
  rownames(combined_metab_tests) <- NULL
  write.csv(combined_metab_tests, file.path(out_root, "combined_switcher_metabolite_tests.csv"), row.names = FALSE)
  saveRDS(combined_metab_tests, file.path(out_root, "combined_switcher_metabolite_tests.rds"))
}

if (length(branch_overviews) > 0) {
  overview <- do.call(rbind, branch_overviews)
  rownames(overview) <- NULL
  write.csv(overview, file.path(out_root, "switcher_overview.csv"), row.names = FALSE)
  saveRDS(overview, file.path(out_root, "switcher_overview.rds"))
}

if (length(margin_proxy_rows) > 0) {
  combined_margin_proxy <- do.call(rbind, margin_proxy_rows)
  rownames(combined_margin_proxy) <- NULL
  write.csv(combined_margin_proxy, file.path(out_root, "combined_switcher_margin_proxy.csv"), row.names = FALSE)
  saveRDS(combined_margin_proxy, file.path(out_root, "combined_switcher_margin_proxy.rds"))

  p <- ggplot2::ggplot(
    combined_margin_proxy,
    ggplot2::aes(x = factor(switched, levels = c(FALSE, TRUE), labels = c("Non-switcher", "Switcher")),
                 y = margin, fill = factor(switched))
  ) +
    ggplot2::geom_boxplot(outlier.alpha = 0.2, width = 0.6) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "none") +
    ggplot2::labs(x = NULL, y = "Best minus second-best fused-network similarity",
                  title = "Borderline proxy: cluster separation margin")
  ggplot2::ggsave(file.path(plot_dir, "combined_cluster_margin_boxplot.png"), p, width = 6.4, height = 4.6, dpi = 300)
}

cat("Switch characterisation complete.\n")
cat("Output directory: ", out_root, "\n", sep = "")
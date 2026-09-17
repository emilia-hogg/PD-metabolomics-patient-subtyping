# 06_posthoc_analysis/05_cluster_switchers.R
# Re-runnable from a fresh R session.
#
# - compare the primary integrated SNF cluster solution to each sensitivity branch
# - identify patients whose assignments change after optimal label permutation
# - summarise switch frequencies, overlaps, and agreement metrics


# -----------------------------------------------------------------------------
# 0. Load config, packages, and shared helpers
# -----------------------------------------------------------------------------
project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))
source(file.path(project_root, "06_posthoc_analysis", "00_shared_characterisation_helpers.R"))

required_pkgs <- c("ggplot2")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required package: ", pkg, call. = FALSE)
  }
}

use_ggalluvial <- requireNamespace("ggalluvial", quietly = TRUE)


# -----------------------------------------------------------------------------
# 1. Small utilities
# -----------------------------------------------------------------------------
safe_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

clean_sample_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[nchar(x) == 0] <- NA_character_
  x
}

parse_cluster_id <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.numeric(x) || is.integer(x)) {
    out <- as.integer(x)
    if (anyNA(out)) {
      stop("Cluster labels contain NA after numeric conversion.", call. = FALSE)
    }
    return(out)
  }

  x <- as.character(x)
  if (length(x) == 0) return(integer(0))

  extracted <- suppressWarnings(as.integer(gsub("^.*?(-?[0-9]+).*$", "\\1", x)))
  if (!anyNA(extracted)) {
    return(as.integer(extracted))
  }

  levs <- unique(x)
  map <- setNames(seq_along(levs), levs)
  out <- unname(map[x])
  if (anyNA(out)) {
    stop("Some cluster labels could not be converted to integers.", call. = FALSE)
  }
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
  if (length(hit) == 0L || is.na(hit) || !nzchar(hit)) {
    return(NA_character_)
  }
  hit
}

align_branch_to_primary <- function(primary_df, branch_df, branch_name) {
  if (!all(primary_df$sample_id %in% branch_df$sample_id)) {
    missing_ids <- setdiff(primary_df$sample_id, branch_df$sample_id)
    stop(
      "Branch '", branch_name, "' is missing ", length(missing_ids),
      " primary sample IDs. First missing ID: ", missing_ids[1],
      call. = FALSE
    )
  }

  extra_ids <- setdiff(branch_df$sample_id, primary_df$sample_id)
  if (length(extra_ids) > 0) {
    stop(
      "Branch '", branch_name, "' contains sample IDs not present in the primary cohort. First extra ID: ",
      extra_ids[1],
      call. = FALSE
    )
  }

  branch_df <- branch_df[match(primary_df$sample_id, branch_df$sample_id), , drop = FALSE]
  if (!identical(branch_df$sample_id, primary_df$sample_id)) {
    stop("Failed to align branch '", branch_name, "' to primary sample order.", call. = FALSE)
  }
  branch_df
}

safe_prop <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

make_switch_summary <- function(df, branch_name) {
  data.frame(
    branch = branch_name,
    n_samples = nrow(df),
    n_switchers = sum(df$switched, na.rm = TRUE),
    prop_switchers = safe_prop(df$switched),
    n_primary_clusters = length(unique(df$main_cluster)),
    n_branch_clusters = length(unique(df$branch_cluster)),
    stringsAsFactors = FALSE
  )
}

make_overlap_table <- function(main_labels, aligned_branch_labels, branch_name) {
  tab <- table(main_labels, aligned_branch_labels)
  out <- as.data.frame.matrix(tab)
  out$primary_cluster <- rownames(out)
  rownames(out) <- NULL
  out <- out[, c("primary_cluster", setdiff(names(out), "primary_cluster")), drop = FALSE]
  attr(out, "branch_name") <- branch_name
  out
}

make_jaccard_matrix <- function(branch_sets) {
  branch_names <- names(branch_sets)
  mat <- matrix(
    NA_real_,
    nrow = length(branch_names),
    ncol = length(branch_names),
    dimnames = list(branch_names, branch_names)
  )
  for (i in seq_along(branch_names)) {
    for (j in seq_along(branch_names)) {
      a <- branch_sets[[branch_names[i]]]
      b <- branch_sets[[branch_names[j]]]
      uni <- length(union(a, b))
      mat[i, j] <- if (uni == 0) NA_real_ else length(intersect(a, b)) / uni
    }
  }
  mat
}

plot_switcher_bar <- function(summary_df, out_file, title = "Switcher proportion by sensitivity branch") {
  p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = branch, y = prop_switchers, fill = branch)) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
    ggplot2::labs(x = NULL, y = "Patients switching cluster", title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(out_file, p, width = 7.2, height = 4.6, dpi = 300)
  invisible(p)
}

plot_switcher_counts <- function(summary_df, out_file, title = "Switcher counts by branch") {
  plot_df <- data.frame(
    branch = rep(summary_df$branch, times = 2),
    status = rep(c("Stayed", "Switched"), each = nrow(summary_df)),
    count = c(summary_df$n_samples - summary_df$n_switchers, summary_df$n_switchers),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = branch, y = count, fill = status)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::labs(x = NULL, y = "Number of patients", fill = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(out_file, p, width = 7.2, height = 4.6, dpi = 300)
  invisible(p)
}

plot_heatmap <- function(mat, out_file, title, fill_title = "Count") {
  # If mat is a data.frame (e.g. overlap table with a row-label column), convert
  # to a pure numeric matrix first so scale_fill_gradient gets continuous values.
  if (is.data.frame(mat)) {
    row_label_col <- which(sapply(mat, function(x) !is.numeric(x)))
    if (length(row_label_col) > 0) {
      row_labels <- mat[[row_label_col[1]]]
      num_mat <- as.matrix(mat[, -row_label_col, drop = FALSE])
      mode(num_mat) <- "numeric"
      rownames(num_mat) <- row_labels
    } else {
      num_mat <- as.matrix(mat)
      mode(num_mat) <- "numeric"
    }
    mat <- num_mat
  } else {
    mat <- as.matrix(mat)
    mode(mat) <- "numeric"
  }

  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(df) <- c("row_var", "col_var", "value")
  df$row_var <- factor(df$row_var, levels = unique(df$row_var))
  df$col_var <- factor(df$col_var, levels = unique(df$col_var))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = col_var, y = row_var, fill = value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(value), "", value)), size = 3) +
    ggplot2::scale_fill_gradient(low = "grey95", high = "steelblue", na.value = "grey95") +
    ggplot2::labs(x = NULL, y = NULL, fill = fill_title, title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  ggplot2::ggsave(out_file, p, width = 7.8, height = 6.0, dpi = 300)
  invisible(p)
}

plot_alluvial_switches <- function(df, branch_name, out_file) {
  df <- df[, c("sample_id", "main_cluster", "aligned_cluster", "switched"), drop = FALSE]
  df$main_cluster <- factor(df$main_cluster, levels = sort(unique(c(df$main_cluster, df$aligned_cluster))))
  df$aligned_cluster <- factor(df$aligned_cluster, levels = levels(df$main_cluster))
  df$switched <- factor(ifelse(df$switched, "Switched", "Stayed"), levels = c("Stayed", "Switched"))

  if (use_ggalluvial) {
    p <- ggplot2::ggplot(df, ggplot2::aes(axis1 = main_cluster, axis2 = aligned_cluster, y = 1)) +
      ggalluvial::geom_alluvium(ggplot2::aes(fill = switched), width = 0.18, alpha = 0.8) +
      ggalluvial::geom_stratum(width = 0.18, color = "grey35", fill = "grey92") +
      ggalluvial::stat_stratum(geom = "text", ggplot2::aes(label = after_stat(stratum)), size = 3) +
      ggplot2::scale_x_discrete(limits = c("Primary", branch_name), expand = c(0.08, 0.08)) +
      ggplot2::labs(
        x = NULL,
        y = "Patients",
        fill = NULL,
        title = paste0("Primary vs ", branch_name, " cluster flow")
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  } else {
    plot_df <- data.frame(
      main_cluster = df$main_cluster,
      aligned_cluster = df$aligned_cluster,
      switched = df$switched,
      stringsAsFactors = FALSE
    )
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = main_cluster, fill = switched)) +
      ggplot2::geom_bar(position = "fill") +
      ggplot2::facet_wrap(~ aligned_cluster) +
      ggplot2::labs(
        x = "Primary cluster",
        y = "Proportion",
        title = paste0("Primary vs ", branch_name, " cluster flow")
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }

  ggplot2::ggsave(out_file, p, width = 8.5, height = 5.8, dpi = 300)
  invisible(p)
}

plot_alluvial_panel <- function(combined_df, out_file) {
  if (!use_ggalluvial) {
    return(invisible(NULL))
  }

  df <- combined_df[, c("branch", "main_cluster", "aligned_cluster", "switched"), drop = FALSE]
  df$branch <- factor(df$branch, levels = unique(df$branch))
  df$main_cluster <- factor(df$main_cluster, levels = sort(unique(c(df$main_cluster, df$aligned_cluster))))
  df$aligned_cluster <- factor(df$aligned_cluster, levels = levels(df$main_cluster))
  df$switched <- factor(ifelse(df$switched, "Switched", "Stayed"), levels = c("Stayed", "Switched"))

  p <- ggplot2::ggplot(df, ggplot2::aes(axis1 = main_cluster, axis2 = aligned_cluster, y = 1)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = switched), width = 0.18, alpha = 0.75) +
    ggalluvial::geom_stratum(width = 0.18, color = "grey35", fill = "grey92") +
    ggalluvial::stat_stratum(geom = "text", ggplot2::aes(label = after_stat(stratum)), size = 2.7) +
    ggplot2::facet_wrap(~ branch, nrow = 1) +
    ggplot2::labs(
      x = NULL,
      y = "Patients",
      fill = NULL,
      title = "Primary vs sensitivity-branch cluster flow"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(out_file, p, width = 14, height = 5.8, dpi = 300)
  invisible(p)
}

load_branch_result <- function(branch_name, result_candidates) {
  result_file <- find_first_existing(result_candidates)
  if (is.na(result_file)) {
    stop(
      "Could not find sensitivity result file for branch '", branch_name, "'.\nLooked in:\n",
      paste(result_candidates, collapse = "\n"),
      call. = FALSE
    )
  }

  obj <- readRDS(result_file)

  sample_ids <- NULL
  cluster_obj <- NULL

  if (is.list(obj) && !is.null(obj$cluster_labels) && !is.null(obj$sample_ids)) {
    cluster_obj <- obj$cluster_labels
    sample_ids <- clean_sample_ids(obj$sample_ids)
  } else if (is.list(obj) && !is.null(obj$best_k_cluster_labels)) {
    cluster_obj <- obj$best_k_cluster_labels
    if (!is.null(obj$sample_ids)) {
      sample_ids <- clean_sample_ids(obj$sample_ids)
    }
  } else if (branch_name == "representation_comparison" && is.list(obj)) {
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
        "representation_comparison result was found, but no cluster labels field was recognisable.",
        call. = FALSE
      )
    }
    if (!is.null(obj$sample_ids)) {
      sample_ids <- clean_sample_ids(obj$sample_ids)
    }
  } else if (is.data.frame(obj) || (is.vector(obj) && !is.null(names(obj)))) {
    cluster_obj <- obj
  } else {
    stop(
      "Unrecognised cluster result structure for branch '", branch_name, "' in file: ", result_file,
      call. = FALSE
    )
  }

  cluster_df <- standardise_cluster_df(
    obj = cluster_obj,
    fallback_sample_ids = sample_ids,
    label = paste0(branch_name, " cluster labels")
  )

  list(
    branch_name = branch_name,
    result_file = result_file,
    cluster_df = cluster_df
  )
}


# -----------------------------------------------------------------------------
# 2. Load cohort and primary labels
# -----------------------------------------------------------------------------
objs <- load_analysis_objects()
master_ids <- as.character(objs$sample_ids)

primary <- load_primary_snf_outputs()
primary_df <- standardise_cluster_df(primary$cluster_labels, label = "primary cluster labels")

primary_df <- primary_df[match(primary$sample_ids, primary_df$sample_id), , drop = FALSE]
if (!identical(primary_df$sample_id, primary$sample_ids)) {
  stop("Failed to align primary labels to primary SNF sample order.", call. = FALSE)
}

primary_df <- primary_df[match(master_ids, primary_df$sample_id), , drop = FALSE]
if (!identical(primary_df$sample_id, master_ids)) {
  stop("Failed to align primary labels to master sample order.", call. = FALSE)
}

if (anyDuplicated(primary_df$sample_id)) {
  stop("Primary labels contain duplicated sample IDs after alignment.", call. = FALSE)
}


# -----------------------------------------------------------------------------
# 3. Branch definitions
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
# 4. Output directories
# -----------------------------------------------------------------------------
out_root <- get_posthoc_output_dir("cluster_switchers")
plot_dir <- file.path(out_root, "plots")
table_dir <- file.path(out_root, "tables")
diag_dir <- file.path(out_root, "diagnostics")
safe_dir(plot_dir)
safe_dir(table_dir)
safe_dir(diag_dir)


# -----------------------------------------------------------------------------
# 5. Run comparisons branch by branch
# -----------------------------------------------------------------------------
branch_results <- list()
comparison_rows <- list()
overlap_tables <- list()
switch_tables <- list()
switch_sets <- list()
primary_cluster_profiles <- list()
branch_cluster_profiles <- list()

for (branch_name in names(branch_specs)) {
  message("Comparing primary clusters to branch: ", branch_name)

  branch <- load_branch_result(
    branch_name = branch_name,
    result_candidates = branch_specs[[branch_name]]$result_candidates
  )

  branch_df <- align_branch_to_primary(primary_df, branch$cluster_df, branch_name)

  comparison <- compare_cluster_solutions(primary_df$cluster, branch_df$cluster)
  alignment <- best_label_permutation(primary_df$cluster, branch_df$cluster)

  switch_df <- data.frame(
    sample_id = primary_df$sample_id,
    main_cluster = primary_df$cluster,
    branch_cluster = branch_df$cluster,
    aligned_cluster = alignment$aligned_labels,
    switched = primary_df$cluster != alignment$aligned_labels,
    stringsAsFactors = FALSE
  )
  switch_df$branch_name <- branch_name

  overlap <- make_overlap_table(
    main_labels = switch_df$main_cluster,
    aligned_branch_labels = switch_df$aligned_cluster,
    branch_name = branch_name
  )

  summary_df <- make_switch_summary(switch_df, branch_name)
  primary_profile <- summarise_switchers_by_cluster(switch_df, "main_cluster", "primary_cluster")
  branch_profile <- summarise_switchers_by_cluster(switch_df, "aligned_cluster", "branch_cluster")

  comparison_row <- data.frame(
    branch = branch_name,
    n_samples = nrow(switch_df),
    n_switchers = sum(switch_df$switched, na.rm = TRUE),
    prop_switchers = mean(switch_df$switched, na.rm = TRUE),
    n_primary_clusters = length(unique(switch_df$main_cluster)),
    n_branch_clusters = length(unique(switch_df$branch_cluster)),
    ari = comparison$ari,
    nmi = comparison$nmi,
    vi = comparison$vi,
    exact_direct = comparison$exact_direct,
    exact_best = comparison$exact_best,
    moved_n = comparison$moved_n,
    moved_prop = comparison$moved_prop,
    stringsAsFactors = FALSE
  )

  comparison_rows[[branch_name]] <- comparison_row
  overlap_tables[[branch_name]] <- overlap
  switch_tables[[branch_name]] <- switch_df
  switch_sets[[branch_name]] <- switch_df$sample_id[switch_df$switched]
  primary_cluster_profiles[[branch_name]] <- primary_profile
  branch_cluster_profiles[[branch_name]] <- branch_profile

  branch_results[[branch_name]] <- list(
    branch_name = branch_name,
    result_file = branch$result_file,
    comparison = comparison,
    switch_df = switch_df,
    overlap = overlap,
    summary = summary_df,
    primary_profile = primary_profile,
    branch_profile = branch_profile
  )

  saveRDS(switch_df, file.path(table_dir, paste0(branch_name, "_switchers.rds")))
  write.csv(switch_df, file.path(table_dir, paste0(branch_name, "_switchers.csv")), row.names = FALSE)

  saveRDS(overlap, file.path(table_dir, paste0(branch_name, "_overlap_table.rds")))
  write.csv(overlap, file.path(table_dir, paste0(branch_name, "_overlap_table.csv")), row.names = FALSE)

  saveRDS(summary_df, file.path(table_dir, paste0(branch_name, "_switch_summary.rds")))
  write.csv(summary_df, file.path(table_dir, paste0(branch_name, "_switch_summary.csv")), row.names = FALSE)

  saveRDS(primary_profile, file.path(table_dir, paste0(branch_name, "_primary_cluster_switch_profile.rds")))
  write.csv(primary_profile, file.path(table_dir, paste0(branch_name, "_primary_cluster_switch_profile.csv")), row.names = FALSE)

  saveRDS(branch_profile, file.path(table_dir, paste0(branch_name, "_branch_cluster_switch_profile.rds")))
  write.csv(branch_profile, file.path(table_dir, paste0(branch_name, "_branch_cluster_switch_profile.csv")), row.names = FALSE)

  write.csv(comparison_row, file.path(diag_dir, paste0(branch_name, "_comparison_metrics.csv")), row.names = FALSE)
  saveRDS(comparison_row, file.path(diag_dir, paste0(branch_name, "_comparison_metrics.rds")))

  plot_switcher_bar(
    summary_df,
    file.path(plot_dir, paste0(branch_name, "_switcher_proportion.png")),
    title = paste0(branch_name, ": switcher proportion")
  )

  plot_switcher_counts(
    summary_df,
    file.path(plot_dir, paste0(branch_name, "_switcher_counts.png")),
    title = paste0(branch_name, ": switcher counts")
  )

  plot_alluvial_switches(
    switch_df,
    branch_name,
    file.path(plot_dir, paste0(branch_name, "_alluvial.png"))
  )

  plot_heatmap(
    overlap,
    file.path(plot_dir, paste0(branch_name, "_overlap_heatmap.png")),
    title = paste0("Primary x ", branch_name, " overlap"),
    fill_title = "Patients"
  )
}


# -----------------------------------------------------------------------------
# 6. Combined outputs
# -----------------------------------------------------------------------------
comparison_df <- do.call(rbind, comparison_rows)
rownames(comparison_df) <- NULL
write.csv(comparison_df, file.path(out_root, "primary_vs_sensitivity_comparison_metrics.csv"), row.names = FALSE)
saveRDS(comparison_df, file.path(out_root, "primary_vs_sensitivity_comparison_metrics.rds"))

switch_summary_df <- do.call(rbind, lapply(names(switch_tables), function(branch_name) {
  make_switch_summary(switch_tables[[branch_name]], branch_name)
}))
rownames(switch_summary_df) <- NULL
write.csv(switch_summary_df, file.path(out_root, "switch_summary_by_branch.csv"), row.names = FALSE)
saveRDS(switch_summary_df, file.path(out_root, "switch_summary_by_branch.rds"))

combined_switchers <- do.call(rbind, switch_tables)
rownames(combined_switchers) <- NULL
write.csv(combined_switchers, file.path(out_root, "all_branch_switchers_long.csv"), row.names = FALSE)
saveRDS(combined_switchers, file.path(out_root, "all_branch_switchers_long.rds"))

combined_profiles <- list(
  primary_cluster_profiles = primary_cluster_profiles,
  branch_cluster_profiles = branch_cluster_profiles
)
saveRDS(combined_profiles, file.path(out_root, "switcher_cluster_profiles.rds"))

combined_objects <- list(
  branch_overlaps = overlap_tables,
  branch_results = branch_results,
  branch_switch_sets = switch_sets
)
saveRDS(combined_objects, file.path(out_root, "branch_overlap_objects.rds"))

branch_names <- names(switch_sets)
switch_count_mat <- matrix(
  0,
  nrow = length(branch_names),
  ncol = length(branch_names),
  dimnames = list(branch_names, branch_names)
)
for (i in seq_along(branch_names)) {
  for (j in seq_along(branch_names)) {
    a <- switch_sets[[branch_names[i]]]
    b <- switch_sets[[branch_names[j]]]
    switch_count_mat[i, j] <- length(intersect(a, b))
  }
}
switch_jaccard_mat <- make_jaccard_matrix(switch_sets)

write.csv(switch_count_mat, file.path(out_root, "switcher_overlap_counts.csv"))
write.csv(switch_jaccard_mat, file.path(out_root, "switcher_overlap_jaccard.csv"))
saveRDS(switch_count_mat, file.path(out_root, "switcher_overlap_counts.rds"))
saveRDS(switch_jaccard_mat, file.path(out_root, "switcher_overlap_jaccard.rds"))

plot_heatmap(
  switch_count_mat,
  file.path(plot_dir, "switcher_overlap_counts_heatmap.png"),
  title = "Overlap of switchers across branches",
  fill_title = "Shared switchers"
)

plot_heatmap(
  switch_jaccard_mat,
  file.path(plot_dir, "switcher_overlap_jaccard_heatmap.png"),
  title = "Jaccard overlap of switchers across branches",
  fill_title = "Jaccard"
)

combined_alluvial_df <- do.call(rbind, lapply(names(switch_tables), function(branch_name) {
  df <- switch_tables[[branch_name]][, c("sample_id", "main_cluster", "aligned_cluster", "switched"), drop = FALSE]
  df$branch <- branch_name
  df
}))
rownames(combined_alluvial_df) <- NULL
saveRDS(combined_alluvial_df, file.path(out_root, "combined_alluvial_df.rds"))
write.csv(combined_alluvial_df, file.path(out_root, "combined_alluvial_df.csv"), row.names = FALSE)
plot_alluvial_panel(combined_alluvial_df, file.path(plot_dir, "all_branches_alluvial_panel.png"))


# -----------------------------------------------------------------------------
# 7. Final checks and messages
# -----------------------------------------------------------------------------
cat("Primary sample count: ", nrow(primary_df), "\n", sep = "")
cat("Branches analysed: ", paste(names(branch_specs), collapse = ", "), "\n", sep = "")
cat("Output directory: ", out_root, "\n", sep = "")
message("Cluster switcher comparison complete.")
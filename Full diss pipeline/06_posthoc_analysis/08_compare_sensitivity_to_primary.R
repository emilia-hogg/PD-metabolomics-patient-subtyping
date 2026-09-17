# 06_posthoc_analysis/08_compare_sensitivity_to_primary.R
# Re-runnable from a fresh R session.
#
# - summarise how each sensitivity solution compares to the primary integrated SNF solution
# - provide a concise robustness overview for the dissertation chapter
# - report agreement metrics, switch proportions, and contingency tables
#
# -----------------------------------------------------------------------------

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
  if (is.factor(x)) x <- as.character(x)

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
  if (length(hit) == 0L || is.na(hit) || !nzchar(hit)) return(NA_character_)
  hit
}

plot_metric_summary <- function(df, metric, out_file, title) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = branch, y = .data[[metric]], fill = branch)) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")) +
    ggplot2::labs(x = NULL, y = metric, title = title)
  ggplot2::ggsave(out_file, p, width = 7.2, height = 4.6, dpi = 300)
  invisible(p)
}

plot_contingency_heatmap <- function(tab_df, out_file, title) {
  p <- ggplot2::ggplot(tab_df, ggplot2::aes(x = threshold_cluster, y = truth_cluster, fill = count)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = count), size = 3) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")) +
    ggplot2::scale_fill_gradient(low = "grey95", high = "steelblue") +
    ggplot2::labs(x = "Sensitivity cluster", y = "Primary cluster", fill = "Count", title = title)
  ggplot2::ggsave(out_file, p, width = 6.5, height = 5.2, dpi = 300)
  invisible(p)
}

make_contingency_table <- function(truth, pred) {
  truth <- as.integer(truth)
  pred <- as.integer(pred)
  keep <- is.finite(truth) & is.finite(pred)
  truth <- truth[keep]
  pred <- pred[keep]

  tab <- table(truth, pred)
  df <- as.data.frame(tab)
  names(df) <- c("truth_cluster", "threshold_cluster", "count")
  df$truth_cluster <- as.integer(as.character(df$truth_cluster))
  df$threshold_cluster <- as.integer(as.character(df$threshold_cluster))
  df$prop_of_total <- df$count / sum(df$count)
  df
}

load_switcher_table <- function(branch_name) {
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

# -----------------------------------------------------------------------------
# 2. Load primary solution
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
primary_labels <- primary_df$cluster
primary_k <- length(unique(primary_labels))

# -----------------------------------------------------------------------------
# 3. Output directories
# -----------------------------------------------------------------------------
out_root <- get_posthoc_output_dir("compare_sensitivity_to_primary")
plot_dir <- file.path(out_root, "plots")
table_dir <- file.path(out_root, "tables")
diag_dir <- file.path(out_root, "diagnostics")
safe_dir(plot_dir)
safe_dir(table_dir)
safe_dir(diag_dir)

# -----------------------------------------------------------------------------
# 4. Branch definitions
# -----------------------------------------------------------------------------
branch_names <- c(
  "clinical_only",
  "no_LEDD",
  "representation_comparison",
  "top20_combined",
  "top20_23network"
)

# -----------------------------------------------------------------------------
# 5. Compare each sensitivity solution to primary
# -----------------------------------------------------------------------------
comparison_rows <- list()
contingency_objects <- list()
switcher_tables <- list()

for (branch_name in branch_names) {
  message("Comparing primary to: ", branch_name)

  switch_df <- load_switcher_table(branch_name)
  branch_labels <- switch_df$branch_cluster

  comp <- compare_cluster_solutions(primary_labels, branch_labels)
  aligned <- best_label_permutation(primary_labels, branch_labels)$aligned_labels

  cont_df <- make_contingency_table(primary_labels, aligned)
  contingency_objects[[branch_name]] <- cont_df
  switcher_tables[[branch_name]] <- switch_df

  comparison_rows[[branch_name]] <- data.frame(
    branch = branch_name,
    primary_k = primary_k,
    branch_k = length(unique(branch_labels)),
    n_samples = length(primary_labels),
    n_switchers = sum(switch_df$switched, na.rm = TRUE),
    prop_switchers = mean(switch_df$switched, na.rm = TRUE),
    ari = comp$ari,
    nmi = comp$nmi,
    vi = comp$vi,
    exact_direct = comp$exact_direct,
    exact_best = comp$exact_best,
    moved_n = comp$moved_n,
    moved_prop = comp$moved_prop,
    stringsAsFactors = FALSE
  )

  write.csv(cont_df, file.path(table_dir, paste0(branch_name, "_contingency_table.csv")), row.names = FALSE)
  saveRDS(cont_df, file.path(table_dir, paste0(branch_name, "_contingency_table.rds")))
  write.csv(switch_df, file.path(table_dir, paste0(branch_name, "_switchers_reloaded.csv")), row.names = FALSE)
  saveRDS(switch_df, file.path(table_dir, paste0(branch_name, "_switchers_reloaded.rds")))

  plot_contingency_heatmap(
    cont_df,
    file.path(plot_dir, paste0(branch_name, "_contingency_heatmap.png")),
    title = paste0(branch_name, " contingency: primary vs sensitivity")
  )
}

comparison_df <- do.call(rbind, comparison_rows)
rownames(comparison_df) <- NULL
comparison_df <- comparison_df[order(-comparison_df$ari, -comparison_df$exact_best, -comparison_df$n_switchers), ]

write.csv(comparison_df, file.path(out_root, "sensitivity_vs_primary_comparison_summary.csv"), row.names = FALSE)
saveRDS(comparison_df, file.path(out_root, "sensitivity_vs_primary_comparison_summary.rds"))

# -----------------------------------------------------------------------------
# 6. Switcher overlap across branches
# -----------------------------------------------------------------------------
branch_names <- names(switcher_tables)
switch_sets <- lapply(switcher_tables, function(df) df$sample_id[df$switched])

branch_overlap_counts <- matrix(
  0,
  nrow = length(branch_names),
  ncol = length(branch_names),
  dimnames = list(branch_names, branch_names)
)

branch_overlap_jaccard <- matrix(
  NA_real_,
  nrow = length(branch_names),
  ncol = length(branch_names),
  dimnames = list(branch_names, branch_names)
)

for (i in seq_along(branch_names)) {
  for (j in seq_along(branch_names)) {
    a <- switch_sets[[branch_names[i]]]
    b <- switch_sets[[branch_names[j]]]
    branch_overlap_counts[i, j] <- length(intersect(a, b))
    u <- length(union(a, b))
    branch_overlap_jaccard[i, j] <- if (u == 0) NA_real_ else length(intersect(a, b)) / u
  }
}

write.csv(branch_overlap_counts, file.path(out_root, "switcher_overlap_counts.csv"))
write.csv(branch_overlap_jaccard, file.path(out_root, "switcher_overlap_jaccard.csv"))
saveRDS(branch_overlap_counts, file.path(out_root, "switcher_overlap_counts.rds"))
saveRDS(branch_overlap_jaccard, file.path(out_root, "switcher_overlap_jaccard.rds"))

# plot heatmaps from matrices
plot_counts_df <- as.data.frame(as.table(branch_overlap_counts), stringsAsFactors = FALSE)
names(plot_counts_df) <- c("truth_cluster", "threshold_cluster", "count")
plot_contingency_heatmap(
  plot_counts_df,
  file.path(plot_dir, "switcher_overlap_counts_heatmap.png"),
  title = "Overlap of switchers across branches"
)

plot_jacc_df <- as.data.frame(as.table(branch_overlap_jaccard), stringsAsFactors = FALSE)
names(plot_jacc_df) <- c("truth_cluster", "threshold_cluster", "count")
plot_contingency_heatmap(
  plot_jacc_df,
  file.path(plot_dir, "switcher_overlap_jaccard_heatmap.png"),
  title = "Jaccard overlap of switchers across branches"
)

# -----------------------------------------------------------------------------
# 7. Summary plots and export tables
# -----------------------------------------------------------------------------
summary_table <- comparison_df[, c(
  "branch", "primary_k", "branch_k", "n_samples", "n_switchers", "prop_switchers",
  "ari", "nmi", "vi", "exact_direct", "exact_best", "moved_n", "moved_prop"
)]

write.csv(summary_table, file.path(out_root, "comparison_summary_table.csv"), row.names = FALSE)
saveRDS(summary_table, file.path(out_root, "comparison_summary_table.rds"))
saveRDS(contingency_objects, file.path(out_root, "comparison_contingency_objects.rds"))
saveRDS(switcher_tables, file.path(out_root, "comparison_switcher_tables.rds"))

plot_metric_summary(summary_table, "ari", file.path(plot_dir, "comparison_ari.png"), "Adjusted Rand index by branch")
plot_metric_summary(summary_table, "exact_best", file.path(plot_dir, "comparison_exact_best.png"), "Best exact agreement by branch")
plot_metric_summary(summary_table, "prop_switchers", file.path(plot_dir, "comparison_switcher_proportion.png"), "Switcher proportion by branch")

# -----------------------------------------------------------------------------
# 8. Console output
# -----------------------------------------------------------------------------
cat("\nSensitivity-vs-primary comparison summary\n")
print(summary_table)

cat("\nFinished comparison script.\n")
cat("Output directory: ", out_root, "\n", sep = "")
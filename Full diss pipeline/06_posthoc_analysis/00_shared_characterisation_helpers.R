# 06_posthoc_analysis/00_shared_characterisation_helpers.R
# Shared helper functions for the post-hoc analyses.
#
# These functions support loading and validating analysis data and SNF outputs,
# comparing primary and sensitivity networks and cluster solutions, identifying
# cluster switchers, summarising clinical variables by cluster, and summarising
# metabolite classes.
#
# Before running:
# - Set `project_root` to the local project directory.
# - The required analysis, SNF, and sensitivity outputs must already exist for
#   whichever downstream post-hoc script is using these functions.

# -----------------------------------------------------------------------------
# 0. Load config and packages
# -----------------------------------------------------------------------------

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

# -----------------------------------------------------------------------------
# 1. Basic file and directory utilities
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

assert_dir_exists <- function(path, label = basename(path)) {
  if (!dir.exists(path)) {
    stop("Missing required directory for ", label, ": ", path, call. = FALSE)
  }
  invisible(TRUE)
}

save_csv <- function(x, path, row.names = FALSE) {
  write.csv(as.data.frame(x), path, row.names = row.names)
  invisible(path)
}

upper_tri_values <- function(mat) {
  mat <- as.matrix(mat)
  if (nrow(mat) != ncol(mat)) {
    stop("upper_tri_values() requires a square matrix.", call. = FALSE)
  }
  mat[upper.tri(mat)]
}

# -----------------------------------------------------------------------------
# 2. Shared data loading helpers
# -----------------------------------------------------------------------------

load_analysis_objects <- function() {
  assert_file_exists(files$analysis_cohort, "analysis_cohort")
  assert_file_exists(files$clinical_processed, "clinical_processed")
  assert_file_exists(files$metabolomics_processed, "metabolomics_processed")

  analysis_cohort <- readRDS(files$analysis_cohort)
  clinical_obj <- readRDS(files$clinical_processed)
  metab_obj <- readRDS(files$metabolomics_processed)

  if (!all(c("pd_model", "metab_pd", "sample_ids") %in% names(analysis_cohort))) {
    stop(
      "analysis_cohort.rds must contain pd_model, metab_pd, and sample_ids.",
      call. = FALSE
    )
  }

  pd_model <- as.data.frame(analysis_cohort$pd_model)
  metab_pd <- as.data.frame(analysis_cohort$metab_pd)
  sample_ids <- as.character(analysis_cohort$sample_ids)

  if (anyDuplicated(sample_ids)) {
    stop("analysis_cohort contains duplicated sample_ids.", call. = FALSE)
  }

  if (!identical(rownames(pd_model), sample_ids)) {
    stop("Row order mismatch: pd_model does not match sample_ids.", call. = FALSE)
  }
  if (!identical(rownames(metab_pd), sample_ids)) {
    stop("Row order mismatch: metab_pd does not match sample_ids.", call. = FALSE)
  }

  list(
    analysis_cohort = analysis_cohort,
    clinical_obj = clinical_obj,
    metab_obj = metab_obj,
    pd_model = pd_model,
    metab_pd = metab_pd,
    sample_ids = sample_ids
  )
}

# -----------------------------------------------------------------------------
# 3. SNF output loaders
# -----------------------------------------------------------------------------

load_primary_snf_outputs <- function() {
  main_snf_dir <- file.path(paths$snf, "main_four_network")
  cluster_dir <- file.path(main_snf_dir, "cluster_solution")

  fused_file <- file.path(main_snf_dir, "primary_four_network_fused_network.rds")
  sample_file <- file.path(main_snf_dir, "primary_four_network_sample_ids.rds")
  labels_file <- file.path(cluster_dir, "best_k_cluster_labels.rds")

  assert_file_exists(fused_file, "primary fused network")
  assert_file_exists(sample_file, "primary sample ids")
  assert_file_exists(labels_file, "primary cluster labels")

  list(
    main_snf_dir = main_snf_dir,
    cluster_dir = cluster_dir,
    fused_network = readRDS(fused_file),
    sample_ids = as.character(readRDS(sample_file)),
    cluster_labels = readRDS(labels_file)
  )
}

load_sensitivity_solution <- function(branch_dir) {
  branch_dir <- normalizePath(branch_dir, mustWork = FALSE)
  assert_dir_exists(branch_dir, basename(branch_dir))

  fused_candidates <- c(
    file.path(branch_dir, "fused_network.rds"),
    file.path(branch_dir, "branch_fused_network.rds"),
    file.path(branch_dir, "network.rds")
  )
  labels_candidates <- c(
    file.path(branch_dir, "best_k_cluster_labels.rds"),
    file.path(branch_dir, "cluster_solution", "best_k_cluster_labels.rds"),
    file.path(branch_dir, "cluster_labels.rds")
  )
  sample_candidates <- c(
    file.path(branch_dir, "sample_ids.rds"),
    file.path(branch_dir, "branch_sample_ids.rds"),
    file.path(branch_dir, "cluster_solution", "sample_ids.rds")
  )

  fused_file <- fused_candidates[file.exists(fused_candidates)][1]
  labels_file <- labels_candidates[file.exists(labels_candidates)][1]
  sample_file <- sample_candidates[file.exists(sample_candidates)][1]

  if (is.na(fused_file) || is.null(fused_file)) {
    stop("Could not find a fused-network RDS in branch dir: ", branch_dir, call. = FALSE)
  }
  if (is.na(labels_file) || is.null(labels_file)) {
    stop("Could not find cluster labels in branch dir: ", branch_dir, call. = FALSE)
  }
  if (is.na(sample_file) || is.null(sample_file)) {
    stop("Could not find sample IDs in branch dir: ", branch_dir, call. = FALSE)
  }

  list(
    branch_dir = branch_dir,
    fused_network = readRDS(fused_file),
    sample_ids = as.character(readRDS(sample_file)),
    cluster_labels = readRDS(labels_file)
  )
}

# -----------------------------------------------------------------------------
# 4. Matrix alignment and validation
# -----------------------------------------------------------------------------

assert_identical_order <- function(x, y, label_x = deparse(substitute(x)), label_y = deparse(substitute(y))) {
  if (!identical(x, y)) {
    stop("Order mismatch between ", label_x, " and ", label_y, ".", call. = FALSE)
  }
  invisible(TRUE)
}

check_symmetric_matrix <- function(mat, label = deparse(substitute(mat))) {
  mat <- as.matrix(mat)
  if (nrow(mat) != ncol(mat)) {
    stop(label, " must be square.", call. = FALSE)
  }
  if (anyNA(mat)) {
    stop(label, " contains NA values.", call. = FALSE)
  }
  if (!isTRUE(all.equal(mat, t(mat)))) {
    stop(label, " is not symmetric.", call. = FALSE)
  }
  invisible(TRUE)
}

reorder_matrix_to_ids <- function(mat, ids) {
  mat <- as.matrix(mat)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("Matrix must have row and column names to reorder.", call. = FALSE)
  }
  ids <- as.character(ids)
  missing_ids <- setdiff(ids, rownames(mat))
  if (length(missing_ids) > 0) {
    stop(
      "Matrix is missing required IDs. First missing ID: ",
      missing_ids[1],
      call. = FALSE
    )
  }
  mat <- mat[ids, ids, drop = FALSE]
  if (!identical(rownames(mat), ids) || !identical(colnames(mat), ids)) {
    stop("Failed to reorder matrix to requested IDs.", call. = FALSE)
  }
  mat
}

# -----------------------------------------------------------------------------
# 5. Permutation / label alignment helpers
# -----------------------------------------------------------------------------

all_permutations <- function(x) {
  x <- as.character(x)
  n <- length(x)

  if (n <= 1) {
    return(matrix(x, nrow = 1))
  }

  out <- list()

  permute_rec <- function(vec, prefix = character(0)) {
    if (length(vec) == 0) {
      out[[length(out) + 1L]] <<- prefix
      return(invisible(NULL))
    }
    for (i in seq_along(vec)) {
      permute_rec(vec[-i], c(prefix, vec[i]))
    }
  }

  permute_rec(x)
  do.call(rbind, out)
}

contingency_table <- function(main_labels, branch_labels) {
  main_labels <- as.integer(main_labels)
  branch_labels <- as.integer(branch_labels)
  table(main_labels, branch_labels)
}

best_label_permutation <- function(main_labels, branch_labels) {
  main_labels <- as.integer(main_labels)
  branch_labels <- as.integer(branch_labels)

  if (length(main_labels) != length(branch_labels)) {
    stop("Label vectors must have the same length.", call. = FALSE)
  }

  tab <- table(main_labels, branch_labels)
  main_lvls <- rownames(tab)
  branch_lvls <- colnames(tab)

  if (length(branch_lvls) > 8) {
    warning(
      "best_label_permutation() uses brute force and may be slow for > 8 clusters."
    )
  }

  branch_perms <- all_permutations(branch_lvls)

  best_score <- -Inf
  best_map <- NULL

  for (i in seq_len(nrow(branch_perms))) {
    perm <- branch_perms[i, ]
    names(perm) <- branch_lvls

    aligned <- perm[as.character(branch_labels)]
    aligned <- suppressWarnings(as.integer(aligned))
    if (anyNA(aligned)) {
      next
    }

    score <- sum(main_labels == aligned)
    if (score > best_score) {
      best_score <- score
      best_map <- perm
    }
  }

  if (is.null(best_map)) {
    stop("Could not determine a valid label permutation.", call. = FALSE)
  }

  aligned_labels <- suppressWarnings(as.integer(best_map[as.character(branch_labels)]))
  if (anyNA(aligned_labels)) {
    stop("Aligned labels contain NA values.", call. = FALSE)
  }

  list(
    aligned_labels = aligned_labels,
    mapping = best_map,
    contingency = tab,
    exact_best = mean(main_labels == aligned_labels)
  )
}

# -----------------------------------------------------------------------------
# 6. Base-R clustering comparison metrics
# -----------------------------------------------------------------------------

cluster_entropy <- function(labels) {
  labels <- as.integer(labels)
  p <- prop.table(table(labels))
  p <- p[p > 0]
  -sum(p * log(p))
}

mutual_information_from_table <- function(tab) {
  tab <- as.matrix(tab)
  n <- sum(tab)
  if (n <= 0) {
    return(0)
  }

  row_sums <- rowSums(tab)
  col_sums <- colSums(tab)
  mi <- 0

  for (i in seq_len(nrow(tab))) {
    for (j in seq_len(ncol(tab))) {
      nij <- tab[i, j]
      if (nij > 0) {
        mi <- mi + (nij / n) * log((nij * n) / (row_sums[i] * col_sums[j]))
      }
    }
  }

  mi
}

normalized_mutual_information <- function(main_labels, branch_labels) {
  tab <- contingency_table(main_labels, branch_labels)
  mi <- mutual_information_from_table(tab)
  hx <- cluster_entropy(main_labels)
  hy <- cluster_entropy(branch_labels)

  denom <- sqrt(hx * hy)
  if (!is.finite(denom) || denom <= 0) {
    return(0)
  }
  mi / denom
}

variation_of_information <- function(main_labels, branch_labels) {
  tab <- contingency_table(main_labels, branch_labels)
  h_main <- cluster_entropy(main_labels)
  h_branch <- cluster_entropy(branch_labels)
  mi <- mutual_information_from_table(tab)
  h_main + h_branch - 2 * mi
}

adjusted_rand_index <- function(main_labels, branch_labels) {
  tab <- contingency_table(main_labels, branch_labels)
  tab <- as.matrix(tab)
  n <- sum(tab)

  if (n <= 1) {
    return(0)
  }

  choose2 <- function(x) if (x < 2) 0 else choose(x, 2)

  sum_ij <- sum(vapply(tab, choose2, numeric(1)))
  sum_i <- sum(vapply(rowSums(tab), choose2, numeric(1)))
  sum_j <- sum(vapply(colSums(tab), choose2, numeric(1)))
  total <- choose2(n)

  expected <- (sum_i * sum_j) / total
  max_index <- 0.5 * (sum_i + sum_j)
  denom <- max_index - expected

  if (denom == 0) {
    return(0)
  }

  (sum_ij - expected) / denom
}

compare_cluster_solutions <- function(main_labels, branch_labels) {
  main_labels <- as.integer(main_labels)
  branch_labels <- as.integer(branch_labels)

  if (length(main_labels) != length(branch_labels)) {
    stop("Cluster vectors must have the same length.", call. = FALSE)
  }

  alignment <- best_label_permutation(main_labels, branch_labels)
  aligned_branch <- alignment$aligned_labels

  data.frame(
    k_main = length(unique(main_labels)),
    k_branch = length(unique(branch_labels)),
    ari = adjusted_rand_index(main_labels, branch_labels),
    nmi = normalized_mutual_information(main_labels, branch_labels),
    vi = variation_of_information(main_labels, branch_labels),
    exact_direct = mean(main_labels == branch_labels),
    exact_best = mean(main_labels == aligned_branch),
    moved_n = sum(main_labels != aligned_branch),
    moved_prop = mean(main_labels != aligned_branch),
    stringsAsFactors = FALSE
  )
}

cluster_switcher_table <- function(sample_ids, main_labels, branch_labels) {
  sample_ids <- as.character(sample_ids)
  main_labels <- as.integer(main_labels)
  branch_labels <- as.integer(branch_labels)

  if (length(sample_ids) != length(main_labels) || length(main_labels) != length(branch_labels)) {
    stop("sample_ids, main_labels, and branch_labels must have the same length.", call. = FALSE)
  }

  alignment <- best_label_permutation(main_labels, branch_labels)
  aligned_branch <- alignment$aligned_labels

  data.frame(
    sample_id = sample_ids,
    main_cluster = main_labels,
    branch_cluster = branch_labels,
    aligned_cluster = aligned_branch,
    switched = main_labels != aligned_branch,
    stringsAsFactors = FALSE
  )
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

summarise_switchers_by_cluster <- function(df, cluster_col, label_prefix) {
  if (!cluster_col %in% names(df)) {
    stop("cluster_col not found in df: ", cluster_col, call. = FALSE)
  }
  if (!"switched" %in% names(df)) {
    stop("df must contain a switched column.", call. = FALSE)
  }

  cl <- as.integer(as.character(df[[cluster_col]]))
  if (anyNA(cl)) {
    stop("Cluster column contains values that could not be converted to integers: ", cluster_col, call. = FALSE)
  }

  tab_n <- table(cl)
  tab_sw <- tapply(df$switched, cl, function(x) sum(x, na.rm = TRUE))
  tab_sw <- tab_sw[names(tab_n)]
  tab_sw[is.na(tab_sw)] <- 0L

  out <- data.frame(
    cluster = as.integer(names(tab_n)),
    n_samples = as.integer(tab_n),
    n_switchers = as.integer(tab_sw),
    prop_switchers = as.numeric(tab_sw / as.integer(tab_n)),
    cluster_label = paste0(label_prefix, "_", as.integer(names(tab_n))),
    stringsAsFactors = FALSE
  )

  out[order(out$cluster), , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 7. Base-R network comparison metrics
# -----------------------------------------------------------------------------

compare_networks <- function(W_main, W_branch) {
  W_main <- as.matrix(W_main)
  W_branch <- as.matrix(W_branch)

  if (!identical(rownames(W_main), rownames(W_branch))) {
    stop("Network row names do not match.", call. = FALSE)
  }
  if (!identical(colnames(W_main), colnames(W_branch))) {
    stop("Network column names do not match.", call. = FALSE)
  }

  check_symmetric_matrix(W_main, "W_main")
  check_symmetric_matrix(W_branch, "W_branch")

  a <- upper_tri_values(W_main)
  b <- upper_tri_values(W_branch)

  data.frame(
    n_samples = nrow(W_main),
    n_edges = length(a),
    pearson = cor(a, b, method = "pearson"),
    spearman = cor(a, b, method = "spearman"),
    rmse = sqrt(mean((a - b)^2)),
    mae = mean(abs(a - b)),
    max_abs_diff = max(abs(a - b)),
    median_abs_diff = stats::median(abs(a - b)),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 8. Simple cluster-wise summaries for later scripts
# -----------------------------------------------------------------------------

summarise_continuous_by_cluster <- function(data, cluster_col, vars, digits = 3) {
  data <- as.data.frame(data)
  if (!cluster_col %in% names(data)) {
    stop("cluster_col not found in data: ", cluster_col, call. = FALSE)
  }

  vars <- intersect(vars, names(data))
  if (length(vars) == 0) {
    stop("No requested variables found in data.", call. = FALSE)
  }

  cl <- factor(data[[cluster_col]])
  out <- list()

  for (v in vars) {
    x <- data[[v]]
    if (!is.numeric(x)) {
      next
    }

    tmp <- aggregate(x, list(cluster = cl), function(z) {
      c(
        n = sum(is.finite(z)),
        mean = mean(z, na.rm = TRUE),
        sd = stats::sd(z, na.rm = TRUE),
        median = stats::median(z, na.rm = TRUE),
        q25 = stats::quantile(z, 0.25, na.rm = TRUE, names = FALSE),
        q75 = stats::quantile(z, 0.75, na.rm = TRUE, names = FALSE)
      )
    })

    tmp <- do.call(data.frame, tmp)
    tmp$variable <- v
    out[[v]] <- tmp
  }

  res <- do.call(rbind, out)
  rownames(res) <- NULL

  if (!is.null(digits)) {
    num_cols <- vapply(res, is.numeric, logical(1))
    res[num_cols] <- lapply(res[num_cols], function(z) round(z, digits))
  }

  res
}

cohens_d <- function(x, g) {
  x <- as.numeric(x)
  g <- as.factor(g)
  lv <- levels(g)

  if (length(lv) != 2) {
    return(NA_real_)
  }

  x1 <- x[g == lv[1]]
  x2 <- x[g == lv[2]]

  m1 <- mean(x1, na.rm = TRUE)
  m2 <- mean(x2, na.rm = TRUE)
  s1 <- stats::sd(x1, na.rm = TRUE)
  s2 <- stats::sd(x2, na.rm = TRUE)
  n1 <- sum(is.finite(x1))
  n2 <- sum(is.finite(x2))

  pooled <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(pooled) || pooled == 0) {
    return(NA_real_)
  }

  (m2 - m1) / pooled
}

# -----------------------------------------------------------------------------
# 9. Output directory helpers
# -----------------------------------------------------------------------------

get_posthoc_output_dir <- function(subdir = NULL) {
  base <- paths$posthoc
  safe_dir(base)

  if (is.null(subdir) || identical(subdir, "")) {
    return(base)
  }

  out <- file.path(base, subdir)
  safe_dir(out)
  out
}

# -----------------------------------------------------------------------------
# 10. Small convenience wrapper for writing paired outputs
# -----------------------------------------------------------------------------

write_comparison_outputs <- function(out_dir, prefix, metrics_df, switcher_df = NULL, contingency = NULL) {
  safe_dir(out_dir)
  save_csv(metrics_df, file.path(out_dir, paste0(prefix, "_comparison_metrics.csv")))
  if (!is.null(switcher_df)) {
    save_csv(switcher_df, file.path(out_dir, paste0(prefix, "_switchers.csv")))
  }
  if (!is.null(contingency)) {
    save_csv(as.data.frame.matrix(contingency), file.path(out_dir, paste0(prefix, "_contingency_table.csv")))
  }
  invisible(TRUE)
}

# -------------------------------------------------------------------
# Class-level summary with annotation carried through.
# Replaces summarise_classes() in 04 and provides the primary equivalent for 02.
#
# direction = "lowest"  -> strongest classes first (smallest median p-adj)
# direction = "highest" -> reproduces the original (weakest first)
# -------------------------------------------------------------------
summarise_classes_annotated <- function(feature_table,
                                        class_col   = "SUB_PATHWAY",
                                        p_col       = "welch_p_adj",
                                        diff_col    = "mean_diff",
                                        super_col   = "SUPER_PATHWAY",
                                        alpha       = 0.05,
                                        direction   = c("lowest", "highest"),
                                        top_n       = 15L) {
  direction <- match.arg(direction)
  df <- as.data.frame(feature_table, stringsAsFactors = FALSE)

  if (!all(c(class_col, p_col) %in% names(df))) {
    stop("feature_table must contain ", class_col, " and ", p_col, call. = FALSE)
  }

  keep <- !is.na(df[[class_col]]) &
    nzchar(as.character(df[[class_col]])) &
    is.finite(suppressWarnings(as.numeric(df[[p_col]])))
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0L) return(data.frame())

  has_super <- super_col %in% names(df)
  has_diff  <- diff_col  %in% names(df)

  sp <- split(seq_len(nrow(df)), as.character(df[[class_col]]))

  out <- do.call(rbind, lapply(names(sp), function(g) {
    sub <- df[sp[[g]], , drop = FALSE]
    p   <- suppressWarnings(as.numeric(sub[[p_col]]))
    n_sig <- sum(p < alpha, na.rm = TRUE)

    # Sub-pathways map one-to-one onto super-pathways in the Metabolon
    # annotation, but take the most frequent non-missing value defensively
    # and flag any class where that assumption fails.
    super_val <- NA_character_
    n_super   <- NA_integer_
    if (has_super) {
      s <- as.character(sub[[super_col]])
      s <- s[!is.na(s) & nzchar(s)]
      n_super <- length(unique(s))
      if (length(s) > 0L) {
        super_val <- names(sort(table(s), decreasing = TRUE))[1]
      }
    }

    data.frame(
      class               = g,
      super_pathway       = super_val,
      n_super_pathways    = n_super,
      n_metabolites       = nrow(sub),
      n_significant       = n_sig,
      metabolites_display = if (n_sig > 0) paste0(nrow(sub), " (", n_sig, ")") else as.character(nrow(sub)),
      median_p_adj        = stats::median(p, na.rm = TRUE),
      mean_p_adj          = mean(p, na.rm = TRUE),
      min_p_adj           = min(p, na.rm = TRUE),
      median_abs_effect   = if (has_diff) stats::median(abs(suppressWarnings(as.numeric(sub[[diff_col]]))), na.rm = TRUE) else NA_real_,
      stringsAsFactors    = FALSE
    )
  }))

  ord <- if (direction == "lowest") order(out$median_p_adj) else order(-out$median_p_adj)
  out <- out[ord, , drop = FALSE]

  out$median_p_display <- format_p_display(out$median_p_adj)
  out$mean_p_display   <- format_p_display(out$mean_p_adj)

  rownames(out) <- NULL
  if (!is.null(top_n) && nrow(out) > top_n) out <- out[seq_len(top_n), , drop = FALSE]
  out
}

format_p_display <- function(p, digits = 3) {
  ifelse(is.na(p), "-",
         ifelse(p >= 0.001,
                formatC(p, format = "f", digits = digits),
                formatC(p, format = "e", digits = 1)))
}

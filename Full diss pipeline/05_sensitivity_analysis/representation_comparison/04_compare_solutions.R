# 04_compare_solutions.R
# This script compares the two-network representation-comparison solution with
# the primary four-network SNF solution. It compares the fused networks and
# cluster assignments, and identifies patients whose cluster differs between
# the two solutions.
#
# Before running:
# - Set `project_root` to the local project directory.
# - The primary four-network SNF and cluster solution must already exist.
# - The representation-comparison SNF and cluster solution must already exist.
#
# Outputs are saved under
# `05_sensitivity_analysis/outputs/representation_comparison/comparison/`,
# including network similarity, cluster agreement, matched labels, and
# switcher tables.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

if (!requireNamespace("SNFtool", quietly = TRUE)) {
  stop("SNFtool is required but is not available.")
}

# -------------------------------------------------------------------
# Input / output locations
# -------------------------------------------------------------------
main_snf_dir <- file.path(paths$snf, "main_four_network")
main_cluster_dir <- file.path(main_snf_dir, "cluster_solution")

branch_dir <- file.path(paths$sensitivity, "representation_comparison")
branch_cluster_dir <- file.path(branch_dir, "cluster_solution")
out_dir <- file.path(branch_dir, "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

main_fused_file <- file.path(main_snf_dir, "primary_four_network_fused_network.rds")
main_best_labels_file <- file.path(main_cluster_dir, "best_k_cluster_labels.rds")
main_run_record_file <- file.path(main_cluster_dir, "cluster_run_record.rds")

branch_results_file <- file.path(branch_cluster_dir, "representation_2network_cluster_results.rds")
branch_best_labels_file <- file.path(branch_cluster_dir, "best_k_cluster_labels.rds")
branch_run_record_file <- file.path(branch_cluster_dir, "cluster_run_record.rds")

required_files <- c(
  main_fused_file,
  main_best_labels_file,
  main_run_record_file,
  branch_results_file,
  branch_best_labels_file,
  branch_run_record_file
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required file(s): ", paste(missing_files, collapse = "; "))
}

# -------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------
compare_networks_upper_triangle <- function(W_main, W_branch) {
  ids <- intersect(rownames(W_main), rownames(W_branch))
  if (length(ids) < 2L) stop("Not enough shared samples to compare fused networks.")

  W_main <- as.matrix(W_main[ids, ids, drop = FALSE])
  W_branch <- as.matrix(W_branch[ids, ids, drop = FALSE])

  tri_main <- W_main[upper.tri(W_main)]
  tri_branch <- W_branch[upper.tri(W_branch)]
  ok <- is.finite(tri_main) & is.finite(tri_branch)
  tri_main <- tri_main[ok]
  tri_branch <- tri_branch[ok]

  data.frame(
    metric = c(
      "n_samples",
      "n_edges",
      "pearson",
      "spearman",
      "rmse",
      "mae",
      "max_abs_diff",
      "median_abs_diff"
    ),
    value = c(
      length(ids),
      length(tri_main),
      suppressWarnings(cor(tri_main, tri_branch, method = "pearson")),
      suppressWarnings(cor(tri_main, tri_branch, method = "spearman")),
      sqrt(mean((tri_main - tri_branch)^2)),
      mean(abs(tri_main - tri_branch)),
      max(abs(tri_main - tri_branch)),
      median(abs(tri_main - tri_branch))
    ),
    stringsAsFactors = FALSE
  )
}

entropy_base <- function(x) {
  p <- prop.table(table(x))
  p <- p[p > 0]
  -sum(p * log(p))
}

mutual_information_base <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  if (n == 0) return(0)

  px <- rowSums(tab) / n
  py <- colSums(tab) / n

  mi <- 0
  for (i in seq_len(nrow(tab))) {
    for (j in seq_len(ncol(tab))) {
      nij <- tab[i, j]
      if (nij > 0) {
        pij <- nij / n
        mi <- mi + pij * log(pij / (px[i] * py[j]))
      }
    }
  }
  mi
}

adjusted_rand_index_base <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  if (n < 2) return(NA_real_)

  comb2 <- function(z) ifelse(z < 2, 0, choose(z, 2))

  sum_ij <- sum(comb2(tab))
  sum_i <- sum(comb2(rowSums(tab)))
  sum_j <- sum(comb2(colSums(tab)))
  total <- comb2(n)
  expected <- (sum_i * sum_j) / total
  max_index <- 0.5 * (sum_i + sum_j)

  if (max_index == expected) return(0)
  (sum_ij - expected) / (max_index - expected)
}

normalized_mutual_information_base <- function(x, y) {
  hx <- entropy_base(x)
  hy <- entropy_base(y)
  mi <- mutual_information_base(x, y)
  if ((hx + hy) == 0) return(1)
  2 * mi / (hx + hy)
}

variation_of_information_base <- function(x, y) {
  hx <- entropy_base(x)
  hy <- entropy_base(y)
  mi <- mutual_information_base(x, y)
  hx + hy - 2 * mi
}

align_labels_best_permutation <- function(main_df, branch_df) {
  merged <- merge(main_df, branch_df, by = "sample_id", suffixes = c("_main", "_branch"), sort = FALSE)
  if (nrow(merged) == 0L) stop("No overlapping sample IDs between main and branch cluster labels.")

  main_labels <- as.integer(factor(merged$cluster_main, levels = sort(unique(merged$cluster_main))))
  branch_labels <- as.integer(factor(merged$cluster_branch, levels = sort(unique(merged$cluster_branch))))

  k_main <- length(unique(main_labels))
  k_branch <- length(unique(branch_labels))

  make_summary <- function(ari_val, nmi_val, vi_val, exact_direct, exact_best, moved_n, moved_prop) {
    data.frame(
      metric = c("k_main", "k_branch", "ari", "nmi", "vi", "exact_direct", "exact_best", "moved_n", "moved_prop"),
      value = c(k_main, k_branch, ari_val, nmi_val, vi_val, exact_direct, exact_best, moved_n, moved_prop),
      stringsAsFactors = FALSE
    )
  }

  ari <- function(x, y) {
    tab <- table(x, y)
    n <- sum(tab)
    if (n < 2) return(NA_real_)
    comb2 <- function(z) ifelse(z < 2, 0, choose(z, 2))
    sum_ij <- sum(comb2(tab))
    sum_i <- sum(comb2(rowSums(tab)))
    sum_j <- sum(comb2(colSums(tab)))
    total <- comb2(n)
    expected <- (sum_i * sum_j) / total
    max_index <- 0.5 * (sum_i + sum_j)
    if (max_index == expected) return(0)
    (sum_ij - expected) / (max_index - expected)
  }

  entropy <- function(x) {
    p <- prop.table(table(x))
    p <- p[p > 0]
    -sum(p * log(p))
  }

  mi <- function(x, y) {
    tab <- table(x, y)
    n <- sum(tab)
    if (n == 0) return(0)
    px <- rowSums(tab) / n
    py <- colSums(tab) / n
    ans <- 0
    for (i in seq_len(nrow(tab))) {
      for (j in seq_len(ncol(tab))) {
        nij <- tab[i, j]
        if (nij > 0) {
          pij <- nij / n
          ans <- ans + pij * log(pij / (px[i] * py[j]))
        }
      }
    }
    ans
  }

  nmi <- function(x, y) {
    hx <- entropy(x)
    hy <- entropy(y)
    m <- mi(x, y)
    if ((hx + hy) == 0) return(1)
    2 * m / (hx + hy)
  }

  vi <- function(x, y) {
    hx <- entropy(x)
    hy <- entropy(y)
    m <- mi(x, y)
    hx + hy - 2 * m
  }

  if (k_main == k_branch && k_main %in% c(2L, 3L)) {
    perms <- if (k_main == 2L) {
      list(c(1L, 2L), c(2L, 1L))
    } else {
      list(
        c(1L, 2L, 3L),
        c(1L, 3L, 2L),
        c(2L, 1L, 3L),
        c(2L, 3L, 1L),
        c(3L, 1L, 2L),
        c(3L, 2L, 1L)
      )
    }

    perm_tbl <- do.call(rbind, lapply(seq_along(perms), function(i) {
      perm <- perms[[i]]
      aligned <- perm[branch_labels]
      data.frame(
        permutation_id = i,
        mapping = paste(seq_len(k_main), "->", perm, collapse = "; "),
        exact = mean(main_labels == aligned),
        stringsAsFactors = FALSE
      )
    }))

    best_idx <- which.max(perm_tbl$exact)
    best_perm <- perms[[best_idx]]
    aligned_branch <- best_perm[branch_labels]

    summary_table <- make_summary(
      ari_val = ari(main_labels, aligned_branch),
      nmi_val = nmi(main_labels, aligned_branch),
      vi_val = vi(main_labels, aligned_branch),
      exact_direct = mean(main_labels == branch_labels),
      exact_best = mean(main_labels == aligned_branch),
      moved_n = sum(main_labels != aligned_branch),
      moved_prop = mean(main_labels != aligned_branch)
    )

    switchers <- merged[main_labels != aligned_branch, c("sample_id", "cluster_main", "cluster_branch"), drop = FALSE]
    names(switchers) <- c("sample_id", "main_cluster", "branch_cluster_raw")
    switchers$branch_cluster_aligned <- aligned_branch[main_labels != aligned_branch]

    return(list(
      merged = data.frame(
        sample_id = merged$sample_id,
        main_cluster = main_labels,
        branch_cluster_raw = branch_labels,
        branch_cluster_aligned = aligned_branch,
        moved = main_labels != aligned_branch,
        stringsAsFactors = FALSE
      ),
      summary = summary_table,
      permutation_check = within(perm_tbl, is_best <- seq_len(nrow(perm_tbl)) == best_idx),
      switchers = switchers
    ))
  }

  summary_table <- make_summary(
    ari_val = ari(main_labels, branch_labels),
    nmi_val = nmi(main_labels, branch_labels),
    vi_val = vi(main_labels, branch_labels),
    exact_direct = NA_real_,
    exact_best = NA_real_,
    moved_n = NA_real_,
    moved_prop = NA_real_
  )

  list(
    merged = data.frame(
      sample_id = merged$sample_id,
      main_cluster = main_labels,
      branch_cluster_raw = branch_labels,
      branch_cluster_aligned = NA_integer_,
      moved = NA,
      stringsAsFactors = FALSE
    ),
    summary = summary_table,
    permutation_check = data.frame(),
    switchers = data.frame()
  )
}

# -------------------------------------------------------------------
# Load main and branch outputs
# -------------------------------------------------------------------
main_fused <- readRDS(main_fused_file)
main_labels <- readRDS(main_best_labels_file)
main_run_record <- readRDS(main_run_record_file)

branch_obj <- readRDS(branch_results_file)
branch_best_labels <- readRDS(branch_best_labels_file)
branch_run_record <- readRDS(branch_run_record_file)

if (!all(c("sample_id", "cluster") %in% names(main_labels))) {
  stop("Main best_k_cluster_labels.rds must contain sample_id and cluster.")
}
if (!all(c("sample_id", "cluster") %in% names(branch_best_labels))) {
  stop("Branch best_k_cluster_labels.rds must contain sample_id and cluster.")
}

main_labels$sample_id <- as.character(main_labels$sample_id)
main_labels$cluster <- as.integer(main_labels$cluster)
branch_best_labels$sample_id <- as.character(branch_best_labels$sample_id)
branch_best_labels$cluster <- as.integer(branch_best_labels$cluster)

main_fused <- as.matrix(main_fused)
branch_fused <- as.matrix(branch_obj$fused_network)

if (is.null(rownames(main_fused)) || is.null(colnames(main_fused))) {
  stop("Main fused network must have dimnames.")
}
if (is.null(rownames(branch_fused)) || is.null(colnames(branch_fused))) {
  branch_ids <- as.character(branch_obj$sample_ids)
  rownames(branch_fused) <- branch_ids
  colnames(branch_fused) <- branch_ids
}

# -------------------------------------------------------------------
# Compare fused networks
# -------------------------------------------------------------------
network_summary <- compare_networks_upper_triangle(main_fused, branch_fused)

# -------------------------------------------------------------------
# Compare cluster labels
# -------------------------------------------------------------------
main_cluster_vec <- setNames(main_labels$cluster, main_labels$sample_id)
branch_cluster_vec <- setNames(branch_best_labels$cluster, branch_best_labels$sample_id)

shared_ids <- intersect(names(main_cluster_vec), names(branch_cluster_vec))
if (length(shared_ids) < 2L) {
  stop("Not enough shared samples between main and branch cluster labels.")
}

main_cluster_vec <- main_cluster_vec[shared_ids]
branch_cluster_vec <- branch_cluster_vec[shared_ids]

comparison_obj <- align_labels_best_permutation(
  data.frame(sample_id = shared_ids, cluster_main = as.integer(main_cluster_vec), stringsAsFactors = FALSE),
  data.frame(sample_id = shared_ids, cluster_branch = as.integer(branch_cluster_vec), stringsAsFactors = FALSE)
)

matched_df <- comparison_obj$merged

# -------------------------------------------------------------------
# Save outputs
# -------------------------------------------------------------------
prefix <- "representation_2network_vs_main_four_network"

write.csv(network_summary, file.path(out_dir, paste0(prefix, "_network_similarity.csv")), row.names = FALSE)
write.csv(comparison_obj$summary, file.path(out_dir, paste0(prefix, "_summary.csv")), row.names = FALSE)
write.csv(matched_df, file.path(out_dir, paste0(prefix, "_matched_labels.csv")), row.names = FALSE)
write.csv(comparison_obj$switchers, file.path(out_dir, paste0(prefix, "_switchers.csv")), row.names = FALSE)
write.csv(comparison_obj$permutation_check, file.path(out_dir, paste0(prefix, "_permutation_check.csv")), row.names = FALSE)

saveRDS(network_summary, file.path(out_dir, paste0(prefix, "_network_similarity.rds")))
saveRDS(comparison_obj$summary, file.path(out_dir, paste0(prefix, "_summary.rds")))
saveRDS(matched_df, file.path(out_dir, paste0(prefix, "_matched_labels.rds")))
saveRDS(comparison_obj$switchers, file.path(out_dir, paste0(prefix, "_switchers.rds")))
saveRDS(comparison_obj$permutation_check, file.path(out_dir, paste0(prefix, "_permutation_check.rds")))

message("Representation comparison complete.")
message("Saved outputs to: ", out_dir)
message("Main best k: ", unique(main_run_record$value[main_run_record$metric == "best_k"]))
message("Branch best k: ", unique(branch_run_record$value[branch_run_record$metric == "best_k"]))

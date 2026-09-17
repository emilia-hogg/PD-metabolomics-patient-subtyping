# =========================================================
# 03_compare_solutions.R
# clinical_only branch vs saved main four-network solution
# Compares the 3-network clinical-only solution to the main 4-network solution
# =========================================================

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

main_snf_dir <- file.path(project_root, "04_snf", "outputs", "main_four_network")
main_cluster_dir <- file.path(main_snf_dir, "cluster_solution")
branch_dir <- file.path(project_root, "05_sensitivity_analysis", "outputs", "clinical_only")
out_dir <- file.path(branch_dir, "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

main_fused_file <- file.path(main_snf_dir, "primary_four_network_fused_network.rds")
main_best_labels_file <- file.path(main_cluster_dir, "best_k_cluster_labels.rds")
branch_results_file <- file.path(branch_dir, "clinical_only_cluster_results.rds")

if (!file.exists(main_fused_file)) stop("Missing main fused network file: ", main_fused_file)
if (!file.exists(main_best_labels_file)) stop("Missing main best-k labels file: ", main_best_labels_file)
if (!file.exists(branch_results_file)) stop("Missing branch cluster results file: ", branch_results_file)

main_fused <- readRDS(main_fused_file)
main_labels <- readRDS(main_best_labels_file)
branch_obj <- readRDS(branch_results_file)

if (!all(c("fused_network", "sample_ids", "cluster_labels") %in% names(branch_obj))) {
  stop("clinical_only_cluster_results.rds must contain fused_network, sample_ids, and cluster_labels.")
}

branch_fused <- as.matrix(branch_obj$fused_network)
branch_sample_ids <- as.character(branch_obj$sample_ids)
branch_labels <- as.data.frame(branch_obj$cluster_labels, stringsAsFactors = FALSE)

if (!all(c("sample_id", "cluster") %in% names(main_labels))) {
  stop("Main best_k_cluster_labels.rds must contain sample_id and cluster.")
}
if (!all(c("sample_id", "cluster") %in% names(branch_labels))) {
  stop("Branch cluster_labels must contain sample_id and cluster.")
}

main_labels$sample_id <- as.character(main_labels$sample_id)
main_labels$cluster <- as.integer(main_labels$cluster)
branch_labels$sample_id <- as.character(branch_labels$sample_id)
branch_labels$cluster <- as.integer(branch_labels$cluster)

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
    metric = c("n_samples", "n_edges", "pearson", "spearman", "rmse", "mae", "max_abs_diff", "median_abs_diff"),
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

align_and_compare_labels <- function(main_df, branch_df) {
  merged <- merge(main_df, branch_df, by = "sample_id", suffixes = c("_main", "_branch"), sort = FALSE)
  if (nrow(merged) == 0L) stop("No overlapping sample IDs between main and branch cluster labels.")

  k_main <- length(unique(merged$cluster_main))
  k_branch <- length(unique(merged$cluster_branch))

  out <- list(
    merged = merged,
    summary = NULL,
    permutation_check = NULL,
    switchers = data.frame()
  )

  # ARI/NMI/VI do not require relabelling, but exact-match and switchers do.
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

  # Only do permutation alignment when k matches and is 2 or 3.
  if (k_main == k_branch && k_main %in% c(2L, 3L)) {
    labs_main <- as.integer(factor(merged$cluster_main, levels = sort(unique(merged$cluster_main))))
    labs_branch <- as.integer(factor(merged$cluster_branch, levels = sort(unique(merged$cluster_branch))))

    perms <- if (k_main == 2L) {
      list(c(1L, 2L), c(2L, 1L))
    } else {
      list(
        c(1L, 2L, 3L), c(1L, 3L, 2L), c(2L, 1L, 3L),
        c(2L, 3L, 1L), c(3L, 1L, 2L), c(3L, 2L, 1L)
      )
    }

    perm_tbl <- do.call(rbind, lapply(seq_along(perms), function(i) {
      perm <- perms[[i]]
      aligned <- perm[labs_branch]
      data.frame(
        permutation_id = i,
        mapping = paste(seq_len(k_main), "->", perm, collapse = "; "),
        exact = mean(labs_main == aligned),
        stringsAsFactors = FALSE
      )
    }))
    best_idx <- which.max(perm_tbl$exact)
    best_perm <- perms[[best_idx]]
    aligned_branch <- best_perm[labs_branch]

    out$summary <- data.frame(
      metric = c("k_main", "k_branch", "ari", "nmi", "vi", "exact_direct", "exact_best", "moved_n", "moved_prop"),
      value = c(
        k_main, k_branch,
        ari(labs_main, aligned_branch),
        nmi(labs_main, aligned_branch),
        vi(labs_main, aligned_branch),
        mean(labs_main == labs_branch),
        mean(labs_main == aligned_branch),
        sum(labs_main != aligned_branch),
        mean(labs_main != aligned_branch)
      ),
      stringsAsFactors = FALSE
    )

    out$permutation_check <- within(perm_tbl, is_best <- seq_len(nrow(perm_tbl)) == best_idx)
    out$switchers <- data.frame(
      sample_id = merged$sample_id[labs_main != aligned_branch],
      main_cluster = labs_main[labs_main != aligned_branch],
      branch_cluster_raw = labs_branch[labs_main != aligned_branch],
      branch_cluster_aligned = aligned_branch[labs_main != aligned_branch],
      stringsAsFactors = FALSE
    )

    merged$branch_cluster_aligned <- aligned_branch
    merged$moved <- labs_main != aligned_branch
  } else {
    out$summary <- data.frame(
      metric = c("k_main", "k_branch", "ari", "nmi", "vi", "exact_direct", "exact_best", "moved_n", "moved_prop"),
      value = c(
        k_main, k_branch,
        ari(merged$cluster_main, merged$cluster_branch),
        nmi(merged$cluster_main, merged$cluster_branch),
        vi(merged$cluster_main, merged$cluster_branch),
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_
      ),
      stringsAsFactors = FALSE
    )
    out$permutation_check <- data.frame()
    out$switchers <- data.frame()
    merged$branch_cluster_aligned <- NA_integer_
    merged$moved <- NA
  }

  out$merged <- merged
  out
}

network_summary <- compare_networks_upper_triangle(main_fused, branch_fused)
comparison <- align_and_compare_labels(main_labels, branch_labels)

prefix <- "clinical_only_vs_main_four_network"
write.csv(network_summary, file.path(out_dir, paste0(prefix, "_network_similarity.csv")), row.names = FALSE)
write.csv(comparison$summary, file.path(out_dir, paste0(prefix, "_summary.csv")), row.names = FALSE)
write.csv(comparison$merged, file.path(out_dir, paste0(prefix, "_matched_labels.csv")), row.names = FALSE)
write.csv(comparison$switchers, file.path(out_dir, paste0(prefix, "_switchers.csv")), row.names = FALSE)
write.csv(comparison$permutation_check, file.path(out_dir, paste0(prefix, "_permutation_check.csv")), row.names = FALSE)

saveRDS(network_summary, file.path(out_dir, paste0(prefix, "_network_similarity.rds")))
saveRDS(comparison$summary, file.path(out_dir, paste0(prefix, "_summary.rds")))
saveRDS(comparison$merged, file.path(out_dir, paste0(prefix, "_matched_labels.rds")))
saveRDS(comparison$switchers, file.path(out_dir, paste0(prefix, "_switchers.rds")))
saveRDS(comparison$permutation_check, file.path(out_dir, paste0(prefix, "_permutation_check.rds")))

message("Clinical-only comparison complete.")
message("Saved outputs to: ", out_dir)
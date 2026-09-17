# 06_posthoc_analysis/07_threshold_scans.R
# Re-runnable from a fresh R session.
#
# - scan simple clinical thresholds against the primary SNF cluster solution
# - ask whether a single clinical threshold (or two thresholds) can approximate
#   the current primary clusters
# - report ARI, NMI, VI, and exact agreement

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
# 1. Utilities
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
    if (anyDuplicated(out$sample_id)) stop(label, " contains duplicated sample IDs.", call. = FALSE)
    return(out)
  }

  if ((is.vector(obj) || is.factor(obj)) && !is.null(names(obj))) {
    out <- data.frame(
      sample_id = clean_sample_ids(names(obj)),
      cluster = parse_cluster_id(obj),
      stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$sample_id), , drop = FALSE]
    if (anyDuplicated(out$sample_id)) stop(label, " contains duplicated sample IDs.", call. = FALSE)
    return(out)
  }

  if (!is.null(fallback_sample_ids) && length(fallback_sample_ids) == length(obj)) {
    out <- data.frame(
      sample_id = clean_sample_ids(fallback_sample_ids),
      cluster = parse_cluster_id(obj),
      stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$sample_id), , drop = FALSE]
    if (anyDuplicated(out$sample_id)) stop(label, " contains duplicated sample IDs.", call. = FALSE)
    return(out)
  }

  stop("Could not standardise ", label, ".", call. = FALSE)
}

make_midpoints <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  ux <- sort(unique(x))
  if (length(ux) < 2) return(numeric(0))
  (ux[-1] + ux[-length(ux)]) / 2
}

choose_threshold_grid <- function(x, max_thresholds = 200L) {
  mids <- make_midpoints(x)
  if (length(mids) <= max_thresholds) {
    return(list(
      thresholds = mids,
      exhaustive = TRUE,
      n_full = length(mids),
      n_used = length(mids)
    ))
  }

  idx <- unique(round(seq(1, length(mids), length.out = max_thresholds)))
  idx <- idx[idx >= 1 & idx <= length(mids)]
  candidate_thresholds <- unique(mids[idx])

  if (length(candidate_thresholds) < 2) {
    stop("Unable to construct a usable reduced threshold grid.", call. = FALSE)
  }

  list(
    thresholds = candidate_thresholds,
    exhaustive = FALSE,
    n_full = length(mids),
    n_used = length(candidate_thresholds)
  )
}

make_binary_split <- function(x, threshold, orientation) {
  x <- suppressWarnings(as.numeric(x))
  if (orientation == "x <= threshold -> cluster 1") {
    ifelse(x <= threshold, 1L, 2L)
  } else if (orientation == "x <= threshold -> cluster 2") {
    ifelse(x <= threshold, 2L, 1L)
  } else {
    stop("Unknown orientation: ", orientation, call. = FALSE)
  }
}

make_three_way_split <- function(x, threshold_1, threshold_2) {
  x <- suppressWarnings(as.numeric(x))
  pred <- rep(NA_integer_, length(x))
  pred[is.finite(x) & x <= threshold_1] <- 1L
  pred[is.finite(x) & x > threshold_1 & x <= threshold_2] <- 2L
  pred[is.finite(x) & x > threshold_2] <- 3L
  pred
}

permute_labels_123 <- function(pred, perm) {
  pred <- as.integer(pred)
  out <- pred
  out[pred == 1L] <- perm[1]
  out[pred == 2L] <- perm[2]
  out[pred == 3L] <- perm[3]
  out
}

evaluate_binary_split <- function(truth, pred) {
  truth <- as.integer(truth)
  pred <- as.integer(pred)
  exact_direct <- mean(truth == pred)
  exact_flipped <- mean(truth == ifelse(pred == 1L, 2L, 1L))
  data.frame(
    ari = adjusted_rand_index(truth, pred),
    nmi = normalized_mutual_information(truth, pred),
    vi = variation_of_information(truth, pred),
    exact_direct = exact_direct,
    exact_flipped = exact_flipped,
    exact_best = max(exact_direct, exact_flipped),
    stringsAsFactors = FALSE
  )
}

evaluate_three_way_split <- function(truth, pred) {
  truth <- as.integer(truth)
  pred <- as.integer(pred)

  exact_direct <- mean(truth == pred)
  perms <- list(
    c(1L, 2L, 3L), c(1L, 3L, 2L), c(2L, 1L, 3L),
    c(2L, 3L, 1L), c(3L, 1L, 2L), c(3L, 2L, 1L)
  )

  best_exact <- -Inf
  best_perm <- NA_character_
  for (perm in perms) {
    perm_pred <- permute_labels_123(pred, perm)
    ex <- mean(truth == perm_pred)
    if (is.finite(ex) && ex > best_exact) {
      best_exact <- ex
      best_perm <- paste0(perm, collapse = "-")
    }
  }

  data.frame(
    ari = adjusted_rand_index(truth, pred),
    nmi = normalized_mutual_information(truth, pred),
    vi = variation_of_information(truth, pred),
    exact_direct = exact_direct,
    exact_best = best_exact,
    best_perm = best_perm,
    stringsAsFactors = FALSE
  )
}

make_contingency_table <- function(truth, pred, truth_levels = NULL, pred_levels = NULL) {
  truth <- as.integer(truth)
  pred <- as.integer(pred)
  keep <- is.finite(truth) & is.finite(pred)
  truth <- truth[keep]
  pred <- pred[keep]

  if (is.null(truth_levels)) truth_levels <- sort(unique(truth))
  if (is.null(pred_levels)) pred_levels <- sort(unique(pred))

  tab <- table(factor(truth, levels = truth_levels), factor(pred, levels = pred_levels))
  tab_df <- as.data.frame(tab)
  names(tab_df) <- c("truth_cluster", "threshold_cluster", "count")
  tab_df$truth_cluster <- as.integer(as.character(tab_df$truth_cluster))
  tab_df$threshold_cluster <- as.integer(as.character(tab_df$threshold_cluster))
  tab_df$proportion_of_all <- tab_df$count / sum(tab_df$count)

  row_prop <- prop.table(tab, margin = 1)
  row_prop_df <- as.data.frame(row_prop)
  names(row_prop_df) <- c("truth_cluster", "threshold_cluster", "row_proportion")
  row_prop_df$truth_cluster <- as.integer(as.character(row_prop_df$truth_cluster))
  row_prop_df$threshold_cluster <- as.integer(as.character(row_prop_df$threshold_cluster))

  col_prop <- prop.table(tab, margin = 2)
  col_prop_df <- as.data.frame(col_prop)
  names(col_prop_df) <- c("truth_cluster", "threshold_cluster", "col_proportion")
  col_prop_df$truth_cluster <- as.integer(as.character(col_prop_df$truth_cluster))
  col_prop_df$threshold_cluster <- as.integer(as.character(col_prop_df$threshold_cluster))

  list(
    table = tab,
    table_df = tab_df,
    row_prop_df = row_prop_df,
    col_prop_df = col_prop_df
  )
}

plot_binary_threshold_scan <- function(scan_df, best_row, out_file, var_name) {
  plot_df <- aggregate(ari ~ threshold_1 + orientation, data = scan_df, FUN = max, na.rm = TRUE)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = threshold_1, y = ari, color = orientation)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::geom_point(
      data = best_row,
      ggplot2::aes(x = threshold_1, y = ari),
      inherit.aes = FALSE,
      size = 2.4,
      shape = 21,
      fill = "white",
      color = "black"
    ) +
    ggplot2::geom_vline(
      data = best_row,
      ggplot2::aes(xintercept = threshold_1),
      inherit.aes = FALSE,
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = paste0(var_name, " threshold scan"),
      x = "Threshold",
      y = "ARI",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  ggplot2::ggsave(out_file, p, width = 9, height = 6.5, dpi = 300)
  invisible(p)
}

plot_three_way_threshold_scan <- function(scan_df, best_row, out_file, var_name) {
  p <- ggplot2::ggplot(scan_df, ggplot2::aes(x = threshold_1, y = threshold_2, fill = ari)) +
    ggplot2::geom_tile() +
    ggplot2::geom_point(
      data = best_row,
      ggplot2::aes(x = threshold_1, y = threshold_2),
      inherit.aes = FALSE,
      shape = 21,
      size = 3,
      stroke = 0.7,
      fill = NA,
      color = "black"
    ) +
    ggplot2::labs(
      title = paste0(var_name, " two-threshold scan"),
      x = "Lower threshold",
      y = "Upper threshold",
      fill = "ARI"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  ggplot2::ggsave(out_file, p, width = 9, height = 7, dpi = 300)
  invisible(p)
}

# -----------------------------------------------------------------------------
# 2. Load analysis cohort and primary cluster labels
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

truth <- as.integer(primary_df$cluster)
k_value <- length(unique(truth))
if (!k_value %in% c(2L, 3L)) {
  stop(
    "This threshold-scan script is set up for a primary solution with 2 or 3 clusters. ",
    "The current primary solution has k = ", k_value, ".",
    call. = FALSE
  )
}

# Variables to scan: keep this conservative and explicit.
vars_to_scan <- c("UPDRS_III_num", "MOCA_total_num", "LEDD_total_num", "AGE_num", "disease_duration")
vars_to_scan <- intersect(vars_to_scan, names(pd_model))
if (length(vars_to_scan) == 0) {
  stop("None of the requested clinical variables are present in pd_model.", call. = FALSE)
}

cat("Primary cluster count: ", k_value, "\n", sep = "")
cat("Variables to scan:\n")
print(vars_to_scan)

# -----------------------------------------------------------------------------
# 3. Output directories
# -----------------------------------------------------------------------------
out_root <- get_posthoc_output_dir("threshold_scans")
plot_dir <- file.path(out_root, "plots")
table_dir <- file.path(out_root, "tables")
diag_dir <- file.path(out_root, "diagnostics")
safe_dir(plot_dir)
safe_dir(table_dir)
safe_dir(diag_dir)

# -----------------------------------------------------------------------------
# 4. Threshold scans
# -----------------------------------------------------------------------------
max_k3_threshold_candidates <- 200L

scan_threshold_variable <- function(x, truth, var_name, k_value, max_k3_threshold_candidates = 200L) {
  x <- suppressWarnings(as.numeric(x))
  keep <- is.finite(x) & !is.na(truth)
  x <- x[keep]
  truth <- truth[keep]

  if (length(unique(truth)) < k_value) {
    stop(
      "Truth labels do not contain at least ", k_value,
      " clusters for variable: ", var_name,
      " (k = ", k_value, ")",
      call. = FALSE
    )
  }

  if (k_value == 2L) {
    thresholds <- make_midpoints(x)
    if (length(thresholds) < 1L) {
      stop("Variable has fewer than 2 unique values: ", var_name, call. = FALSE)
    }

    scan_rows <- vector("list", length(thresholds) * 2L)
    row_idx <- 1L

    for (thr in thresholds) {
      pred_a <- ifelse(x <= thr, 1L, 2L)
      met_a <- evaluate_binary_split(truth, pred_a)
      scan_rows[[row_idx]] <- data.frame(
        k = 2L,
        variable = var_name,
        threshold_1 = thr,
        threshold_2 = NA_real_,
        orientation = "x <= threshold -> cluster 1",
        n_used = length(x),
        scan_method = "exhaustive",
        n_candidate_thresholds_full = length(thresholds),
        n_candidate_thresholds_used = length(thresholds),
        ari = met_a$ari,
        nmi = met_a$nmi,
        vi = met_a$vi,
        exact_direct = met_a$exact_direct,
        exact_flipped = met_a$exact_flipped,
        exact_best = met_a$exact_best,
        best_perm = NA_character_,
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1L

      pred_b <- ifelse(x <= thr, 2L, 1L)
      met_b <- evaluate_binary_split(truth, pred_b)
      scan_rows[[row_idx]] <- data.frame(
        k = 2L,
        variable = var_name,
        threshold_1 = thr,
        threshold_2 = NA_real_,
        orientation = "x <= threshold -> cluster 2",
        n_used = length(x),
        scan_method = "exhaustive",
        n_candidate_thresholds_full = length(thresholds),
        n_candidate_thresholds_used = length(thresholds),
        ari = met_b$ari,
        nmi = met_b$nmi,
        vi = met_b$vi,
        exact_direct = met_b$exact_direct,
        exact_flipped = met_b$exact_flipped,
        exact_best = met_b$exact_best,
        best_perm = NA_character_,
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1L
    }

    scan_df <- do.call(rbind, scan_rows)
    ord <- order(-scan_df$ari, -scan_df$exact_best, -scan_df$nmi, scan_df$vi)
    best_row <- scan_df[ord[1], , drop = FALSE]

    list(scan = scan_df, best = best_row)
  } else if (k_value == 3L) {
    candidate_info <- choose_threshold_grid(x, max_thresholds = max_k3_threshold_candidates)
    thresholds <- candidate_info$thresholds

    if (length(thresholds) < 2L) {
      stop("Variable has fewer than 2 candidate thresholds: ", var_name, call. = FALSE)
    }

    n_pairs <- choose(length(thresholds), 2)
    scan_rows <- vector("list", n_pairs)
    row_idx <- 1L

    for (i in seq_len(length(thresholds) - 1L)) {
      t1 <- thresholds[i]
      for (j in seq.int(i + 1L, length(thresholds))) {
        t2 <- thresholds[j]
        pred <- make_three_way_split(x, t1, t2)
        met <- evaluate_three_way_split(truth, pred)
        scan_rows[[row_idx]] <- data.frame(
          k = 3L,
          variable = var_name,
          threshold_1 = t1,
          threshold_2 = t2,
          orientation = "low-mid-high",
          n_used = length(x),
          scan_method = if (candidate_info$exhaustive) "exhaustive" else "thinned",
          n_candidate_thresholds_full = candidate_info$n_full,
          n_candidate_thresholds_used = candidate_info$n_used,
          ari = met$ari,
          nmi = met$nmi,
          vi = met$vi,
          exact_direct = met$exact_direct,
          exact_flipped = NA_real_,
          exact_best = met$exact_best,
          best_perm = met$best_perm,
          stringsAsFactors = FALSE
        )
        row_idx <- row_idx + 1L
      }
    }

    scan_df <- do.call(rbind, scan_rows)
    ord <- order(-scan_df$ari, -scan_df$exact_best, -scan_df$nmi, scan_df$vi)
    best_row <- scan_df[ord[1], , drop = FALSE]

    list(scan = scan_df, best = best_row, candidate_info = candidate_info)
  } else {
    stop("k_value must be 2 or 3.", call. = FALSE)
  }
}

threshold_results <- list()
best_rows <- list()
contingency_results <- list()
contingency_rows <- list()
scan_notes <- list()

for (v in vars_to_scan) {
  message("Scanning thresholds for: ", v)
  res <- scan_threshold_variable(
    x = pd_model[[v]],
    truth = truth,
    var_name = v,
    k_value = k_value,
    max_k3_threshold_candidates = max_k3_threshold_candidates
  )

  file_prefix <- paste0("k", k_value, "_", v)
  scan_df <- res$scan
  best_row <- res$best

  threshold_results[[file_prefix]] <- scan_df
  best_rows[[file_prefix]] <- best_row

  write.csv(scan_df, file.path(table_dir, paste0(file_prefix, "_threshold_scan_full.csv")), row.names = FALSE)
  saveRDS(scan_df, file.path(table_dir, paste0(file_prefix, "_threshold_scan_full.rds")))

  if (k_value == 2L) {
    pred_best <- make_binary_split(pd_model[[v]], best_row$threshold_1, best_row$orientation)
    cont <- make_contingency_table(
      truth = truth,
      pred = pred_best,
      truth_levels = sort(unique(truth)),
      pred_levels = c(1L, 2L)
    )
    contingency_results[[file_prefix]] <- cont
    write.csv(cont$table_df, file.path(table_dir, paste0(file_prefix, "_contingency_counts.csv")), row.names = FALSE)
    write.csv(cont$row_prop_df, file.path(table_dir, paste0(file_prefix, "_contingency_row_proportions.csv")), row.names = FALSE)
    write.csv(cont$col_prop_df, file.path(table_dir, paste0(file_prefix, "_contingency_column_proportions.csv")), row.names = FALSE)

    contingency_rows[[file_prefix]] <- data.frame(
      k = 2L,
      variable = v,
      threshold_1 = best_row$threshold_1,
      threshold_2 = NA_real_,
      orientation = best_row$orientation,
      scan_method = "exhaustive",
      n_candidate_thresholds_full = length(make_midpoints(pd_model[[v]])),
      n_candidate_thresholds_used = length(make_midpoints(pd_model[[v]])),
      ari = best_row$ari,
      nmi = best_row$nmi,
      vi = best_row$vi,
      exact_direct = best_row$exact_direct,
      exact_flipped = best_row$exact_flipped,
      exact_best = best_row$exact_best,
      best_perm = NA_character_,
      total_n = sum(cont$table),
      stringsAsFactors = FALSE
    )

    plot_binary_threshold_scan(
      scan_df = scan_df,
      best_row = best_row,
      out_file = file.path(plot_dir, paste0(file_prefix, "_threshold_scan.png")),
      var_name = v
    )

    scan_notes[[file_prefix]] <- data.frame(
      k = 2L,
      variable = v,
      scan_method = "exhaustive",
      n_candidate_thresholds_full = length(make_midpoints(pd_model[[v]])),
      n_candidate_thresholds_used = length(make_midpoints(pd_model[[v]])),
      stringsAsFactors = FALSE
    )
  } else {
    pred_best <- make_three_way_split(pd_model[[v]], best_row$threshold_1, best_row$threshold_2)
    if (!is.na(best_row$best_perm)) {
      perm_vals <- as.integer(strsplit(best_row$best_perm, "-")[[1]])
      pred_best <- permute_labels_123(pred_best, perm_vals)
    }

    cont <- make_contingency_table(
      truth = truth,
      pred = pred_best,
      truth_levels = sort(unique(truth)),
      pred_levels = c(1L, 2L, 3L)
    )
    contingency_results[[file_prefix]] <- cont
    write.csv(cont$table_df, file.path(table_dir, paste0(file_prefix, "_contingency_counts.csv")), row.names = FALSE)
    write.csv(cont$row_prop_df, file.path(table_dir, paste0(file_prefix, "_contingency_row_proportions.csv")), row.names = FALSE)
    write.csv(cont$col_prop_df, file.path(table_dir, paste0(file_prefix, "_contingency_column_proportions.csv")), row.names = FALSE)

    contingency_rows[[file_prefix]] <- data.frame(
      k = 3L,
      variable = v,
      threshold_1 = best_row$threshold_1,
      threshold_2 = best_row$threshold_2,
      orientation = best_row$orientation,
      scan_method = best_row$scan_method,
      n_candidate_thresholds_full = res$candidate_info$n_full,
      n_candidate_thresholds_used = res$candidate_info$n_used,
      ari = best_row$ari,
      nmi = best_row$nmi,
      vi = best_row$vi,
      exact_direct = best_row$exact_direct,
      exact_flipped = NA_real_,
      exact_best = best_row$exact_best,
      best_perm = best_row$best_perm,
      total_n = sum(cont$table),
      stringsAsFactors = FALSE
    )

    plot_three_way_threshold_scan(
      scan_df = scan_df,
      best_row = best_row,
      out_file = file.path(plot_dir, paste0(file_prefix, "_threshold_scan.png")),
      var_name = v
    )

    scan_notes[[file_prefix]] <- data.frame(
      k = 3L,
      variable = v,
      scan_method = best_row$scan_method,
      n_candidate_thresholds_full = res$candidate_info$n_full,
      n_candidate_thresholds_used = res$candidate_info$n_used,
      exhaustive = res$candidate_info$exhaustive,
      stringsAsFactors = FALSE
    )

    if (!res$candidate_info$exhaustive) {
      message(
        "k = 3 scan for ", v,
        " used a reduced candidate grid: ",
        res$candidate_info$n_used,
        " of ",
        res$candidate_info$n_full,
        " midpoint thresholds."
      )
    }
  }
}

# -----------------------------------------------------------------------------
# 5. Combined tables
# -----------------------------------------------------------------------------
full_scan_table <- do.call(rbind, threshold_results)
best_scan_table <- do.call(rbind, best_rows)
contingency_summary_table <- do.call(rbind, contingency_rows)
scan_notes_table <- do.call(rbind, scan_notes)

best_scan_table <- best_scan_table[order(best_scan_table$k, -best_scan_table$ari, -best_scan_table$exact_best, -best_scan_table$nmi, best_scan_table$vi), ]
contingency_summary_table <- contingency_summary_table[order(contingency_summary_table$k, -contingency_summary_table$ari, -contingency_summary_table$nmi, -contingency_summary_table$exact_best), ]
scan_notes_table <- scan_notes_table[order(scan_notes_table$k, scan_notes_table$variable), ]

rownames(full_scan_table) <- NULL
rownames(best_scan_table) <- NULL
rownames(contingency_summary_table) <- NULL
rownames(scan_notes_table) <- NULL

write.csv(full_scan_table, file.path(out_root, "threshold_scan_all_variables_full.csv"), row.names = FALSE)
write.csv(best_scan_table, file.path(out_root, "threshold_scan_best_per_variable.csv"), row.names = FALSE)
write.csv(contingency_summary_table, file.path(out_root, "threshold_scan_contingency_summary.csv"), row.names = FALSE)
write.csv(scan_notes_table, file.path(out_root, "threshold_scan_notes.csv"), row.names = FALSE)

saveRDS(
  list(
    k_value = k_value,
    primary_model = primary$branch_name %||% "primary",
    vars_to_scan = vars_to_scan,
    full_scan_table = full_scan_table,
    best_scan_table = best_scan_table,
    contingency_summary_table = contingency_summary_table,
    scan_notes_table = scan_notes_table,
    contingency_results = contingency_results
  ),
  file.path(out_root, "threshold_scan_results.rds")
)

saveRDS(contingency_results, file.path(out_root, "threshold_scan_contingency_results.rds"))

# -----------------------------------------------------------------------------
# 6. Console summary
# -----------------------------------------------------------------------------
cat("\nBest threshold scan rows\n")
print(best_scan_table[, c(
  "k", "variable", "threshold_1", "threshold_2", "orientation",
  "scan_method", "n_candidate_thresholds_full", "n_candidate_thresholds_used",
  "ari", "nmi", "vi", "exact_best"
)])

cat("\nContingency summary for best thresholds\n")
print(contingency_summary_table[, c(
  "k", "variable", "threshold_1", "threshold_2", "orientation",
  "scan_method", "n_candidate_thresholds_full", "n_candidate_thresholds_used",
  "ari", "nmi", "vi", "exact_best", "total_n"
)])

writeLines(capture.output(sessionInfo()), file.path(out_root, "session_info.txt"))

cat("\nFinished threshold scan.\n")
cat("Outputs saved to: ", out_root, "\n", sep = "")
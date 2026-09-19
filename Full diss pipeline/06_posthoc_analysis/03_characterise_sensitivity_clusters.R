# 06_posthoc_analysis/03_characterise_sensitivity_clusters.R

# This script characterises the clusters produced by each sensitivity analysis
# using clinical variables. It summarises cluster sizes and clinical
# characteristics, tests for differences between clusters, and produces
# clinical plots for each sensitivity branch.
# - load sensitivity branch cluster labels
# - merge labels onto the clinical cohort
# - mirror the primary-cluster characterisation logic
# - use raw clinical variables where available, with safe numeric coercion
# - compute cluster sizes, descriptive stats, Wilcoxon/Kruskal tests, and effect sizes
#
# Sensitivity branches:
# - Clinical-only
# - No-LEDD
# - Representation comparison
# - Top-20 combined
# - Top-20 23-network
#
# Before running:
# - Set `project_root` to the local project directory.
# - `analysis_cohort.rds` and the cluster solutions for all sensitivity
#   branches must already exist.
#
# Outputs are saved under `paths$posthoc/sensitivity_clusters/`, with separate
# results for each branch and combined summary tables across branches.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"

source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

# -------------------------------------------------------------------
# Basic checks
# -------------------------------------------------------------------
if (!file.exists(files$analysis_cohort)) {
  stop("Missing analysis cohort file: ", files$analysis_cohort, call. = FALSE)
}

if (!dir.exists(paths$posthoc)) {
  dir.create(paths$posthoc, recursive = TRUE, showWarnings = FALSE)
}

out_root <- file.path(paths$posthoc, "sensitivity_clusters")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------
clean_sample_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[nchar(x) == 0] <- NA_character_
  x
}

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

to_numeric_safe <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  suppressWarnings(as.numeric(x))
}

pick_first_available <- function(candidates, data) {
  found <- candidates[candidates %in% names(data)]
  if (length(found) == 0) return(NA_character_)
  found[1]
}

parse_cluster_id <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.numeric(x) || is.integer(x)) {
    out <- as.integer(x)
    if (anyNA(out)) stop("Cluster labels contain NA values after numeric conversion.", call. = FALSE)
    return(out)
  }

  x <- as.character(x)
  digits <- stringr::str_extract(x, "-?[0-9]+")
  out <- suppressWarnings(as.integer(digits))

  if (all(!is.na(out))) {
    return(out)
  }

  # If labels are non-numeric text, map to stable integers in observed order.
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

safe_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

safe_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::IQR(x, na.rm = TRUE)
}

cohens_d <- function(x, g) {
  x <- to_numeric_safe(x)
  g <- as.factor(g)
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

epsilon_squared_kw <- function(H, k, n) {
  if (!is.finite(H) || !is.finite(k) || !is.finite(n) || n <= k) return(NA_real_)
  (H - k + 1) / (n - k)
}

cramers_v <- function(tab) {
  tab <- as.matrix(tab)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  n <- sum(tab)
  if (!is.finite(chi$statistic) || !is.finite(n) || n <= 0) return(NA_real_)
  sqrt(as.numeric(chi$statistic) / (n * min(nrow(tab) - 1, ncol(tab) - 1)))
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
        stop("Could not identify sample_id and cluster columns in labels data frame.", call. = FALSE)
      }
    } else {
      out <- data.frame(
        sample_id = as.character(obj[[sid_col[1]]]),
        cluster = obj[[cl_col[1]]],
        stringsAsFactors = FALSE
      )
    }
  } else if (is.list(obj)) {
    candidate_fields <- c(
      "best_k_cluster_labels",
      "cluster_labels",
      "labels",
      "cluster_df",
      "best_labels"
    )

    for (fld in candidate_fields) {
      if (!is.null(obj[[fld]])) {
        return(extract_cluster_labels(obj[[fld]]))
      }
    }

    if (!is.null(obj$sample_ids) && !is.null(obj$cluster_labels)) {
      return(extract_cluster_labels(obj$cluster_labels))
    }

    stop("Could not interpret the sensitivity cluster object.", call. = FALSE)
  } else if (is.vector(obj) && !is.null(names(obj))) {
    out <- data.frame(
      sample_id = as.character(names(obj)),
      cluster = as.vector(obj),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Could not interpret the sensitivity cluster labels object.", call. = FALSE)
  }

  out$sample_id <- clean_sample_id(out$sample_id)
  out <- out[!is.na(out$sample_id), , drop = FALSE]
  out$cluster <- parse_cluster_id(out$cluster)

  if (anyDuplicated(out$sample_id)) {
    stop("Cluster labels contain duplicated sample IDs.", call. = FALSE)
  }

  out[order(out$sample_id), , drop = FALSE]
}

save_cluster_barplot <- function(cluster_sizes, out_file, title_text, xlab_text = "Cluster", ylab_text = "Patients") {
  png(out_file, width = 1800, height = 1200, res = 250)
  on.exit(dev.off(), add = TRUE)
  barplot(
    height = cluster_sizes$n,
    names.arg = cluster_sizes$cluster,
    xlab = xlab_text,
    ylab = ylab_text,
    main = title_text
  )
}

save_boxplot <- function(df, var, cluster_col, out_file, title_text) {
  x <- df[[var]]
  g <- as.factor(df[[cluster_col]])

  png(out_file, width = 2000, height = 1500, res = 250)
  on.exit(dev.off(), add = TRUE)
  boxplot(
    x ~ g,
    xlab = "Cluster",
    ylab = var,
    main = title_text
  )
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

# -------------------------------------------------------------------
# Sensitivity branch definitions
# -------------------------------------------------------------------
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

# -------------------------------------------------------------------
# Load source clinical cohort
# -------------------------------------------------------------------
analysis_cohort <- readRDS(files$analysis_cohort)
pd_model <- as.data.frame(analysis_cohort$pd_model, check.names = FALSE)

if (is.null(rownames(pd_model)) && "Anonymised_sampleID" %in% names(pd_model)) {
  rownames(pd_model) <- as.character(pd_model$Anonymised_sampleID)
}

sample_id_col <- pick_first_available(c("Anonymised_sampleID", "sample_id"), pd_model)
if (is.na(sample_id_col)) {
  stop("Could not find a sample ID column in pd_model.", call. = FALSE)
}

pd_model$sample_id <- clean_sample_id(pd_model[[sample_id_col]])
if (anyDuplicated(pd_model$sample_id)) {
  stop("pd_model contains duplicated sample IDs.", call. = FALSE)
}

# -------------------------------------------------------------------
# Branch processing functions
# -------------------------------------------------------------------
make_continuous_summary <- function(df, continuous_vars, cluster_col) {
  out_list <- list()

  for (v in continuous_vars) {
    for (cl in sort(unique(df[[cluster_col]]))) {
      x <- df[df[[cluster_col]] == cl, v]
      out_list[[length(out_list) + 1L]] <- data.frame(
        variable = v,
        cluster = as.integer(cl),
        n = sum(is.finite(x)),
        mean = safe_mean(x),
        sd = safe_sd(x),
        median = safe_median(x),
        iqr = safe_iqr(x),
        min = suppressWarnings(min(x, na.rm = TRUE)),
        q25 = suppressWarnings(stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE)),
        q75 = suppressWarnings(stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)),
        max = suppressWarnings(max(x, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, out_list)
}

run_continuous_tests <- function(df, continuous_vars, cluster_col) {
  clusters <- as.factor(df[[cluster_col]])
  k <- nlevels(clusters)
  out_list <- list()

  for (v in continuous_vars) {
    x <- df[[v]]

    if (k == 2L) {
      levs <- levels(clusters)
      x1 <- x[clusters == levs[1]]
      x2 <- x[clusters == levs[2]]

      wt <- tryCatch(
        stats::wilcox.test(x ~ clusters, exact = FALSE, conf.int = TRUE),
        error = function(e) NULL
      )

      out_list[[length(out_list) + 1L]] <- data.frame(
        variable = v,
        test_name = "wilcox.test",
        n_clusters = 2L,
        statistic = if (is.null(wt)) NA_real_ else unname(wt$statistic),
        p_value = if (is.null(wt)) NA_real_ else wt$p.value,
        conf_low = if (is.null(wt) || is.null(wt$conf.int)) NA_real_ else wt$conf.int[1],
        conf_high = if (is.null(wt) || is.null(wt$conf.int)) NA_real_ else wt$conf.int[2],
        effect_size_type = "Cohen's d",
        effect_size = cohens_d(x, clusters),
        group1 = levs[1],
        group2 = levs[2],
        group1_mean = mean(x1, na.rm = TRUE),
        group2_mean = mean(x2, na.rm = TRUE),
        group1_sd = stats::sd(x1, na.rm = TRUE),
        group2_sd = stats::sd(x2, na.rm = TRUE),
        group1_n = sum(is.finite(x1)),
        group2_n = sum(is.finite(x2)),
        stringsAsFactors = FALSE
      )
    } else {
      kt <- tryCatch(stats::kruskal.test(x ~ clusters), error = function(e) NULL)
      H <- if (is.null(kt)) NA_real_ else unname(kt$statistic)

      out_list[[length(out_list) + 1L]] <- data.frame(
        variable = v,
        test_name = "kruskal.test",
        n_clusters = k,
        statistic = H,
        p_value = if (is.null(kt)) NA_real_ else kt$p.value,
        conf_low = NA_real_,
        conf_high = NA_real_,
        effect_size_type = "epsilon-squared",
        effect_size = epsilon_squared_kw(H, k = k, n = sum(is.finite(x))),
        group1 = NA_character_,
        group2 = NA_character_,
        group1_mean = NA_real_,
        group2_mean = NA_real_,
        group1_sd = NA_real_,
        group2_sd = NA_real_,
        group1_n = NA_integer_,
        group2_n = NA_integer_,
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, out_list)
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  out
}

run_categorical_tests <- function(df, categorical_vars, cluster_col) {
  clusters <- as.factor(df[[cluster_col]])
  out_list <- list()

  for (v in categorical_vars) {
    x <- as.factor(df[[v]])
    tab <- table(clusters, x)

    if (nrow(tab) < 2 || ncol(tab) < 2) {
      next
    }

    chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
    effect <- cramers_v(tab)

    out_list[[length(out_list) + 1L]] <- data.frame(
      variable = v,
      test_name = "chisq.test",
      statistic = unname(chi$statistic),
      p_value = chi$p.value,
      effect_size_type = "Cramer's V",
      effect_size = effect,
      n_clusters = nlevels(clusters),
      n_levels = nlevels(x),
      stringsAsFactors = FALSE
    )
  }

  if (length(out_list) == 0) {
    return(data.frame())
  }

  out <- do.call(rbind, out_list)
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  out
}

make_gender_table <- function(df, cluster_col, gender_col) {
  tab <- table(df[[cluster_col]], df[[gender_col]])
  tab_df <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(tab_df) <- c("cluster", "gender", "count")
  tab_df$prop_within_cluster <- tab_df$count / ave(tab_df$count, tab_df$cluster, FUN = sum)
  tab_df
}

save_gender_plot <- function(df, cluster_col, gender_col, out_file, title_text) {
  gender_tab <- table(df[[cluster_col]], df[[gender_col]])
  png(out_file, width = 2000, height = 1500, res = 250)
  on.exit(dev.off(), add = TRUE)
  barplot(
    gender_tab,
    beside = TRUE,
    legend.text = TRUE,
    xlab = "Cluster",
    ylab = "Count",
    main = title_text
  )
}

process_branch <- function(branch_name, branch_dir, result_candidates) {
  out_dir <- file.path(out_root, branch_name)
  plot_dir <- file.path(out_dir, "plots")
  table_dir <- file.path(out_dir, "tables")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  result_file <- find_first_existing(result_candidates)
  if (is.na(result_file) || !nzchar(result_file)) {
    stop(
      "Could not find a sensitivity cluster result file for branch '",
      branch_name,
      "'. Looked in:\n",
      paste(result_candidates, collapse = "\n"),
      call. = FALSE
    )
  }

  branch_obj <- readRDS(result_file)
  branch_labels <- extract_cluster_labels(branch_obj)

  branch_labels$sample_id <- clean_sample_id(branch_labels$sample_id)

  branch_model_name <- if (is.list(branch_obj) && !is.null(branch_obj$model_name)) {
    as.character(branch_obj$model_name)
  } else {
    branch_name
  }

  if (!all(branch_labels$sample_id %in% pd_model$sample_id)) {
    missing_ids <- setdiff(branch_labels$sample_id, pd_model$sample_id)
    stop(
      "Some branch sample IDs are missing from the clinical cohort for branch '",
      branch_name,
      "'. First missing ID: ",
      missing_ids[1],
      call. = FALSE
    )
  }

  pd_branch <- pd_model[match(branch_labels$sample_id, pd_model$sample_id), , drop = FALSE]
  if (!identical(pd_branch$sample_id, branch_labels$sample_id)) {
    stop("Failed to align branch clinical data to cluster labels for branch '", branch_name, "'.", call. = FALSE)
  }

  pd_branch$cluster <- factor(branch_labels$cluster, levels = sort(unique(branch_labels$cluster)))
  pd_branch$cluster_id <- as.integer(pd_branch$cluster)
  pd_branch$branch <- branch_name

  # ---------------------------------------------------------------
  # Select continuous and categorical variables
  # ---------------------------------------------------------------
  raw_continuous_vars <- c(
    pick_first_available(c("AGE", "AGE_num"), pd_branch),
    pick_first_available(c("UPDRS_III", "UPDRS_III_num"), pd_branch),
    pick_first_available(c("MOCA_total", "MOCA_total_num"), pd_branch),
    pick_first_available(c("LEDD_total", "LEDD_total_num"), pd_branch),
    pick_first_available(c("disease_duration"), pd_branch)
  )
  raw_continuous_vars <- raw_continuous_vars[!is.na(raw_continuous_vars)]

  categorical_vars <- c(
    pick_first_available(c("GENDER"), pd_branch)
  )
  categorical_vars <- categorical_vars[!is.na(categorical_vars)]

  continuous_data <- pd_branch[, c("sample_id", raw_continuous_vars, "cluster", "cluster_id"), drop = FALSE]
  for (v in raw_continuous_vars) {
    continuous_data[[v]] <- to_numeric_safe(continuous_data[[v]])
  }

  continuous_vars <- raw_continuous_vars
  keep_continuous <- vapply(continuous_data[continuous_vars], function(z) any(is.finite(z)), logical(1))
  continuous_vars <- continuous_vars[keep_continuous]
  continuous_data <- continuous_data[, c("sample_id", continuous_vars, "cluster", "cluster_id"), drop = FALSE]

  if (length(continuous_vars) == 0) {
    stop("No usable continuous clinical variables were found for branch '", branch_name, "'.", call. = FALSE)
  }

  if (length(categorical_vars) > 0) {
    for (v in categorical_vars) {
      continuous_data[[v]] <- pd_branch[[v]]
    }
  }

  # ---------------------------------------------------------------
  # Cluster sizes
  # ---------------------------------------------------------------
  cluster_size_table <- as.data.frame(table(pd_branch$cluster), stringsAsFactors = FALSE)
  names(cluster_size_table) <- c("cluster", "n")
  cluster_size_table$prop <- cluster_size_table$n / sum(cluster_size_table$n)

  write.csv(cluster_size_table, file.path(table_dir, paste0(branch_name, "_cluster_sizes.csv")), row.names = FALSE)
  saveRDS(cluster_size_table, file.path(table_dir, paste0(branch_name, "_cluster_sizes.rds")))

  save_cluster_barplot(
    cluster_sizes = cluster_size_table,
    out_file = file.path(plot_dir, paste0(branch_name, "_cluster_sizes.png")),
    title_text = paste0(branch_name, " cluster sizes")
  )

  # ---------------------------------------------------------------
  # Continuous summaries and tests
  # ---------------------------------------------------------------
  cluster_means <- aggregate(
    continuous_data[, continuous_vars, drop = FALSE],
    by = list(cluster = continuous_data$cluster),
    FUN = function(z) mean(z, na.rm = TRUE)
  )
  write.csv(cluster_means, file.path(table_dir, paste0(branch_name, "_cluster_means.csv")), row.names = FALSE)
  saveRDS(cluster_means, file.path(table_dir, paste0(branch_name, "_cluster_means.rds")))

  continuous_summary <- make_continuous_summary(
    df = continuous_data,
    continuous_vars = continuous_vars,
    cluster_col = "cluster"
  )
  write.csv(
    continuous_summary,
    file.path(table_dir, paste0(branch_name, "_continuous_summary.csv")),
    row.names = FALSE
  )
  saveRDS(
    continuous_summary,
    file.path(table_dir, paste0(branch_name, "_continuous_summary.rds"))
  )

  continuous_tests_df <- run_continuous_tests(
    df = continuous_data,
    continuous_vars = continuous_vars,
    cluster_col = "cluster"
  )
  write.csv(
    continuous_tests_df,
    file.path(table_dir, paste0(branch_name, "_continuous_tests.csv")),
    row.names = FALSE
  )
  saveRDS(
    continuous_tests_df,
    file.path(table_dir, paste0(branch_name, "_continuous_tests.rds"))
  )

  # Plot per continuous variable
  for (v in continuous_vars) {
    safe_v <- safe_filename(v)
    save_boxplot(
      df = continuous_data,
      var = v,
      cluster_col = "cluster",
      out_file = file.path(plot_dir, paste0(branch_name, "_", safe_v, "_by_cluster.png")),
      title_text = paste0(v, " by ", branch_name, " cluster")
    )
  }

  # ---------------------------------------------------------------
  # Categorical summaries and tests
  # ---------------------------------------------------------------
  categorical_summary <- data.frame()
  categorical_tests_df <- data.frame()

  if (length(categorical_vars) > 0) {
    cat_sum_list <- list()

    for (v in categorical_vars) {
      x <- as.factor(continuous_data[[v]])
      tab <- table(continuous_data$cluster, x)
      tab_df <- as.data.frame(tab, stringsAsFactors = FALSE)
      names(tab_df) <- c("cluster", "category", "count")
      tab_df$prop_within_cluster <- tab_df$count / ave(tab_df$count, tab_df$cluster, FUN = sum)
      tab_df$variable <- v
      cat_sum_list[[length(cat_sum_list) + 1L]] <- tab_df

      safe_v <- safe_filename(v)
      save_gender_plot(
        df = continuous_data,
        cluster_col = "cluster",
        gender_col = v,
        out_file = file.path(plot_dir, paste0(branch_name, "_", safe_v, "_by_cluster.png")),
        title_text = paste0(v, " distribution by ", branch_name, " cluster")
      )
    }

    categorical_summary <- do.call(rbind, cat_sum_list)
    write.csv(
      categorical_summary,
      file.path(table_dir, paste0(branch_name, "_categorical_summary.csv")),
      row.names = FALSE
    )
    saveRDS(
      categorical_summary,
      file.path(table_dir, paste0(branch_name, "_categorical_summary.rds"))
    )

    categorical_tests_df <- run_categorical_tests(
      df = continuous_data,
      categorical_vars = categorical_vars,
      cluster_col = "cluster"
    )
    write.csv(
      categorical_tests_df,
      file.path(table_dir, paste0(branch_name, "_categorical_tests.csv")),
      row.names = FALSE
    )
    saveRDS(
      categorical_tests_df,
      file.path(table_dir, paste0(branch_name, "_categorical_tests.rds"))
    )
  }

  # ---------------------------------------------------------------
  # Annotated data bundle
  # ---------------------------------------------------------------
  pd_annotated <- continuous_data
  pd_annotated$cluster_label <- pd_annotated$cluster
  pd_annotated$branch <- branch_name
  pd_annotated$cluster_numeric <- pd_annotated$cluster_id

  saveRDS(pd_annotated, file.path(out_dir, paste0(branch_name, "_annotated_data.rds")))
  write.csv(pd_annotated, file.path(out_dir, paste0(branch_name, "_annotated_data.csv")), row.names = FALSE)

  summary_list <- list(
    branch_name = branch_name,
    branch_model_name = branch_model_name,
    branch_result_file = result_file,
    cluster_sizes = cluster_size_table,
    cluster_means = cluster_means,
    continuous_summary = continuous_summary,
    continuous_tests = continuous_tests_df
  )

  if (exists("categorical_summary") && nrow(categorical_summary) > 0) {
    summary_list$categorical_summary <- categorical_summary
  }
  if (exists("categorical_tests_df") && nrow(categorical_tests_df) > 0) {
    summary_list$categorical_tests <- categorical_tests_df
  }

  saveRDS(summary_list, file.path(out_dir, paste0(branch_name, "_characterisation_summary_list.rds")))

  # ---------------------------------------------------------------
  # Compact branch summary table
  # ---------------------------------------------------------------
  branch_overview <- data.frame(
    branch = branch_name,
    branch_model_name = branch_model_name,
    n_samples = nrow(pd_annotated),
    n_clusters = nlevels(pd_annotated$cluster),
    n_continuous_vars = length(continuous_vars),
    n_categorical_vars = length(categorical_vars),
    stringsAsFactors = FALSE
  )

  write.csv(
    branch_overview,
    file.path(table_dir, paste0(branch_name, "_branch_overview.csv")),
    row.names = FALSE
  )
  saveRDS(branch_overview, file.path(table_dir, paste0(branch_name, "_branch_overview.rds")))

  # Session info for reproducibility
  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, paste0(branch_name, "_session_info.txt"))
  )

  message("Finished sensitivity cluster characterisation for: ", branch_name)
  message("Outputs saved to: ", out_dir)

  invisible(list(
    branch_name = branch_name,
    branch_model_name = branch_model_name,
    result_file = result_file,
    pd_annotated = pd_annotated,
    cluster_sizes = cluster_size_table,
    cluster_means = cluster_means,
    continuous_summary = continuous_summary,
    continuous_tests = continuous_tests_df,
    categorical_summary = categorical_summary,
    categorical_tests = categorical_tests_df,
    summary_list = summary_list
  ))
}

# -------------------------------------------------------------------
# Run all branches
# -------------------------------------------------------------------
branch_results <- list()
combined_overview <- list()
combined_continuous_tests <- list()
combined_categorical_tests <- list()

for (branch_name in names(branch_specs)) {
  spec <- branch_specs[[branch_name]]
  message("Running sensitivity characterisation for branch: ", branch_name)
  res <- process_branch(
    branch_name = branch_name,
    branch_dir = spec$branch_dir,
    result_candidates = spec$result_candidates
  )
  branch_results[[branch_name]] <- res

  if (!is.null(res$continuous_tests) && nrow(res$continuous_tests) > 0) {
    tmp <- res$continuous_tests
    tmp$branch <- branch_name
    combined_continuous_tests[[branch_name]] <- tmp
  }

  if (!is.null(res$categorical_tests) && nrow(res$categorical_tests) > 0) {
    tmp <- res$categorical_tests
    tmp$branch <- branch_name
    combined_categorical_tests[[branch_name]] <- tmp
  }

  combined_overview[[branch_name]] <- data.frame(
    branch = branch_name,
    branch_model_name = res$branch_model_name,
    n_samples = nrow(res$pd_annotated),
    n_clusters = nlevels(res$pd_annotated$cluster),
    n_continuous_vars = nrow(res$continuous_tests),
    n_categorical_vars = if (!is.null(res$categorical_tests)) nrow(res$categorical_tests) else 0L,
    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------------------------
# Combined outputs across branches
# -------------------------------------------------------------------
combined_overview_df <- do.call(rbind, combined_overview)
write.csv(
  combined_overview_df,
  file.path(out_root, "sensitivity_branch_overview.csv"),
  row.names = FALSE
)
saveRDS(combined_overview_df, file.path(out_root, "sensitivity_branch_overview.rds"))

if (length(combined_continuous_tests) > 0) {
  combined_continuous_tests_df <- do.call(rbind, combined_continuous_tests)
  write.csv(
    combined_continuous_tests_df,
    file.path(out_root, "sensitivity_continuous_tests_combined.csv"),
    row.names = FALSE
  )
  saveRDS(
    combined_continuous_tests_df,
    file.path(out_root, "sensitivity_continuous_tests_combined.rds")
  )
}

if (length(combined_categorical_tests) > 0) {
  combined_categorical_tests_df <- do.call(rbind, combined_categorical_tests)
  write.csv(
    combined_categorical_tests_df,
    file.path(out_root, "sensitivity_categorical_tests_combined.csv"),
    row.names = FALSE
  )
  saveRDS(
    combined_categorical_tests_df,
    file.path(out_root, "sensitivity_categorical_tests_combined.rds")
  )
}

saveRDS(branch_results, file.path(out_root, "sensitivity_cluster_characterisation_results.rds"))

message("All sensitivity cluster characterisation complete.")
message("Combined outputs saved to: ", out_root)

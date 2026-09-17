# 06_posthoc_analysis/01_characterise_primary_clusters.R
#
# Characterise the primary integrated SNF clusters.
#

# -----------------------------------------------------------------------------
# 0. Load helpers
# -----------------------------------------------------------------------------

source(file.path(project_root, "06_posthoc_analysis", "00_shared_characterisation_helpers.R"))

# -----------------------------------------------------------------------------
# 1. Load analysis objects
# -----------------------------------------------------------------------------

objs <- load_analysis_objects()
pd_model <- as.data.frame(objs$pd_model)
sample_ids <- as.character(objs$sample_ids)

if (!"Anonymised_sampleID" %in% names(pd_model)) {
  stop("pd_model must contain Anonymised_sampleID.", call. = FALSE)
}

pd_model$sample_id <- as.character(pd_model$Anonymised_sampleID)

# -----------------------------------------------------------------------------
# 2. Load primary SNF outputs
# -----------------------------------------------------------------------------

primary <- load_primary_snf_outputs()

main_network <- as.matrix(primary$fused_network)
main_sample_ids <- as.character(primary$sample_ids)
main_labels_obj <- primary$cluster_labels

check_symmetric_matrix(main_network, "primary fused network")

# Extract cluster labels robustly
if (is.data.frame(main_labels_obj)) {
  label_names <- names(main_labels_obj)

  if (all(c("sample_id", "cluster") %in% label_names)) {
    main_labels_df <- main_labels_obj[, c("sample_id", "cluster"), drop = FALSE]
    names(main_labels_df) <- c("sample_id", "cluster")
  } else if (all(c("sample_id", "label") %in% label_names)) {
    main_labels_df <- main_labels_obj[, c("sample_id", "label"), drop = FALSE]
    names(main_labels_df) <- c("sample_id", "cluster")
  } else if (ncol(main_labels_obj) == 2) {
    main_labels_df <- main_labels_obj
    names(main_labels_df) <- c("sample_id", "cluster")
  } else {
    stop(
      "Primary cluster labels data frame has no obvious sample_id / cluster columns.",
      call. = FALSE
    )
  }

  main_labels_df$sample_id <- as.character(main_labels_df$sample_id)
  main_labels_df$cluster <- as.integer(main_labels_df$cluster)

} else if (is.list(main_labels_obj) && !is.null(main_labels_obj$cluster) && !is.null(main_labels_obj$sample_id)) {
  main_labels_df <- data.frame(
    sample_id = as.character(main_labels_obj$sample_id),
    cluster = as.integer(main_labels_obj$cluster),
    stringsAsFactors = FALSE
  )
} else {
  if (length(main_labels_obj) != length(main_sample_ids)) {
    stop(
      "Primary cluster labels object is not a data frame and does not match sample_ids in length.",
      call. = FALSE
    )
  }
  main_labels_df <- data.frame(
    sample_id = main_sample_ids,
    cluster = as.integer(main_labels_obj),
    stringsAsFactors = FALSE
  )
}

if (anyDuplicated(main_labels_df$sample_id)) {
  stop("Primary cluster label file contains duplicated sample IDs.", call. = FALSE)
}

if (!all(main_labels_df$sample_id %in% sample_ids)) {
  stop("Some primary cluster sample IDs are not present in the analysis cohort.", call. = FALSE)
}

# Reorder clinical data to match cluster-label order
pd_model <- pd_model[match(main_labels_df$sample_id, pd_model$sample_id), , drop = FALSE]

if (!identical(pd_model$sample_id, main_labels_df$sample_id)) {
  stop("Failed to align pd_model to the primary cluster sample order.", call. = FALSE)
}

main_labels <- main_labels_df$cluster
if (length(main_labels) != nrow(pd_model)) {
  stop("Primary cluster labels do not match the number of rows in pd_model.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 3. Output directory
# -----------------------------------------------------------------------------

out_dir <- get_posthoc_output_dir("primary_clusters")
plot_dir <- file.path(out_dir, "plots")
safe_dir(plot_dir)

# -----------------------------------------------------------------------------
# 4. Merge cluster labels onto clinical data
# -----------------------------------------------------------------------------

cluster_membership <- data.frame(
  sample_id = main_labels_df$sample_id,
  primary_cluster = main_labels_df$cluster,
  stringsAsFactors = FALSE
)

saveRDS(cluster_membership, file.path(out_dir, "primary_cluster_membership.rds"))
write.csv(cluster_membership, file.path(out_dir, "primary_cluster_membership.csv"), row.names = FALSE)

pd_annotated <- merge(
  pd_model,
  cluster_membership,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(pd_annotated$primary_cluster)) {
  stop("Some primary cluster labels are missing after merge.", call. = FALSE)
}

pd_annotated$primary_cluster <- factor(pd_annotated$primary_cluster, levels = sort(unique(main_labels)))

# -----------------------------------------------------------------------------
# 5. Helper functions for variable handling
# -----------------------------------------------------------------------------

pick_first_available <- function(candidates, data) {
  found <- candidates[candidates %in% names(data)]
  if (length(found) == 0) return(NA_character_)
  found[1]
}

to_numeric_safe <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  suppressWarnings(as.numeric(x))
}

# Use raw variables where available, otherwise fall back to numeric versions
raw_continuous_vars <- c(
  pick_first_available(c("AGE", "AGE_num"), pd_annotated),
  pick_first_available(c("UPDRS_III", "UPDRS_III_num"), pd_annotated),
  pick_first_available(c("MOCA_total", "MOCA_total_num"), pd_annotated),
  pick_first_available(c("LEDD_total", "LEDD_total_num"), pd_annotated),
  pick_first_available(c("disease_duration"), pd_annotated)
)
raw_continuous_vars <- raw_continuous_vars[!is.na(raw_continuous_vars)]

categorical_vars <- c(
  pick_first_available(c("GENDER"), pd_annotated)
)
categorical_vars <- categorical_vars[!is.na(categorical_vars)]

continuous_data <- pd_annotated[, c("sample_id", raw_continuous_vars, "primary_cluster"), drop = FALSE]
for (v in raw_continuous_vars) {
  continuous_data[[v]] <- to_numeric_safe(continuous_data[[v]])
}

continuous_vars <- raw_continuous_vars

# Remove any variables that are entirely missing after coercion
keep_continuous <- vapply(continuous_data[continuous_vars], function(z) any(is.finite(z)), logical(1))
continuous_vars <- continuous_vars[keep_continuous]
continuous_data <- continuous_data[, c("sample_id", continuous_vars, "primary_cluster"), drop = FALSE]

if (length(continuous_vars) == 0) {
  stop("No usable continuous clinical variables were found.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 6. Cluster sizes
# -----------------------------------------------------------------------------

cluster_size_table <- as.data.frame(table(pd_annotated$primary_cluster), stringsAsFactors = FALSE)
names(cluster_size_table) <- c("primary_cluster", "n")
cluster_size_table$proportion <- cluster_size_table$n / sum(cluster_size_table$n)

write.csv(cluster_size_table, file.path(out_dir, "primary_cluster_sizes.csv"), row.names = FALSE)
saveRDS(cluster_size_table, file.path(out_dir, "primary_cluster_sizes.rds"))

# -----------------------------------------------------------------------------
# 7. Continuous-variable summaries
# -----------------------------------------------------------------------------

summarise_continuous_by_cluster <- function(x, g) {
  g <- as.factor(g)
  out <- aggregate(
    x,
    by = list(cluster = g),
    FUN = function(z) {
      c(
        n = sum(is.finite(z)),
        mean = mean(z, na.rm = TRUE),
        sd = stats::sd(z, na.rm = TRUE),
        median = stats::median(z, na.rm = TRUE),
        q25 = stats::quantile(z, 0.25, na.rm = TRUE, names = FALSE),
        q75 = stats::quantile(z, 0.75, na.rm = TRUE, names = FALSE)
      )
    }
  )
  out <- do.call(data.frame, out)
  out
}

continuous_summary_list <- list()

for (v in continuous_vars) {
  tmp <- summarise_continuous_by_cluster(continuous_data[[v]], continuous_data$primary_cluster)
  tmp$variable <- v
  continuous_summary_list[[v]] <- tmp
}

continuous_summary <- do.call(rbind, continuous_summary_list)
rownames(continuous_summary) <- NULL
names(continuous_summary) <- c("cluster", "n", "mean", "sd", "median", "q25", "q75", "variable")
continuous_summary <- continuous_summary[, c("variable", "cluster", "n", "mean", "sd", "median", "q25", "q75")]

write.csv(continuous_summary, file.path(out_dir, "primary_continuous_summary_by_cluster.csv"), row.names = FALSE)
saveRDS(continuous_summary, file.path(out_dir, "primary_continuous_summary_by_cluster.rds"))

# -----------------------------------------------------------------------------
# 8. Continuous-variable tests
# -----------------------------------------------------------------------------

cohens_d <- function(x, g) {
  x <- as.numeric(x)
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

continuous_tests <- list()

for (v in continuous_vars) {
  x <- continuous_data[[v]]
  g <- continuous_data$primary_cluster

  test_name <- NA_character_
  statistic <- NA_real_
  p_value <- NA_real_
  estimate <- NA_real_
  conf_low <- NA_real_
  conf_high <- NA_real_
  effect_size <- NA_real_

  if (nlevels(g) == 2) {
    wt <- suppressWarnings(stats::wilcox.test(x ~ g, exact = FALSE, conf.int = TRUE))
    test_name <- "wilcox.test"
    statistic <- unname(wt$statistic)
    p_value <- wt$p.value
    if (!is.null(wt$estimate)) estimate <- unname(wt$estimate)
    if (!is.null(wt$conf.int) && length(wt$conf.int) == 2) {
      conf_low <- wt$conf.int[1]
      conf_high <- wt$conf.int[2]
    }
    effect_size <- cohens_d(x, g)
  } else if (nlevels(g) > 2) {
    kt <- suppressWarnings(stats::kruskal.test(x ~ g))
    test_name <- "kruskal.test"
    statistic <- unname(kt$statistic)
    p_value <- kt$p.value
  }

  continuous_tests[[v]] <- data.frame(
    variable = v,
    test = test_name,
    statistic = statistic,
    p_value = p_value,
    hl_estimate = estimate,
    hl_conf_low = conf_low,
    hl_conf_high = conf_high,
    cohens_d = effect_size,
    stringsAsFactors = FALSE
  )
}

continuous_tests_df <- do.call(rbind, continuous_tests)
continuous_tests_df$padj_bh <- p.adjust(continuous_tests_df$p_value, method = "BH")

write.csv(continuous_tests_df, file.path(out_dir, "primary_continuous_tests.csv"), row.names = FALSE)
saveRDS(continuous_tests_df, file.path(out_dir, "primary_continuous_tests.rds"))

# -----------------------------------------------------------------------------
# 9. Categorical-variable summaries and tests
# -----------------------------------------------------------------------------

cramers_v <- function(tab) {
  tab <- as.matrix(tab)
  n <- sum(tab)
  if (n == 0) return(NA_real_)
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  r <- nrow(tab)
  c <- ncol(tab)
  denom <- min(r - 1, c - 1)
  if (denom <= 0) return(NA_real_)
  sqrt(as.numeric(chi$statistic) / (n * denom))
}

categorical_summary_list <- list()
categorical_tests_list <- list()

for (v in categorical_vars) {
  x <- as.factor(pd_annotated[[v]])
  g <- pd_annotated$primary_cluster

  tab <- table(cluster = g, level = x)

  tab_df <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(tab_df) <- c("cluster", "level", "n")
  tab_df$variable <- v
  tab_df$proportion_within_cluster <- ave(tab_df$n, tab_df$cluster, FUN = function(z) z / sum(z))

  categorical_summary_list[[v]] <- tab_df

  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  test_name <- "chisq.test"
  test_obj <- chi

  if (any(chi$expected < 5)) {
    test_name <- "fisher.test"
    test_obj <- tryCatch(
      stats::fisher.test(tab),
      error = function(e) chi
    )
  }

  p_value <- if (!is.null(test_obj$p.value)) test_obj$p.value else NA_real_
  statistic <- if (!is.null(test_obj$statistic)) unname(test_obj$statistic) else NA_real_
  vcramer <- cramers_v(tab)

  categorical_tests_list[[v]] <- data.frame(
    variable = v,
    test = test_name,
    statistic = statistic,
    p_value = p_value,
    cramers_v = vcramer,
    stringsAsFactors = FALSE
  )
}

if (length(categorical_summary_list) > 0) {
  categorical_summary <- do.call(rbind, categorical_summary_list)
  rownames(categorical_summary) <- NULL
  write.csv(categorical_summary, file.path(out_dir, "primary_categorical_summary_by_cluster.csv"), row.names = FALSE)
  saveRDS(categorical_summary, file.path(out_dir, "primary_categorical_summary_by_cluster.rds"))
}

if (length(categorical_tests_list) > 0) {
  categorical_tests_df <- do.call(rbind, categorical_tests_list)
  categorical_tests_df$padj_bh <- p.adjust(categorical_tests_df$p_value, method = "BH")
  write.csv(categorical_tests_df, file.path(out_dir, "primary_categorical_tests.csv"), row.names = FALSE)
  saveRDS(categorical_tests_df, file.path(out_dir, "primary_categorical_tests.rds"))
}

# -----------------------------------------------------------------------------
# 10. Cluster means table
# -----------------------------------------------------------------------------

cluster_means <- aggregate(
  continuous_data[, continuous_vars, drop = FALSE],
  by = list(primary_cluster = continuous_data$primary_cluster),
  FUN = function(z) mean(z, na.rm = TRUE)
)

write.csv(cluster_means, file.path(out_dir, "primary_cluster_means.csv"), row.names = FALSE)
saveRDS(cluster_means, file.path(out_dir, "primary_cluster_means.rds"))

# -----------------------------------------------------------------------------
# 11. Plots
# -----------------------------------------------------------------------------

png(file.path(plot_dir, "primary_cluster_sizes.png"), width = 1800, height = 1200, res = 250)
barplot(
  height = cluster_size_table$n,
  names.arg = cluster_size_table$primary_cluster,
  xlab = "Primary cluster",
  ylab = "Number of patients",
  main = "Primary SNF cluster sizes"
)
dev.off()

for (v in continuous_vars) {
  png(file.path(plot_dir, paste0("primary_", v, "_by_cluster.png")), width = 2000, height = 1500, res = 250)
  boxplot(
    continuous_data[[v]] ~ continuous_data$primary_cluster,
    xlab = "Primary cluster",
    ylab = v,
    main = paste0(v, " by primary cluster")
  )
  dev.off()
}

if ("GENDER" %in% categorical_vars) {
  gender_tab <- table(pd_annotated$primary_cluster, pd_annotated$GENDER)
  png(file.path(plot_dir, "primary_gender_by_cluster.png"), width = 2000, height = 1500, res = 250)
  barplot(
    gender_tab,
    beside = TRUE,
    legend.text = TRUE,
    xlab = "Primary cluster",
    ylab = "Count",
    main = "Gender distribution by primary cluster"
  )
  dev.off()
}

# -----------------------------------------------------------------------------
# 12. Save annotated data and summary bundle
# -----------------------------------------------------------------------------

saveRDS(pd_annotated, file.path(out_dir, "primary_cluster_annotated_data.rds"))
write.csv(pd_annotated, file.path(out_dir, "primary_cluster_annotated_data.csv"), row.names = FALSE)

summary_list <- list(
  cluster_sizes = cluster_size_table,
  continuous_summary = continuous_summary,
  continuous_tests = continuous_tests_df
)

if (exists("categorical_summary")) {
  summary_list$categorical_summary <- categorical_summary
}
if (exists("categorical_tests_df")) {
  summary_list$categorical_tests <- categorical_tests_df
}

saveRDS(summary_list, file.path(out_dir, "primary_characterisation_summary_list.rds"))

message("Primary cluster characterisation complete.")
message("Outputs saved to: ", out_dir)
message("Samples: ", nrow(pd_annotated))
message("Clusters: ", nlevels(pd_annotated$primary_cluster))
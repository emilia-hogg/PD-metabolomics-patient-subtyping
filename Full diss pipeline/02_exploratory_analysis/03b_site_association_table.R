# ============================================================
# 03b_site_association_table.R
# Builds from the contingency tables already written by 03_kmeans.R. 
# Produces: <exploratory>/kmeans/site_association_all.csv
# ============================================================

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))

kmeans_root <- file.path(paths$exploratory, "kmeans")

files_counts <- list.files(kmeans_root, pattern = "_site_cluster_counts_k[0-9]+\\.csv$",
                           recursive = TRUE, full.names = TRUE)
if (!length(files_counts)) stop("No site_cluster_counts files found under ", kmeans_root)
message("Found ", length(files_counts), " contingency tables.")

set.seed(1)   # the Monte Carlo p-values below are stochastic

rows <- lapply(files_counts, function(f) {
  bn <- basename(f)
  k  <- as.integer(sub(".*_k([0-9]+)\\.csv$", "\\1", bn))
  lb <- sub("_site_cluster_counts_k[0-9]+\\.csv$", "", bn)

  # tidy the label: PCA_PC05 -> PC5, UMAP -> UMAP
  space <- sub("^PCA_PC0?", "PC", lb)

  d <- read.csv(f, row.names = 1, check.names = FALSE)
  m <- as.matrix(d)
  m <- m[rowSums(m) > 0, colSums(m) > 0, drop = FALSE]   # drop empty sites/clusters
  n <- sum(m)

  cs  <- suppressWarnings(chisq.test(m))
  # Monte Carlo p-value: many sites are very small, so a large share of
  # expected counts fall below 5 and the asymptotic p-value is unreliable
  csm <- suppressWarnings(chisq.test(m, simulate.p.value = TRUE, B = 10000))

  cramers_v <- sqrt(as.numeric(cs$statistic) / (n * (min(dim(m)) - 1)))
  exp_lt5   <- mean(cs$expected < 5) * 100

  data.frame(
    space          = space,
    k              = k,
    n_sites        = nrow(m),
    n_patients     = n,
    chi_squared    = as.numeric(cs$statistic),
    df             = as.numeric(cs$parameter),
    p_asymptotic   = as.numeric(cs$p.value),
    p_monte_carlo  = as.numeric(csm$p.value),
    cramers_v      = cramers_v,
    pct_expected_lt5 = exp_lt5,
    stringsAsFactors = FALSE
  )
})

res <- do.call(rbind, rows)

# order PC5 < PC10 < ... < PC50 < UMAP, then by k
num <- suppressWarnings(as.integer(sub("^PC", "", res$space)))
res$space <- factor(res$space,
  levels = c(paste0("PC", sort(unique(num[!is.na(num)]))), "UMAP"))
res <- res[order(res$space, res$k), ]

# correct across the whole family of tests, not per test
res$p_adj_asymptotic  <- p.adjust(res$p_asymptotic,  method = "BH")
res$p_adj_monte_carlo <- p.adjust(res$p_monte_carlo, method = "BH")

out <- file.path(kmeans_root, "site_association_all.csv")
write.csv(res, out, row.names = FALSE)

# ---- printed summary ----------------------------------------
cat("\n")
cat(sprintf("%-6s %3s %6s %10s %5s %11s %11s %9s %9s %7s\n",
            "space", "k", "sites", "chi-sq", "df", "p", "p-adj", "p (MC)",
            "p-adj MC", "V"))
for (i in seq_len(nrow(res))) {
  cat(sprintf("%-6s %3d %6d %10.2f %5d %11.3g %11.3g %9.4f %9.4f %7.3f\n",
              as.character(res$space[i]), res$k[i], res$n_sites[i],
              res$chi_squared[i], res$df[i],
              res$p_asymptotic[i], res$p_adj_asymptotic[i],
              res$p_monte_carlo[i], res$p_adj_monte_carlo[i],
              res$cramers_v[i]))
}

cat("\n--- summary ---\n")
cat("total tests:                      ", nrow(res), "\n")
cat("sites contributing:               ", paste(sort(unique(res$n_sites)), collapse = ", "), "\n")
cat("significant, uncorrected (0.05):  ", sum(res$p_asymptotic < 0.05), "\n")
cat("significant, BH-adjusted:         ", sum(res$p_adj_asymptotic < 0.05), "\n")
cat("significant, BH-adjusted (MC):    ", sum(res$p_adj_monte_carlo < 0.05), "\n")
cat(sprintf("Cramer's V range:                  %.3f to %.3f\n",
            min(res$cramers_v), max(res$cramers_v)))
cat(sprintf("mean %% of expected counts below 5: %.1f%%\n", mean(res$pct_expected_lt5)))
cat("\nWritten to: ", out, "\n")

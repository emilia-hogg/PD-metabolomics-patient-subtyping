# 04_snf/02_run_snf.R
# Re-runnable from a fresh R session.

project_root <- "/home/ehogg/analysis/Baseline-Tracking-PD-Metabolite-Analysis-Natacha-/Full diss pipeline and outputs"
source(file.path(project_root, "00_config", "config.R"))
source(file.path(project_root, "00_config", "packages.R"))

if (!requireNamespace("SNFtool", quietly = TRUE)) {
  stop("SNFtool is required but is not available.")
}

# -------------------------------------------------------------------
# Basic checks on configuration
# -------------------------------------------------------------------
required_param_names <- c("K", "sigma", "t")
missing_param_names <- setdiff(required_param_names, names(snf_params))
if (length(missing_param_names) > 0) {
  stop("snf_params is missing required element(s): ", paste(missing_param_names, collapse = ", "))
}

if (!is.numeric(snf_params$K) || length(snf_params$K) != 1 || snf_params$K <= 0) {
  stop("snf_params$K must be a single positive number.")
}
if (!is.numeric(snf_params$sigma) || length(snf_params$sigma) != 1 || snf_params$sigma <= 0) {
  stop("snf_params$sigma must be a single positive number.")
}
if (!is.numeric(snf_params$t) || length(snf_params$t) != 1 || snf_params$t <= 0) {
  stop("snf_params$t must be a single positive number.")
}

K <- as.integer(snf_params$K)
sigma <- as.numeric(snf_params$sigma)
t <- as.integer(snf_params$t)

# sigma is not used directly by SNF() here because it is already baked
# into the upstream affinity matrices. Keep it only for metadata/logging.

# -------------------------------------------------------------------
# Locate input folders robustly
# -------------------------------------------------------------------
clinical_base <- file.path(project_root, "03_network_construction", "outputs")
metabolomics_base <- clinical_base
snf_out_base <- paths$snf

find_existing_dir_with_file <- function(base_dir, subdirs, filename) {
  candidates <- c(file.path(base_dir, subdirs), base_dir)
  for (cand in candidates) {
    if (file.exists(file.path(cand, filename))) {
      return(cand)
    }
  }
  stop("Could not find ", filename, " under: ", paste(candidates, collapse = " | "))
}

clinical_dir <- find_existing_dir_with_file(
  base_dir = clinical_base,
  subdirs = c("clinical_network"),
  filename = "updrs_affinity_mat.rds"
)

metabolomics_dir <- find_existing_dir_with_file(
  base_dir = metabolomics_base,
  subdirs = c("metabolomics_network"),
  filename = "metab_affinity_mat.rds"
)

out_dir <- file.path(snf_out_base, "main_four_network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------
read_square_matrix <- function(path, object_name) {
  mat <- readRDS(path)
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop(object_name, " is empty: ", path)
  }
  if (nrow(mat) != ncol(mat)) {
    stop(object_name, " must be square: ", path)
  }
  if (anyNA(mat)) {
    stop(object_name, " contains NA values: ", path)
  }
  if (!isTRUE(all.equal(mat, t(mat)))) {
    stop(object_name, " is not symmetric: ", path)
  }
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop(object_name, " must have row and column names: ", path)
  }
  if (!identical(rownames(mat), colnames(mat))) {
    stop(object_name, " row and column names do not match: ", path)
  }
  mat
}

check_common_sample_order <- function(mats) {
  ref_ids <- rownames(mats[[1]])
  for (nm in names(mats)) {
    if (!identical(rownames(mats[[nm]]), ref_ids)) {
      stop("Row names do not match across matrices. Mismatch found in: ", nm)
    }
    if (!identical(colnames(mats[[nm]]), ref_ids)) {
      stop("Column names do not match across matrices. Mismatch found in: ", nm)
    }
  }
  ref_ids
}

save_fused_qc <- function(W, out_dir, prefix) {
  offdiag <- W[upper.tri(W)]

  saveRDS(W, file.path(out_dir, paste0(prefix, "_fused_network.rds")))
  saveRDS(rownames(W), file.path(out_dir, paste0(prefix, "_sample_ids.rds")))

  write.csv(
    data.frame(
      metric = c(
        "n_samples",
        "n_edges_upper_triangle",
        "min_value",
        "max_value",
        "median_offdiag",
        "mean_offdiag",
        "q95_offdiag",
        "n_above_001",
        "prop_above_001"
      ),
      value = c(
        nrow(W),
        length(offdiag),
        min(W, na.rm = TRUE),
        max(W, na.rm = TRUE),
        median(offdiag, na.rm = TRUE),
        mean(offdiag, na.rm = TRUE),
        as.numeric(quantile(offdiag, 0.95, na.rm = TRUE, names = FALSE)),
        sum(offdiag > 0.01, na.rm = TRUE),
        mean(offdiag > 0.01, na.rm = TRUE)
      ),
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, paste0(prefix, "_summary.csv")),
    row.names = FALSE
  )

  png(file.path(out_dir, paste0(prefix, "_histogram.png")), width = 2000, height = 1500, res = 300)
  hist(
    offdiag,
    breaks = 50,
    main = paste0("Fused network off-diagonal values - ", prefix),
    xlab = "Similarity"
  )
  dev.off()

  invisible(TRUE)
}

# -------------------------------------------------------------------
# Load affinity matrices
# -------------------------------------------------------------------
message("Loading affinity matrices...")

updrs_affinity_mat <- read_square_matrix(
  file.path(clinical_dir, "updrs_affinity_mat.rds"),
  "UPDRS affinity matrix"
)
moca_affinity_mat <- read_square_matrix(
  file.path(clinical_dir, "moca_affinity_mat.rds"),
  "MoCA affinity matrix"
)
ledd_affinity_mat <- read_square_matrix(
  file.path(clinical_dir, "ledd_affinity_mat.rds"),
  "LEDD affinity matrix"
)
metab_affinity_mat <- read_square_matrix(
  file.path(metabolomics_dir, "metab_affinity_mat.rds"),
  "Metabolomics affinity matrix"
)

mats <- list(
  UPDRS = updrs_affinity_mat,
  MoCA = moca_affinity_mat,
  LEDD = ledd_affinity_mat,
  Metabolomics = metab_affinity_mat
)

sample_ids <- check_common_sample_order(mats)
message("All affinity matrices share the same sample order: ", length(sample_ids), " samples.")

# -------------------------------------------------------------------
# Basic sanity checks
# -------------------------------------------------------------------
for (nm in names(mats)) {
  mat <- mats[[nm]]
  if (any(mat < 0, na.rm = TRUE)) {
    stop(nm, " contains negative affinity values.")
  }
}

# -------------------------------------------------------------------
# Run the primary four-network SNF
# -------------------------------------------------------------------
message("Running primary four-network SNF with K = ", K, ", sigma = ", sigma, ", t = ", t, "...")

fused_network <- SNFtool::SNF(
  list(
    updrs_affinity_mat,
    moca_affinity_mat,
    ledd_affinity_mat,
    metab_affinity_mat
  ),
  K = K,
  t = t
)

if (!is.matrix(fused_network)) {
  fused_network <- as.matrix(fused_network)
}
rownames(fused_network) <- sample_ids
colnames(fused_network) <- sample_ids

if (nrow(fused_network) != length(sample_ids) || ncol(fused_network) != length(sample_ids)) {
  stop("Fused network has unexpected dimensions.")
}
if (anyNA(fused_network)) {
  stop("Fused network contains NA values.")
}
if (!isTRUE(all.equal(fused_network, t(fused_network)))) {
  stop("Fused network is not symmetric.")
}

# -------------------------------------------------------------------
# Save outputs
# -------------------------------------------------------------------
prefix <- "primary_four_network"

saveRDS(fused_network, file.path(out_dir, paste0(prefix, "_fused_network.rds")))
saveRDS(sample_ids, file.path(out_dir, paste0(prefix, "_sample_ids.rds")))
saveRDS(
  list(
    K = K,
    sigma = sigma,
    t = t,
    modality_names = names(mats),
    sample_ids = sample_ids
  ),
  file.path(out_dir, paste0(prefix, "_run_metadata.rds"))
)

save_fused_qc(fused_network, out_dir, prefix)

write.csv(
  fused_network,
  file.path(out_dir, paste0(prefix, "_fused_network.csv")),
  row.names = TRUE
)

message("Primary four-network SNF complete.")
message("Outputs saved to: ", out_dir)
message("Samples: ", nrow(fused_network))
message("Parameters: K = ", K, ", sigma = ", sigma, ", t = ", t)
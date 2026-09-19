# 00_config/packages.R
# This script defines the R packages required for the analysis pipeline.
# For each package, it checks whether the package is installed, installs it
# from CRAN if necessary, and then loads it while suppressing startup messages.
#
# No file paths or other parameters need to be set before running this script.

required_packages <- c(
  "SNFtool",
  "ggplot2",
  "dplyr",
  "tibble",
  "cluster",
  "dbscan",
  "readr",
  "purrr",
  "stringr",
  "tidyr",
  "forcats",
  "ggrepel"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

invisible(lapply(required_packages, install_if_missing))

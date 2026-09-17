# 00_config/packages.R

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
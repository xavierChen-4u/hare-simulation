#!/usr/bin/env Rscript

# Install only packages that are currently missing. This does not pin versions.
# Run scripts/check_environment.R to record the environment used for a run.
packages <- c("gbm", "splitTools", "dplyr", "future", "future.apply", "ggplot2")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

still_missing <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(still_missing)) {
  stop("Packages still unavailable: ", paste(still_missing, collapse = ", "))
}
message("All required packages are available.")

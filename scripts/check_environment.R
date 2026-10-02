#!/usr/bin/env Rscript

# This checks package availability without running a simulation.
# A report is written even when packages are missing.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) {
  stop("Run this script with Rscript scripts/check_environment.R")
}
repo_dir <- normalizePath(
  file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."),
  mustWork = TRUE
)
output_dir <- file.path(repo_dir, "generated")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

packages <- c("gbm", "splitTools", "dplyr", "future", "future.apply", "ggplot2")
available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
versions <- vapply(packages, function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    "MISSING"
  }
}, character(1))

report <- c(
  paste("Environment recorded at:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "",
  "Package versions:",
  paste(names(versions), versions, sep = ": "),
  "",
  capture.output(sessionInfo())
)
report_path <- file.path(output_dir, "session-info.txt")
writeLines(report, report_path)
message("Environment report: ", report_path)
if (any(!available)) {
  stop(
    "Missing packages: ", paste(packages[!available], collapse = ", "),
    ". Run Rscript scripts/install_dependencies.R first."
  )
}
message("All required packages are available; no simulation was run.")

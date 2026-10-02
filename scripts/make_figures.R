#!/usr/bin/env Rscript

# Reconstruct the contents of manuscript Figures 1--3 from a summary CSV.
# Styling is new; these are not claimed to be the manuscript's original files.
# Figure 1 error bars are mean +/- one empirical Monte Carlo SD, not CIs.
# Usage (paths supplied by the caller are relative to the working directory):
#   Rscript scripts/make_figures.R
#   Rscript scripts/make_figures.R path/to/lambda_transport_summary.csv
#   Rscript scripts/make_figures.R path/to/summary.csv path/to/output_directory

required <- c("dplyr", "ggplot2")
missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       ". Run Rscript scripts/install_dependencies.R first.")
}
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) {
  stop("Run this script with Rscript scripts/make_figures.R")
}
repo_dir <- normalizePath(
  file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."),
  mustWork = TRUE
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) {
  stop("Usage: Rscript scripts/make_figures.R [summary.csv] [output_directory]")
}
summary_path <- if (length(args) >= 1L) args[[1L]] else {
  file.path(repo_dir, "data", "reference", "lambda_transport_summary.csv")
}
output_dir <- if (length(args) >= 2L) args[[2L]] else {
  file.path(repo_dir, "figures")
}
summary_path <- normalizePath(summary_path, mustWork = TRUE)
results <- read.csv(summary_path, stringsAsFactors = FALSE, check.names = FALSE)

numeric_columns <- c(
  "lambda", "label_budget", "successful_replications", "mean_estimate",
  "bias", "empirical_sd", "rmse"
)
missing_columns <- setdiff(c("method", numeric_columns), names(results))
if (length(missing_columns)) {
  stop("Summary CSV is missing: ", paste(missing_columns, collapse = ", "))
}
for (column in numeric_columns) {
  results[[column]] <- suppressWarnings(as.numeric(results[[column]]))
  if (any(!is.finite(results[[column]]))) {
    stop("Column contains missing or nonnumeric values: ", column)
  }
}
if (any(results$empirical_sd < 0 | results$rmse <= 0 |
        results$successful_replications < 2)) {
  stop("Expected nonnegative SD, positive RMSE, and at least two replications.")
}

method_names <- c(
  "Bridge-Constant" = "HARE-C",
  "Kallus-GBM" = "SA-EIF-GBM",
  "OI-GBM" = "SI-AIPW-GBM",
  "Kallus Oracle" = "SA-EIF-Oracle"
)
budget_levels <- c(100, 200, 400, 800, 1600)
lambda_levels <- c(0, 0.5, 1)
keys <- c("lambda", "label_budget", "method")
expected <- expand.grid(
  lambda = lambda_levels,
  label_budget = budget_levels,
  method = names(method_names),
  stringsAsFactors = FALSE
)
results <- semi_join(results, expected, by = keys)
if (any(duplicated(results[keys]))) {
  stop("Summary CSV has duplicate (lambda, label_budget, method) rows.")
}
missing_rows <- anti_join(expected, results, by = keys)
if (nrow(missing_rows)) {
  stop("Summary CSV does not contain the complete 3 x 5 x 4 manuscript grid.")
}
tau_true <- 2
if (any(abs(results$mean_estimate - results$bias - tau_true) > 1e-8)) {
  stop("Summary estimates and biases do not imply the manuscript's true ATE of 2.")
}
results$method_paper <- factor(
  unname(method_names[results$method]), levels = unname(method_names)
)
results$budget <- factor(
  results$label_budget, levels = budget_levels,
  labels = c("100", "200", "400", "800", "1,600")
)

theme_paper <- theme_classic(base_size = 11, base_family = "serif") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.caption = element_text(size = 9, hjust = 0),
    plot.margin = margin(8, 10, 8, 8)
  )
method_colors <- c(
  "HARE-C" = "#0072B2", "SA-EIF-GBM" = "#D55E00",
  "SI-AIPW-GBM" = "#009E73", "SA-EIF-Oracle" = "#7A5195"
)

# Figure 1: Monte Carlo mean +/- one empirical SD at lambda = 0, B = 800.
fig1_data <- filter(
  results, lambda == 0, label_budget == 800,
  method %in% c("OI-GBM", "Bridge-Constant")
)
fig1_data$method_paper <- factor(
  fig1_data$method_paper, levels = c("SI-AIPW-GBM", "HARE-C")
)
figure1 <- ggplot(fig1_data, aes(method_paper, mean_estimate, color = method_paper)) +
  geom_hline(yintercept = tau_true, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(
    aes(ymin = mean_estimate - empirical_sd, ymax = mean_estimate + empirical_sd),
    width = 0.15, linewidth = 0.7
  ) +
  geom_point(size = 3) +
  scale_color_manual(values = method_colors, guide = "none") +
  labs(
    x = NULL, y = "ATE estimate",
    title = "Constant discrepancy, expected label budget = 800",
    caption = "Points: Monte Carlo means. Bars: +/- one empirical SD (not confidence intervals).\nDashed line: true ATE = 2."
  ) + theme_paper

# Figure 2: all four methods, with the B = 800 / 1,600 comparison marked.
fig2_data <- filter(results, lambda == 0)
highlight <- filter(
  fig2_data,
  (method == "Bridge-Constant" & label_budget == 800) |
    (method == "Kallus-GBM" & label_budget == 1600)
)
figure2 <- ggplot(
  fig2_data,
  aes(budget, rmse, color = method_paper, linetype = method_paper,
      shape = method_paper, group = method_paper)
) +
  geom_line(linewidth = 0.7) + geom_point(size = 2.5) +
  geom_point(
    data = highlight, aes(budget, rmse), inherit.aes = FALSE,
    shape = 21, fill = NA, color = "black", size = 5, stroke = 0.8
  ) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
  labs(
    x = "Expected target-label budget", y = "ATE RMSE",
    color = NULL, linetype = NULL, shape = NULL,
    title = "RMSE under constant discrepancy",
    caption = "Outlined points: HARE-C at budget 800 and SA-EIF-GBM at budget 1,600."
  ) + theme_paper + guides(color = guide_legend(nrow = 2, byrow = TRUE))

# Figure 3: a ratio of RMSE summaries, not a mean of per-replication ratios.
hare <- select(
  filter(results, method == "Bridge-Constant"),
  lambda, label_budget, budget, hare_rmse = rmse
)
baseline <- select(
  filter(results, method == "Kallus-GBM"),
  lambda, label_budget, baseline_rmse = rmse
)
fig3_data <- inner_join(hare, baseline, by = c("lambda", "label_budget"))
fig3_data$relative_rmse <- fig3_data$hare_rmse / fig3_data$baseline_rmse
figure3 <- ggplot(
  fig3_data,
  aes(lambda, relative_rmse, color = budget, linetype = budget,
      shape = budget, group = budget)
) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 0.7) + geom_point(size = 2.5) +
  scale_x_continuous(breaks = lambda_levels) +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73", "#7A5195", "#333333")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash")) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8)) +
  labs(
    x = "Discrepancy heterogeneity (lambda)",
    y = "RMSE(HARE-C) / RMSE(SA-EIF-GBM)",
    color = "Expected target-label budget", linetype = "Expected target-label budget",
    shape = "Expected target-label budget",
    caption = "Values below 1 indicate lower RMSE for HARE-C."
  ) + theme_paper

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
plots <- list(
  bias_budget_800 = figure1,
  rmse_label_budget_curve = figure2,
  heterogeneity_relative_rmse = figure3
)
for (name in names(plots)) {
  for (extension in c("pdf", "png")) {
    ggsave(
      file.path(output_dir, paste0(name, ".", extension)),
      plot = plots[[name]], width = 7, height = 4.7, units = "in",
      dpi = 300, bg = "white"
    )
  }
}
message("Input: ", summary_path)
message("Figures saved to: ", normalizePath(output_dir, mustWork = TRUE))

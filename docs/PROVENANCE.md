# Simulation provenance and reproducibility

Audit date: 2026-10-02. Manuscript: `Siyu_Chen_Writing_Sample.pdf`, 24 pages, by Siyu Chen. The audit focused on Section 5 (pages 16-23), its numerical claims elsewhere in the paper, and the supplied `04_cloud_run.R` project.

## Finding

The current manuscript's reported simulation numbers are consistent with the archived results associated with `04_cloud_run.R`. All numerical cells in Tables 6-8 match. The earlier exploratory notebooks are not required to produce the reported Section 5 results.

This conclusion is based on the preserved source, archived run artifacts, independent aggregation of replication-level results, workbook comparisons, and inspection of the paper. It is not a claim that the original Monte Carlo run was independently repeated in its original software environment.

## Evidence chain

1. The supplied runner and the runner bundled with `results/lambda_transport_full_results/` are byte-for-byte identical. SHA-256: `67657cb4b7f24c58fc25a41c4b74318a8c21aa5e85d215f62d8f77e833498033`.
2. The archived run contains 750 scenario/budget/replication checkpoint files, 3,000 method-level result rows, and 60 summary rows. Each method/configuration contains exactly 50 successful replications. The cloud log identifies a `full` run using 24 workers. It lists 249 pending bundles; the full results contain 250 budget/replication bundles, consistent with reuse of an existing checkpoint bundle, though the log alone does not establish that bundle's earlier history.
3. The workbook `surrogate_experiment_results_summary.xlsx`, sheet `05_正式异质d`, explicitly names `lambda_transport_summary.csv` and `04_cloud_run.R` as its sources. All 60 x 24 data cells match the CSV: 1,230 numeric, 120 text, and 90 missing-value cells. Missing values are displayed as `n.a.` in the workbook. The workbook's SHA-256 is `1028e1a136bc8d8cb5574b8fc9ccde912f6d0e169bf0addac85bd0f730aa6187`.
4. The legacy `make_hare_figures.R/figures.R` reads that sheet and maps the original method identifiers to the paper's names. Its filename expects a workbook ending in `(2).xlsx`, which was not present beside that script. The workbook copy under the figure directory was locked during this audit; the accessible parent simulation-directory workbook was checked instead.
5. Independently recomputing the archived raw results reproduces all 1,200 summary metrics, all 45 printed numerical result cells in Tables 6-8, eight additional independent text values, and the within-budget and cross-budget RMSE comparisons.

Hashes for the preserved runner and the three included reference CSV files are in `data/reference/SHA256SUMS.json`. Individual simulated participant-level datasets were generated in memory; the saved result CSV is replication-level output, not those datasets.

## Manuscript mapping

| Manuscript material | Source / interpretation |
| --- | --- |
| Table 1 | Conceptual observed-data structure; not a Monte Carlo result table |
| Tables 2-3, Section 5.1 | Runner settings and data-generating equations; Table 3 reports population design values |
| Tables 4-5 | Estimator construction and learner settings in the runner |
| Table 6, p. 19 | `lambda = 0`; biases of `OI-GBM` and `Bridge-Constant`; reduction computed from unrounded bias |
| Table 7, p. 20 | `lambda = 0`; RMSE of `Kallus-GBM` and `Bridge-Constant`; reduction computed from unrounded RMSE |
| Table 8, p. 22 | Budget 800; HARE-C bias, sample SD and RMSE at all three lambda values; comparator RMSE and relative reduction |
| Figure 1, p. 20 | Mean estimates plus/minus empirical Monte Carlo SD, budget 800 and lambda 0 |
| Figure 2, p. 21 | RMSE of all four methods at lambda 0 across the five budgets |
| Figure 3, p. 23 | HARE-C RMSE divided by SA-EIF-GBM RMSE, across lambda for each budget |
| Abstract, introduction, conclusion | Repeated headline values from this same experiment |

The introduction's phrase "larger-label experiments" can suggest a separate study. The current Section 5 reports only the five main budgets up to 1,600. It does not report numerical results from `03b_large_label_robustness.ipynb`. If that phrase is intended to refer to additional experiments, those results and their code should be documented separately; otherwise, it can be tightened to the larger budgets within the main experiment.

## Figure source limitation and repair

The legacy plotting script is not the final rendering script for the supplied PDF. It lacks the empirical-SD bars shown in Figure 1 and plots only two methods in Figure 2, without the highlighted cross-budget comparison. Its plotting style also differs from the manuscript.

The new `scripts/make_figures.R` uses the archived CSV directly and reconstructs all three figures with the correct comparisons. It adds the SD bars, all four curves, and the highlights. It was executed successfully during repository preparation. These are new figure renderings; they are not presented as recovered original figure source files.

## Validation performed

- `python scripts/verify_results.py`: passed, including 1,200 summary metrics and 45 printed table values.
- `python scripts/export_tables.py`: produced Tables 6-8 as CSV files.
- R syntax parsing: original runner and all three R helper files passed.
- `Rscript scripts/make_figures.R`: generated three PNG and three PDF figures successfully; figures were visually inspected.
- `scripts/check_environment.R`: successfully recorded the local environment and reported missing `future` and `future.apply`. This expected failure is an environment readiness check, not a simulation result.
- The original runner's hash was checked after figure generation and remained unchanged.

Local figure-validation environment: R 4.4.3 on Windows; gbm 2.2.3, splitTools 1.0.1, dplyr 1.2.1, ggplot2 4.0.3. These are local validation versions, not recovered versions from the original cloud run.

## Limits on a fresh simulation run

- No full Monte Carlo rerun was performed during repository preparation. `future` and `future.apply` were absent from the inspected local R environment; the provided installation helper includes them.
- The supplied cloud archive does not include an R session report or a package lockfile. Seeds and matching source files alone do not establish bit-for-bit reproduction under changed R/package versions.
- The runner is self-contained apart from its R packages; it does not read an external empirical dataset or source other project scripts.
- Labels are sampled from calibrated probabilities, so label budgets are expected counts rather than fixed counts.
- Summary calculations use the population ATE of 2. Empirical SD uses sample SD across replications; RMSE averages squared errors with denominator 50. Consequently, an RMSE can be slightly smaller than the reported SD without being an error.
- Saved checkpoints are accepted by filename without checking code, tuning parameters, or package versions. Summaries aggregate all matching checkpoints, irrespective of the last requested mode.
- Within a bundle, learner fitting consumes random numbers in scenario order. Resuming only some missing lambda checkpoints may change later learner randomness compared with a fresh complete bundle. For a strict fresh-run comparison, use a new checkout with no generated checkpoints.
- The new-run output directory is intentionally separate from `data/reference/`, so new results do not silently replace the archived manuscript evidence.

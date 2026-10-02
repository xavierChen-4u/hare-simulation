#!/usr/bin/env python3
"""Verify archived Monte Carlo results and the manuscript's numerical tables.

Uses only the Python standard library; does not run R or rerun simulations.
The defaults resolve relative to this script, so it can run from any directory.
Comparison tolerance accommodates CSV serialization, not simulation randomness.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
from pathlib import Path
from statistics import NormalDist, fmean, stdev
import sys


SCENARIOS = {"lambda_0": 0.0, "lambda_0p5": 0.5, "lambda_1": 1.0}
BUDGETS = (100, 200, 400, 800, 1600)
METHODS = ("Bridge-Constant", "Kallus Oracle", "Kallus-GBM", "OI-GBM")
REPETITIONS = 50
TRUE_ATE = 2.0

# Printed values transcribed from Siyu_Chen_Writing_Sample.pdf, pp. 19, 20, 22.
# Table 6: SI bias, HARE bias, absolute-bias reduction (%).
TABLE_6 = {
    100: ("-0.5331", "0.0092", "98.28"),
    200: ("-0.5299", "-0.0548", "89.66"),
    400: ("-0.5230", "-0.0293", "94.40"),
    800: ("-0.5201", "0.0103", "98.01"),
    1600: ("-0.5280", "0.0093", "98.24"),
}
# Table 7: target-only RMSE, HARE RMSE, relative RMSE reduction (%).
TABLE_7 = {
    100: ("0.8691", "0.3593", "58.66"),
    200: ("0.7658", "0.5270", "31.19"),
    400: ("0.4014", "0.2867", "28.59"),
    800: ("0.2392", "0.1416", "40.81"),
    1600: ("0.1549", "0.1175", "24.16"),
}
# Table 8 (budget 800): HARE bias, HARE SD, HARE RMSE, target-only RMSE,
# relative RMSE reduction (%).
TABLE_8 = {
    "lambda_0": ("0.0103", "0.1426", "0.1416", "0.2392", "40.8"),
    "lambda_0p5": ("0.0153", "0.1730", "0.1719", "0.2612", "34.2"),
    "lambda_1": ("0.0203", "0.2255", "0.2241", "0.2905", "22.8"),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def number(value: str) -> float:
    return math.nan if value in ("", "NA", "NaN", "n.a.") else float(value)


def mean(rows: list[dict[str, str]], column: str) -> float:
    values = [number(row[column]) for row in rows]
    available = [value for value in values if not math.isnan(value)]
    return fmean(available) if available else math.nan


def key(row: dict[str, str]) -> tuple[str, int, str]:
    return row["scenario"], int(row["label_budget"]), row["method"]


def recompute(rows: list[dict[str, str]]) -> dict[str, float]:
    estimates = [number(row["estimate"]) for row in rows]
    standard_errors = [number(row["standard_error"]) for row in rows]
    z = NormalDist().inv_cdf(0.975)
    empirical_sd = stdev(estimates)
    result = {
        "successful_replications": sum(math.isfinite(x) for x in estimates),
        "mean_labeled_N": mean(rows, "actual_labeled_N"),
        "mean_estimate": fmean(estimates),
        "bias": fmean(x - TRUE_ATE for x in estimates),
        "empirical_sd": empirical_sd,
        "rmse": math.sqrt(fmean((x - TRUE_ATE) ** 2 for x in estimates)),
        "mean_standard_error": fmean(standard_errors),
        "coverage": fmean(
            x - z * se <= TRUE_ATE <= x + z * se
            for x, se in zip(estimates, standard_errors)
        ),
        "mean_ci_width": fmean(2 * z * se for se in standard_errors),
        "se_to_empirical_sd": fmean(standard_errors) / empirical_sd,
    }
    averaged_columns = {
        "mean_outcome_rmse": "outcome_rmse",
        "mean_correction_rmse": "correction_rmse",
        "mean_history_h_rmse": "history_h_rmse",
        "mean_delta1_hat": "mean_delta1_hat",
        "mean_delta0_hat": "mean_delta0_hat",
        "mean_target_d1_true": "target_mean_d1_true",
        "mean_target_d0_true": "target_mean_d0_true",
        "mean_target_sd_d1": "target_sd_d1_true",
        "mean_target_sd_d0": "target_sd_d0_true",
        "mean_sample_true_ATE": "sample_true_ATE",
    }
    result.update({target: mean(rows, source) for target, source in averaged_columns.items()})
    return result


def main(argv: list[str] | None = None) -> int:
    repository = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw", type=Path,
        default=repository / "data/reference/lambda_transport_results.csv",
    )
    parser.add_argument(
        "--summary", type=Path,
        default=repository / "data/reference/lambda_transport_summary.csv",
    )
    args = parser.parse_args(argv)
    failures: list[str] = []

    def check(condition: bool, label: str) -> None:
        if not condition:
            failures.append(label)

    def near(actual: float, expected: float, label: str) -> None:
        both_missing = math.isnan(actual) and math.isnan(expected)
        check(both_missing or math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-11),
              f"{label}: actual={actual!r}, expected={expected!r}")

    def printed(actual: float, expected: str, label: str) -> None:
        digits = len(expected.partition(".")[2])
        rendered = f"{actual:.{digits}f}"
        check(rendered == expected, f"{label}: calculated {rendered}, manuscript {expected}")

    try:
        raw = read_csv(args.raw)
        summaries = read_csv(args.summary)
        expected_keys = set(itertools.product(SCENARIOS, BUDGETS, METHODS))
        groups: dict[tuple[str, int, str], list[dict[str, str]]] = {}
        for row in raw:
            groups.setdefault(key(row), []).append(row)
        summary_by_key = {key(row): row for row in summaries}
        check(len(raw) == 3000, f"Expected 3,000 raw rows, found {len(raw)}")
        check(len(summaries) == 60, f"Expected 60 summary rows, found {len(summaries)}")
        check(len(summary_by_key) == len(summaries), "Duplicate summary keys")
        check(set(groups) == expected_keys, "Unexpected or missing raw scenario/budget/method keys")
        check(set(summary_by_key) == expected_keys,
              "Unexpected or missing summary scenario/budget/method keys")
        for group_key, rows in groups.items():
            check(len(rows) == REPETITIONS, f"{group_key}: expected 50 rows")
            check({int(r["replication"]) for r in rows} == set(range(1, REPETITIONS + 1)),
                  f"{group_key}: repetition IDs must be exactly 1..50")
        for index, row in enumerate(raw, 2):
            check(math.isfinite(number(row["estimate"])), f"raw line {index}: nonfinite estimate")
            se = number(row["standard_error"])
            check(math.isfinite(se) and se >= 0, f"raw line {index}: invalid standard error")
            check(number(row["true_ATE"]) == TRUE_ATE, f"raw line {index}: true ATE differs from 2")
            check(number(row["lambda"]) == SCENARIOS.get(row["scenario"]),
                  f"raw line {index}: scenario/lambda mismatch")
            check(number(row["history_best_trees"]) == 3000,
                  f"raw line {index}: historical model did not use 3,000 trees")
        if failures:
            for failure in failures:
                print(f"FAIL: {failure}", file=sys.stderr)
            return 1

        metric_count = 0
        for group_key, rows in groups.items():
            summary = summary_by_key[group_key]
            near(number(summary["lambda"]), SCENARIOS[group_key[0]], f"{group_key}: lambda")
            for column, value in recompute(rows).items():
                near(value, number(summary[column]), f"{group_key}: {column}")
                metric_count += 1

        def metric(scenario: str, budget: int, method: str, column: str) -> float:
            return number(summary_by_key[(scenario, budget, method)][column])

        table_cells = 0
        for budget, expected in TABLE_6.items():
            si = metric("lambda_0", budget, "OI-GBM", "bias")
            hare = metric("lambda_0", budget, "Bridge-Constant", "bias")
            actual = (si, hare, 100 * (abs(si) - abs(hare)) / abs(si))
            for column, (value, display) in enumerate(zip(actual, expected), 1):
                printed(value, display, f"Table 6 budget {budget} column {column}")
                table_cells += 1
        for budget, expected in TABLE_7.items():
            sa = metric("lambda_0", budget, "Kallus-GBM", "rmse")
            hare = metric("lambda_0", budget, "Bridge-Constant", "rmse")
            actual = (sa, hare, 100 * (1 - hare / sa))
            for column, (value, display) in enumerate(zip(actual, expected), 1):
                printed(value, display, f"Table 7 budget {budget} column {column}")
                table_cells += 1
        for scenario, expected in TABLE_8.items():
            hare = metric(scenario, 800, "Bridge-Constant", "rmse")
            sa = metric(scenario, 800, "Kallus-GBM", "rmse")
            actual = (
                metric(scenario, 800, "Bridge-Constant", "bias"),
                metric(scenario, 800, "Bridge-Constant", "empirical_sd"),
                hare, sa, 100 * (1 - hare / sa),
            )
            for column, (value, display) in enumerate(zip(actual, expected), 1):
                printed(value, display, f"Table 8 {scenario} column {column}")
                table_cells += 1

        additional_claims = (
            ("lambda_0", 800, "Bridge-Constant", "mean_labeled_N", "806.6"),
            ("lambda_0", 1600, "Kallus-GBM", "mean_labeled_N", "1598.54"),
            ("lambda_0p5", 1600, "Kallus-GBM", "rmse", "0.1688"),
            ("lambda_1", 1600, "Kallus-GBM", "rmse", "0.1854"),
            ("lambda_1", 100, "OI-GBM", "bias", "-0.5331"),
            ("lambda_1", 100, "Bridge-Constant", "bias", "0.0372"),
            ("lambda_1", 100, "Bridge-Constant", "rmse", "0.6709"),
            ("lambda_1", 100, "OI-GBM", "rmse", "0.5353"),
        )
        for scenario, budget, method, column, display in additional_claims:
            printed(metric(scenario, budget, method, column), display,
                    f"Text claim {scenario}, {budget}, {method}, {column}")
        for scenario, budget in itertools.product(SCENARIOS, BUDGETS):
            check(metric(scenario, budget, "Bridge-Constant", "rmse") <
                  metric(scenario, budget, "Kallus-GBM", "rmse"),
                  f"HARE RMSE is not below target-only RMSE at {scenario}, {budget}")
        check(metric("lambda_0", 800, "Bridge-Constant", "rmse") <
              metric("lambda_0", 1600, "Kallus-GBM", "rmse"),
              "Constant-discrepancy cross-budget label-efficiency comparison failed")
        for scenario in ("lambda_0p5", "lambda_1"):
            check(metric(scenario, 800, "Bridge-Constant", "rmse") >
                  metric(scenario, 1600, "Kallus-GBM", "rmse"),
                  f"Heterogeneous cross-budget comparison failed at {scenario}")
        if failures:
            for failure in failures:
                print(f"FAIL: {failure}", file=sys.stderr)
            return 1
        print("PASS: 3,000 raw rows; 3 scenarios x 5 budgets x 4 methods x 50 repetitions.")
        print(f"PASS: all {metric_count} summary metrics recomputed from archived replication data.")
        print(f"PASS: all {table_cells} numerical cells in manuscript Tables 6-8 match printed rounding.")
        print(f"PASS: {len(additional_claims)} additional numerical text claims and all RMSE orderings.")
        print("This verifies archived numerical consistency; it does not rerun simulations or verify figure rendering.")
        return 0
    except (OSError, ValueError, KeyError, ZeroDivisionError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

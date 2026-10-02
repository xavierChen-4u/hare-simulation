#!/usr/bin/env python3
"""Export manuscript Tables 6-8 from the archived summary (Python stdlib only)."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    with (ROOT / "data/reference/lambda_transport_summary.csv").open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    lookup = {(float(r["lambda"]), int(r["label_budget"]), r["method"]): r for r in rows}
    budgets = (100, 200, 400, 800, 1600)
    table6, table7, table8 = [], [], []
    for budget in budgets:
        si = lookup[0, budget, "OI-GBM"]
        hare = lookup[0, budget, "Bridge-Constant"]
        sa = lookup[0, budget, "Kallus-GBM"]
        table6.append(dict(label_budget=budget, SI_AIPW_GBM_bias=float(si["bias"]), HARE_C_bias=float(hare["bias"]), absolute_bias_reduction_pct=100 * (1 - abs(float(hare["bias"])) / abs(float(si["bias"])))))
        table7.append(dict(label_budget=budget, SA_EIF_GBM_RMSE=float(sa["rmse"]), HARE_C_RMSE=float(hare["rmse"]), RMSE_reduction_pct=100 * (1 - float(hare["rmse"]) / float(sa["rmse"]))))
    for lam in (0, 0.5, 1):
        hare = lookup[lam, 800, "Bridge-Constant"]
        sa = lookup[lam, 800, "Kallus-GBM"]
        table8.append(dict(lambda_value=lam, HARE_C_bias=float(hare["bias"]), HARE_C_empirical_SD=float(hare["empirical_sd"]), HARE_C_RMSE=float(hare["rmse"]), SA_EIF_GBM_RMSE=float(sa["rmse"]), RMSE_reduction_pct=100 * (1 - float(hare["rmse"]) / float(sa["rmse"]))))
    (ROOT / "tables").mkdir(exist_ok=True)
    for name, data in (("table6_bias.csv", table6), ("table7_rmse.csv", table7), ("table8_heterogeneity.csv", table8)):
        with (ROOT / "tables" / name).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
        print("Wrote tables/" + name)


if __name__ == "__main__":
    main()

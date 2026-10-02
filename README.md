# HARE: History-Assisted Residual Estimation

**Estimating long-term treatment effects with limited target outcomes**

Research code by Siyu Chen, Renmin University of China.

HARE estimates long-term average treatment effects (ATEs) when short-term measures are widely observed, only a subset of target outcomes is available, and an independent historical sample contains both short- and long-term outcomes. Historical data help learn the outcome relationship; observed target outcomes correct for differences between the historical and target populations.

This repository implements **HARE-C**, a version of HARE with a constant residual correction in each treatment arm, and evaluates its accuracy and target-label efficiency in Monte Carlo simulations.

## Method

Let $X$ denote pretreatment covariates, $S$ short-term measures, and $Y$ the long-term outcome. HARE decomposes the target outcome regression in treatment arm $t$ into a historical predictor and a target-specific discrepancy:

$$
m_t^E(x,s) = h_H(x,s) + d_t(x,s).
$$

When the discrepancy is simpler to learn than the full outcome relationship, a small labeled target sample can be used primarily to estimate that discrepancy. HARE-C approximates it by one constant per treatment arm:

$$
\widehat m_t^E(x,s) = \widehat h_H(x,s) + \widehat\Delta_t.
$$

The implementation proceeds in three steps:

1. Fit the historical outcome predictor $\widehat h_H(X,S)$ using gradient-boosted trees.
2. Within each target training fold, estimate $\widehat\Delta_t$ as the mean of $Y-\widehat h_H(X,S)$ among labeled observations in arm $t$.
3. Combine the corrected predictions with the surrogate-assisted estimating equation of Kallus and Mao. Five-fold cross-fitting separates target model fitting from evaluation. The final estimator retains weighted residual corrections from observed target outcomes and uses self-normalized weights.

The method allows the historical and target outcome relationships to differ and does not require statistical surrogacy or cross-population comparability. Identification still requires consistency, unconfounded treatment assignment, missing-at-random target outcomes, and positivity.

## Simulation results

The main design has 8,000 historical observations and 8,000 target observations, with a true average treatment effect of 2. Under constant historical-target discrepancies and an expected target-label budget of 800, empirical bias and root mean squared error (RMSE) over 50 replications are:

| Method | Empirical bias | ATE RMSE |
| --- | ---: | ---: |
| HARE-C | 0.0103 | 0.1416 |
| SA-EIF-GBM, target-only | 0.0317 | 0.2392 |
| SI-AIPW-GBM, historical index | -0.5201 | 0.5230 |
| SA-EIF-Oracle | -0.0114 | 0.1181 |

HARE-C has **40.8% lower RMSE** than the target-only estimator at the same label budget. Its RMSE at 800 target labels (0.1416) is also below the target-only estimator's RMSE at 1,600 labels (0.1549). This comparison measures savings in target outcome collection, conditional on the availability of the historical sample.

![ATE RMSE across expected target-label budgets under constant discrepancies](figures/rmse_label_budget_curve.png)

*The outlined points compare HARE-C at a budget of 800 with SA-EIF-GBM at 1,600. Each point summarizes 50 Monte Carlo replications; label budgets are expected counts.*

The advantage depends on the discrepancy structure. At a budget of 800, increasing discrepancy heterogeneity from $\lambda=0$ to $\lambda=1$ reduces HARE-C's relative RMSE improvement over the target-only estimator from 40.8% to 22.8%. With an expected budget of only 100 target labels and $\lambda=1$, HARE-C reduces historical-index bias but has higher RMSE than the historical-index estimator. These experiments evaluate the benefits and limits of the constant correction model.

## Experimental design

Target outcomes follow

$$
Y(t) = h_H(X,S(t)) + \Delta_t + \lambda u_t(X,S(t)) + \varepsilon_t,
$$

with $\Delta_1=0.25$ and $\Delta_0=-0.25$. The nonlinear discrepancy functions $u_t$ are centered and approximately scaled to unit standard deviation under the corresponding potential-surrogate distributions. Varying $\lambda$ changes heterogeneity while preserving the population ATE and the historical-index effect.

| Component | Setting |
| --- | --- |
| Historical / target sample size | 8,000 / 8,000 |
| Pretreatment covariates / short-term measures | 6 / 5 |
| Expected target-label budgets | 100, 200, 400, 800, 1,600 |
| Discrepancy heterogeneity, $\lambda$ | 0, 0.5, 1 |
| Replications | 50 per scenario-budget configuration |
| Target cross-fitting | 5 folds, stratified by treatment and label observability |
| Outcome learner | GBM; at most 3,000 trees, depth 3, learning rate 0.03 |
| Population ATE / historical-index effect | 2 / 1.5 |

The feasible estimators use the same target observations and label assignments within each replication. The comparisons separate the contribution of historical prediction from that of observed target outcomes:

| Method | Estimation approach | Identifier in the result files |
| --- | --- | --- |
| HARE-C | Historical prediction with arm-specific target residual corrections | `Bridge-Constant` |
| SA-EIF-GBM | Surrogate-assisted estimation using target data only | `Kallus-GBM` |
| SI-AIPW-GBM | AIPW estimation of the effect on the historical prediction index | `OI-GBM` |
| SA-EIF-Oracle | Reference estimator supplied with the true nuisance functions | `Kallus Oracle` |

## Run the code

Use R 4.1 or later with compatible package versions. From the repository root, install dependencies and run a small pilot:

```sh
Rscript scripts/install_dependencies.R
Rscript scripts/check_environment.R
Rscript 04_cloud_run.R pilot 1
```

The pilot runs one replication at an expected label budget of 200 for all three discrepancy settings. To run the full experiment and plot the resulting summaries:

```sh
Rscript 04_cloud_run.R full 1
Rscript 04_cloud_run.R summary
Rscript scripts/make_figures.R 04_lambda_transport_cloud/results/lambda_transport_summary.csv
```

The final argument to the simulation runner is the worker count. Increase it to use parallel processing; omitting it defaults to 24 workers. A `fast` mode runs five replications per configuration. Results and checkpoints are saved under `04_lambda_transport_cloud/results/`. Use a fresh output directory when changing the simulation code or environment.

### Explore the reported results

The repository includes the 3,000 method-level replication records and their summaries in `data/reference/`. Generate the three result figures without rerunning the simulation:

```sh
Rscript scripts/make_figures.R
```

Python 3.10 or later, with no additional packages, can regenerate Tables 6-8 and check the summaries against the replication-level results:

```sh
python scripts/export_tables.py
python scripts/verify_results.py
```

## Repository contents

```text
04_cloud_run.R       Simulation design, estimators, and parallel runner
data/reference/     Replication-level results and summary statistics
scripts/            Dependency setup and figure/table generation
figures/            Bias, label-efficiency, and heterogeneity figures
tables/             Numerical comparisons from the simulation study
docs/               Design-to-code mapping and reproducibility notes
```

Detailed source mappings and environment information are in [the reproducibility notes](docs/PROVENANCE.md). The original cloud package versions were not recorded, so fresh runs may differ numerically from the included results.

## Associated paper

Siyu Chen. *HARE: History-Assisted Residual Estimation of Long-Term Treatment Effects with Limited Target Outcomes.* Research writing sample, Renmin University of China.

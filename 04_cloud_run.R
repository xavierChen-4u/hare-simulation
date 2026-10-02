#!/usr/bin/env Rscript

# Cloud runner for the heterogeneous transport simulation.
#
# Usage:
#   Rscript 04_cloud_run.R pilot 1
#   Rscript 04_cloud_run.R fast  12
#   Rscript 04_cloud_run.R full  24
#   Rscript 04_cloud_run.R summary
#
# The three scenarios differ only in lambda:
#   d_t(X,S) = delta_t + lambda * u_t(X,S),
# where u_t is a centered, unit-SD nonlinear function as complex as the
# shared historical outcome surface. Thus lambda is also the SD of the
# heterogeneous part of the transport shift.

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

options(
  warn = 1,
  future.globals.maxSize = 4 * 1024^3
)

required_packages <- c(
  "gbm",
  "splitTools",
  "dplyr",
  "future",
  "future.apply"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

script_directory <- function() {
  file_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  if (length(file_argument) == 0L) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }

  normalizePath(
    dirname(sub("^--file=", "", file_argument[[1L]])),
    mustWork = TRUE
  )
}

BASE_DIRECTORY <- script_directory()

PROJECT_DIRECTORY <- file.path(
  BASE_DIRECTORY,
  "04_lambda_transport_cloud"
)

RESULT_DIRECTORY <- file.path(
  PROJECT_DIRECTORY,
  "results"
)

CHECKPOINT_DIRECTORY <- file.path(
  RESULT_DIRECTORY,
  "checkpoints_lambda_v1"
)

ERROR_DIRECTORY <- file.path(
  RESULT_DIRECTORY,
  "errors_lambda_v1"
)

dir.create(
  CHECKPOINT_DIRECTORY,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ERROR_DIRECTORY,
  recursive = TRUE,
  showWarnings = FALSE
)

BASE_SEED <- 20261010L

N_E <- 8000L
N_H <- 8000L

LABEL_BUDGETS <- c(
  100L,
  200L,
  400L,
  800L,
  1600L
)

SCENARIOS <- c(
  "lambda_0",
  "lambda_0p5",
  "lambda_1"
)

LAMBDA_VALUES <- c(
  lambda_0 = 0,
  lambda_0p5 = 0.5,
  lambda_1 = 1
)

TRUE_TAU <- 2
DELTA1_MEAN <- 0.25
DELTA0_MEAN <- -0.25

K_FOLDS <- 5L

GBM_MAX_TREES <- 3000L
GBM_DEPTH <- 3L
GBM_SHRINKAGE <- 0.03
GBM_BAG_FRACTION <- 0.8
GBM_MIN_NODE <- 5L
GBM_MIN_SELECTED_TREES <- 50L

x_terms <- paste0("X", 1:6)
s_terms <- paste0("S", 1:5)

x_square_terms <- paste0(
  "I(",
  x_terms,
  "^2)"
)

NONLINEAR_FORMULAS <- list(
  label = as.formula(
    paste(
      "R ~",
      paste(
        c(x_terms, x_square_terms),
        collapse = " + "
      ),
      "+ offset(offset_pi)"
    )
  ),

  treatment = as.formula(
    paste(
      "T ~ 0 +",
      paste(x_terms, collapse = " + ")
    )
  ),

  target_outcome = as.formula(
    paste(
      "Y_obs ~",
      paste(
        c(x_terms, x_square_terms, s_terms),
        collapse = " + "
      )
    )
  ),

  target_baseline = as.formula(
    paste(
      "Y_obs ~",
      paste(
        c(x_terms, x_square_terms),
        collapse = " + "
      )
    )
  ),

  history_outcome = as.formula(
    paste(
      "Y ~",
      paste(
        c(x_terms, x_square_terms, s_terms),
        collapse = " + "
      )
    )
  ),

  history_baseline = as.formula(
    paste(
      "pseudo_h ~",
      paste(
        c(x_terms, x_square_terms),
        collapse = " + "
      )
    )
  )
)

clip_probability <- function(
    p,
    lower = 0.05,
    upper = 0.95
) {
  pmin(pmax(p, lower), upper)
}

calibrate_label_intercept <- function(
    g_value,
    target_rate
) {
  objective <- function(intercept) {
    mean(plogis(intercept + g_value)) - target_rate
  }

  uniroot(
    objective,
    interval = c(-30, 30)
  )$root
}

calculate_label_shift <- function(X) {
  (
    0.40 * X[, 1] -
      0.40 * X[, 2] +
      0.25 * X[, 3]^2 -
      0.25 * X[, 4]^2
  )
}

calculate_treatment_probability <- function(X) {
  eta <- c(
    1,
    -0.5,
    -0.5,
    -0.5,
    -0.5,
    1
  )

  clip_probability(
    plogis(drop(X %*% eta)),
    lower = 0.05,
    upper = 0.95
  )
}

EXPECTED_TANH_SUM_TREATED <- integrate(
  f = function(z) {
    tanh(z) * dnorm(
      z,
      mean = 2,
      sd = sqrt(2)
    )
  },
  lower = -Inf,
  upper = Inf,
  rel.tol = 1e-10
)$value

TRUE_DIRECT_EFFECT <- DELTA1_MEAN - DELTA0_MEAN
TARGET_SURROGATE_EFFECT <- TRUE_TAU - TRUE_DIRECT_EFFECT

TANH_COEFFICIENT <- (
  TARGET_SURROGATE_EFFECT /
    (2 * EXPECTED_TANH_SUM_TREATED)
)

calculate_complex_h <- function(X, S) {
  (
    1.2 * sin(X[, 1]) +
      0.8 * (X[, 2]^2 - 1) +
      0.6 * X[, 3] * X[, 4] +
      TANH_COEFFICIENT * tanh(S[, 1] + S[, 2]) +
      0.5 * (S[, 3]^2 - 1) +
      0.6 * X[, 1] * S[, 1]
  )
}

calculate_true_mu_h <- function(X, treatment) {
  x_part <- (
    1.2 * sin(X[, 1]) +
      0.8 * (X[, 2]^2 - 1) +
      0.6 * X[, 3] * X[, 4]
  )

  if (treatment == 1) {
    return(
      x_part +
        TANH_COEFFICIENT * EXPECTED_TANH_SUM_TREATED +
        0.5 +
        0.6 * X[, 1]
    )
  }

  x_part -
    TANH_COEFFICIENT * EXPECTED_TANH_SUM_TREATED +
    0.5 -
    0.6 * X[, 1]
}

# This is deliberately as nonlinear as the shared historical surface.
# It is centered under S(1) for treatment 1 and S(0) for treatment 0.
calculate_raw_transport_component <- function(
    X,
    S,
    treatment
) {
  if (treatment == 1) {
    return(
      1.00 * sin(X[, 1]) +
        0.80 * (X[, 2]^2 - 1) +
        0.60 * X[, 3] * X[, 4] +
        0.80 * (
          tanh(S[, 1] + S[, 2]) -
            EXPECTED_TANH_SUM_TREATED
        ) +
        0.50 * (S[, 3]^2 - 2) +
        0.60 * X[, 1] * S[, 1]
    )
  }

  -1.10 * sin(X[, 1]) +
    0.75 * (X[, 2]^2 - 1) -
    0.65 * X[, 3] * X[, 4] +
    0.85 * (
      tanh(S[, 1] + S[, 2]) +
        EXPECTED_TANH_SUM_TREATED
    ) +
    0.55 * (S[, 4]^2 - 2) -
    0.60 * X[, 1] * S[, 2]
}

SHIFT_CALIBRATION_N <- 200000L

set.seed(BASE_SEED + 900L)

X_calibration <- matrix(
  rnorm(SHIFT_CALIBRATION_N * 6L),
  nrow = SHIFT_CALIBRATION_N,
  ncol = 6L
)

S0_calibration <- matrix(
  rnorm(
    SHIFT_CALIBRATION_N * 5L,
    mean = -1,
    sd = 1
  ),
  nrow = SHIFT_CALIBRATION_N,
  ncol = 5L
)

S1_calibration <- matrix(
  rnorm(
    SHIFT_CALIBRATION_N * 5L,
    mean = 1,
    sd = 1
  ),
  nrow = SHIFT_CALIBRATION_N,
  ncol = 5L
)

raw_shift_1 <- calculate_raw_transport_component(
  X_calibration,
  S1_calibration,
  treatment = 1
)

raw_shift_0 <- calculate_raw_transport_component(
  X_calibration,
  S0_calibration,
  treatment = 0
)

SHIFT_SCALE <- c(
  treated = sd(raw_shift_1),
  control = sd(raw_shift_0)
)

SHIFT_CALIBRATION_DIAGNOSTIC <- data.frame(
  arm = c("treated", "control"),
  raw_mean = c(mean(raw_shift_1), mean(raw_shift_0)),
  raw_sd = c(sd(raw_shift_1), sd(raw_shift_0))
)

rm(
  X_calibration,
  S0_calibration,
  S1_calibration,
  raw_shift_1,
  raw_shift_0
)

invisible(gc(full = TRUE))

scenario_lambda <- function(scenario) {
  scenario <- match.arg(scenario, SCENARIOS)
  unname(LAMBDA_VALUES[[scenario]])
}

calculate_transport_shift <- function(
    X,
    S,
    treatment,
    scenario
) {
  lambda <- scenario_lambda(scenario)

  raw_component <- calculate_raw_transport_component(
    X = X,
    S = S,
    treatment = treatment
  )

  scale_value <- if (treatment == 1) {
    SHIFT_SCALE[["treated"]]
  } else {
    SHIFT_SCALE[["control"]]
  }

  mean_shift <- if (treatment == 1) {
    DELTA1_MEAN
  } else {
    DELTA0_MEAN
  }

  mean_shift + lambda * raw_component / scale_value
}

calculate_transport_shift_mu <- function(
    X,
    treatment,
    scenario
) {
  lambda <- scenario_lambda(scenario)

  raw_conditional_mean <- if (treatment == 1) {
    (
      1.00 * sin(X[, 1]) +
        0.80 * (X[, 2]^2 - 1) +
        0.60 * X[, 3] * X[, 4] +
        0.60 * X[, 1]
    )
  } else {
    (
      -1.10 * sin(X[, 1]) +
        0.75 * (X[, 2]^2 - 1) -
        0.65 * X[, 3] * X[, 4] +
        0.60 * X[, 1]
    )
  }

  scale_value <- if (treatment == 1) {
    SHIFT_SCALE[["treated"]]
  } else {
    SHIFT_SCALE[["control"]]
  }

  mean_shift <- if (treatment == 1) {
    DELTA1_MEAN
  } else {
    DELTA0_MEAN
  }

  mean_shift + lambda * raw_conditional_mean / scale_value
}

simulate_heterogeneous_target <- function(
    n_E,
    label_budget,
    scenario,
    seed = NULL
) {
  scenario <- match.arg(scenario, SCENARIOS)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  p_x <- 6L
  p_s <- 5L

  X <- matrix(
    rnorm(n_E * p_x),
    nrow = n_E,
    ncol = p_x
  )

  colnames(X) <- paste0("X", 1:p_x)

  target_label_rate <- label_budget / n_E
  g_label <- calculate_label_shift(X)

  label_intercept <- calibrate_label_intercept(
    g_value = g_label,
    target_rate = target_label_rate
  )

  r_true <- plogis(label_intercept + g_label)

  R <- rbinom(
    n_E,
    size = 1,
    prob = r_true
  )

  e_true <- calculate_treatment_probability(X)

  T <- rbinom(
    n_E,
    size = 1,
    prob = e_true
  )

  S0 <- matrix(
    rnorm(
      n_E * p_s,
      mean = -1,
      sd = 1
    ),
    nrow = n_E,
    ncol = p_s
  )

  S1 <- matrix(
    rnorm(
      n_E * p_s,
      mean = 1,
      sd = 1
    ),
    nrow = n_E,
    ncol = p_s
  )

  colnames(S0) <- paste0("S", 1:p_s)
  colnames(S1) <- paste0("S", 1:p_s)

  S <- S0
  S[T == 1, ] <- S1[T == 1, ]
  colnames(S) <- paste0("S", 1:p_s)

  h0 <- calculate_complex_h(X, S0)
  h1 <- calculate_complex_h(X, S1)
  h_observed <- calculate_complex_h(X, S)

  d0_potential <- calculate_transport_shift(
    X,
    S0,
    treatment = 0,
    scenario = scenario
  )

  d1_potential <- calculate_transport_shift(
    X,
    S1,
    treatment = 1,
    scenario = scenario
  )

  d0_observed <- calculate_transport_shift(
    X,
    S,
    treatment = 0,
    scenario = scenario
  )

  d1_observed <- calculate_transport_shift(
    X,
    S,
    treatment = 1,
    scenario = scenario
  )

  Y0 <- h0 + d0_potential + rnorm(n_E)
  Y1 <- h1 + d1_potential + rnorm(n_E)

  Y_full <- ifelse(T == 1, Y1, Y0)
  Y_obs <- ifelse(R == 1, Y_full, NA_real_)

  m1_true <- h_observed + d1_observed
  m0_true <- h_observed + d0_observed

  mu1_true <- (
    calculate_true_mu_h(X, treatment = 1) +
      calculate_transport_shift_mu(
        X,
        treatment = 1,
        scenario = scenario
      )
  )

  mu0_true <- (
    calculate_true_mu_h(X, treatment = 0) +
      calculate_transport_shift_mu(
        X,
        treatment = 0,
        scenario = scenario
      )
  )

  data.frame(
    id = seq_len(n_E),
    scenario = scenario,
    lambda = scenario_lambda(scenario),
    R = R,
    T = T,
    target_label_rate = target_label_rate,
    label_intercept = label_intercept,
    r_true = r_true,
    e_true = e_true,
    Y_full = Y_full,
    Y_obs = Y_obs,
    Y0 = Y0,
    Y1 = Y1,
    h0_true = h0,
    h1_true = h1,
    h_observed_true = h_observed,
    d0_potential_true = d0_potential,
    d1_potential_true = d1_potential,
    d0_observed_true = d0_observed,
    d1_observed_true = d1_observed,
    m1_true = m1_true,
    m0_true = m0_true,
    mu1_true = mu1_true,
    mu0_true = mu0_true,
    X,
    S,
    check.names = FALSE
  )
}

simulate_heterogeneous_history <- function(
    n_H,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  p_x <- 6L
  p_s <- 5L

  X <- matrix(
    rnorm(n_H * p_x),
    nrow = n_H,
    ncol = p_x
  )

  colnames(X) <- paste0("X", 1:p_x)

  latent_probability <- calculate_treatment_probability(X)

  latent_assignment <- rbinom(
    n_H,
    size = 1,
    prob = latent_probability
  )

  S0 <- matrix(
    rnorm(
      n_H * p_s,
      mean = -1,
      sd = 1
    ),
    nrow = n_H,
    ncol = p_s
  )

  S1 <- matrix(
    rnorm(
      n_H * p_s,
      mean = 1,
      sd = 1
    ),
    nrow = n_H,
    ncol = p_s
  )

  S <- S0
  S[latent_assignment == 1, ] <- S1[latent_assignment == 1, ]
  colnames(S) <- paste0("S", 1:p_s)

  h_value <- calculate_complex_h(X, S)
  Y <- h_value + rnorm(n_H)

  data.frame(
    id_H = seq_len(n_H),
    Y = Y,
    X,
    S,
    check.names = FALSE
  )
}

fit_gbm_learner <- function(
    formula,
    data,
    distribution = "gaussian",
    tree_selection = c("oob", "fixed"),
    max_trees = GBM_MAX_TREES
) {
  tree_selection <- match.arg(tree_selection)

  fitted_model <- gbm::gbm(
    formula = formula,
    data = data,
    distribution = distribution,
    n.trees = max_trees,
    interaction.depth = GBM_DEPTH,
    shrinkage = GBM_SHRINKAGE,
    bag.fraction = GBM_BAG_FRACTION,
    n.minobsinnode = GBM_MIN_NODE,
    train.fraction = 1,
    keep.data = FALSE,
    verbose = FALSE
  )

  if (tree_selection == "fixed") {
    best_trees <- max_trees
  } else {
    best_trees <- suppressMessages(
      suppressWarnings(
        gbm::gbm.perf(
          fitted_model,
          method = "OOB",
          plot.it = FALSE
        )
      )
    )

    if (
      length(best_trees) == 0L ||
        !is.finite(best_trees)
    ) {
      best_trees <- max_trees
    }

    best_trees <- min(
      max_trees,
      max(
        GBM_MIN_SELECTED_TREES,
        as.integer(best_trees)
      )
    )
  }

  list(
    model = fitted_model,
    best_trees = best_trees
  )
}

predict_gbm_learner <- function(
    fitted_object,
    newdata
) {
  predict(
    fitted_object$model,
    newdata = newdata,
    n.trees = fitted_object$best_trees,
    type = "response"
  )
}

fit_predict_release_gbm <- function(
    formula,
    training_data,
    prediction_data,
    tree_selection = "oob",
    max_trees = GBM_MAX_TREES
) {
  fitted_object <- fit_gbm_learner(
    formula = formula,
    data = training_data,
    tree_selection = tree_selection,
    max_trees = max_trees
  )

  prediction <- predict_gbm_learner(
    fitted_object,
    prediction_data
  )

  rm(fitted_object)
  invisible(gc(full = TRUE))

  prediction
}

create_empty_frame <- function(n, column_names) {
  output <- as.data.frame(
    matrix(
      NA_real_,
      nrow = n,
      ncol = length(column_names)
    )
  )

  names(output) <- column_names
  output
}

# Shared across lambda values: r, e, historical h, and OI q_t(X).
fit_predict_shared_fold <- function(
    train_data,
    test_data,
    history_model
) {
  fit_r <- glm(
    NONLINEAR_FORMULAS$label,
    data = train_data,
    family = binomial()
  )

  r_hat <- clip_probability(
    predict(
      fit_r,
      newdata = test_data,
      type = "response"
    ),
    lower = 1e-6,
    upper = 1
  )

  rm(fit_r)
  invisible(gc(full = TRUE))

  fit_e <- glm(
    NONLINEAR_FORMULAS$treatment,
    data = train_data,
    family = binomial()
  )

  e_hat <- clip_probability(
    predict(
      fit_e,
      newdata = test_data,
      type = "response"
    ),
    lower = 0.05,
    upper = 0.95
  )

  rm(fit_e)
  invisible(gc(full = TRUE))

  h_train <- predict_gbm_learner(
    history_model,
    train_data
  )

  h_test <- predict_gbm_learner(
    history_model,
    test_data
  )

  train_augmented <- train_data
  train_augmented$pseudo_h <- h_train

  q1_hat <- fit_predict_release_gbm(
    NONLINEAR_FORMULAS$history_baseline,
    train_augmented[
      train_augmented$T == 1,
      ,
      drop = FALSE
    ],
    test_data
  )

  q0_hat <- fit_predict_release_gbm(
    NONLINEAR_FORMULAS$history_baseline,
    train_augmented[
      train_augmented$T == 0,
      ,
      drop = FALSE
    ],
    test_data
  )

  rm(train_augmented)
  invisible(gc(full = TRUE))

  list(
    h_train = h_train,
    test = data.frame(
      e = e_hat,
      r = r_hat,
      h = h_test,
      q1 = q1_hat,
      q0 = q0_hat
    )
  )
}

# Scenario-specific fits: Kallus learns the entire outcome surface, while
# Bridge learns only two arm-specific constants on top of historical h.
fit_predict_scenario_fold <- function(
    train_data,
    test_data,
    shared_fold
) {
  labeled_treated <- (
    train_data$R == 1 &
      train_data$T == 1
  )

  labeled_control <- (
    train_data$R == 1 &
      train_data$T == 0
  )

  if (
    sum(labeled_treated) < 20L ||
      sum(labeled_control) < 20L
  ) {
    stop("Too few labeled observations in a treatment arm.")
  }

  treated_labeled_data <- train_data[
    labeled_treated,
    ,
    drop = FALSE
  ]

  control_labeled_data <- train_data[
    labeled_control,
    ,
    drop = FALSE
  ]

  kallus_m1_hat <- fit_predict_release_gbm(
    NONLINEAR_FORMULAS$target_outcome,
    treated_labeled_data,
    test_data
  )

  kallus_m0_hat <- fit_predict_release_gbm(
    NONLINEAR_FORMULAS$target_outcome,
    control_labeled_data,
    test_data
  )

  kallus_mu1_hat <- fit_predict_release_gbm(
    NONLINEAR_FORMULAS$target_baseline,
    treated_labeled_data,
    test_data
  )

  kallus_mu0_hat <- fit_predict_release_gbm(
    NONLINEAR_FORMULAS$target_baseline,
    control_labeled_data,
    test_data
  )

  history_residual <- (
    train_data$Y_obs -
      shared_fold$h_train
  )

  delta1_constant <- mean(
    history_residual[labeled_treated],
    na.rm = TRUE
  )

  delta0_constant <- mean(
    history_residual[labeled_control],
    na.rm = TRUE
  )

  n_test <- nrow(test_data)

  list(
    kallus = data.frame(
      e = shared_fold$test$e,
      r = shared_fold$test$r,
      m1 = kallus_m1_hat,
      m0 = kallus_m0_hat,
      mu1 = kallus_mu1_hat,
      mu0 = kallus_mu0_hat
    ),

    bridge_constant = data.frame(
      e = shared_fold$test$e,
      r = shared_fold$test$r,
      m1 = shared_fold$test$h + delta1_constant,
      m0 = shared_fold$test$h + delta0_constant,
      mu1 = shared_fold$test$q1 + delta1_constant,
      mu0 = shared_fold$test$q0 + delta0_constant,
      d1 = rep(delta1_constant, n_test),
      d0 = rep(delta0_constant, n_test)
    )
  )
}

crossfit_lambda_bundle <- function(
    E_data_list,
    history_model,
    K = 5L,
    seed = 1L
) {
  scenario_names <- names(E_data_list)

  for (scenario in scenario_names) {
    E_data_list[[scenario]]$offset_pi <- qlogis(
      mean(E_data_list[[scenario]]$R)
    )
  }

  base_data <- E_data_list[[1L]]

  strata <- interaction(
    base_data$R,
    base_data$T,
    drop = TRUE
  )

  test_folds <- splitTools::create_folds(
    y = strata,
    k = K,
    type = "stratified",
    invert = TRUE,
    shuffle = TRUE,
    seed = seed
  )

  shared_predictions <- create_empty_frame(
    nrow(base_data),
    c("e", "r", "h", "q1", "q0")
  )

  scenario_predictions <- lapply(
    scenario_names,
    function(unused_scenario) {
      list(
        kallus = create_empty_frame(
          nrow(base_data),
          c("e", "r", "m1", "m0", "mu1", "mu0")
        ),

        bridge_constant = create_empty_frame(
          nrow(base_data),
          c(
            "e", "r", "m1", "m0", "mu1", "mu0", "d1", "d0"
          )
        )
      )
    }
  )

  names(scenario_predictions) <- scenario_names
  all_indices <- seq_len(nrow(base_data))

  for (fold_id in seq_along(test_folds)) {
    test_indices <- test_folds[[fold_id]]
    train_indices <- setdiff(all_indices, test_indices)

    shared_fold <- fit_predict_shared_fold(
      train_data = base_data[
        train_indices,
        ,
        drop = FALSE
      ],
      test_data = base_data[
        test_indices,
        ,
        drop = FALSE
      ],
      history_model = history_model
    )

    shared_predictions[test_indices, ] <- shared_fold$test

    for (scenario in scenario_names) {
      scenario_data <- E_data_list[[scenario]]

      scenario_fold <- fit_predict_scenario_fold(
        train_data = scenario_data[
          train_indices,
          ,
          drop = FALSE
        ],
        test_data = scenario_data[
          test_indices,
          ,
          drop = FALSE
        ],
        shared_fold = shared_fold
      )

      scenario_predictions[[scenario]]$kallus[
        test_indices,
      ] <- scenario_fold$kallus

      scenario_predictions[[scenario]]$bridge_constant[
        test_indices,
      ] <- scenario_fold$bridge_constant

      rm(scenario_fold, scenario_data)
      invisible(gc(full = TRUE))
    }

    rm(shared_fold, train_indices, test_indices)
    invisible(gc(full = TRUE))
  }

  if (anyNA(shared_predictions)) {
    stop("Missing shared cross-fitted predictions remain.")
  }

  for (scenario in scenario_names) {
    if (
      anyNA(scenario_predictions[[scenario]]$kallus) ||
        anyNA(scenario_predictions[[scenario]]$bridge_constant)
    ) {
      stop(
        "Missing scenario predictions remain for ",
        scenario,
        "."
      )
    }
  }

  list(
    shared = shared_predictions,
    scenarios = scenario_predictions
  )
}

estimate_kallus_style <- function(data, prediction) {
  e_hat <- clip_probability(
    prediction$e,
    0.001,
    0.999
  )

  r_hat <- clip_probability(
    prediction$r,
    1e-6,
    1
  )

  residual1 <- ifelse(
    data$R == 1,
    data$Y_obs - prediction$m1,
    0
  )

  residual0 <- ifelse(
    data$R == 1,
    data$Y_obs - prediction$m0,
    0
  )

  weight_y1 <- data$T * data$R / (e_hat * r_hat)
  weight_y0 <- (1 - data$T) * data$R / ((1 - e_hat) * r_hat)
  weight_s1 <- data$T / e_hat
  weight_s0 <- (1 - data$T) / (1 - e_hat)

  norm_y1 <- mean(weight_y1)
  norm_y0 <- mean(weight_y0)
  norm_s1 <- mean(weight_s1)
  norm_s0 <- mean(weight_s0)

  score <- (
    prediction$mu1 -
      prediction$mu0 +
      weight_y1 * residual1 / norm_y1 -
      weight_y0 * residual0 / norm_y0 +
      weight_s1 * (prediction$m1 - prediction$mu1) / norm_s1 -
      weight_s0 * (prediction$m0 - prediction$mu0) / norm_s0
  )

  data.frame(
    estimate = mean(score),
    standard_error = sqrt(var(score) / nrow(data))
  )
}

estimate_oi_style <- function(data, prediction) {
  e_hat <- clip_probability(
    prediction$e,
    0.001,
    0.999
  )

  weight1 <- data$T / e_hat
  weight0 <- (1 - data$T) / (1 - e_hat)

  norm1 <- mean(weight1)
  norm0 <- mean(weight0)

  score <- (
    prediction$q1 -
      prediction$q0 +
      weight1 * (prediction$h - prediction$q1) / norm1 -
      weight0 * (prediction$h - prediction$q0) / norm0
  )

  data.frame(
    estimate = mean(score),
    standard_error = sqrt(var(score) / nrow(data))
  )
}

calculate_outcome_rmse <- function(
    data,
    m1_hat,
    m0_hat
) {
  sqrt(
    mean(
      c(
        (m1_hat - data$m1_true)^2,
        (m0_hat - data$m0_true)^2
      )
    )
  )
}

calculate_correction_rmse <- function(
    data,
    d1_hat,
    d0_hat
) {
  sqrt(
    mean(
      c(
        (d1_hat - data$d1_potential_true)^2,
        (d0_hat - data$d0_potential_true)^2
      )
    )
  )
}

safe_mean <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}

format_method_result <- function(
    estimation_result,
    method,
    outcome_rmse,
    correction_rmse = NA_real_,
    mean_delta1_hat = NA_real_,
    mean_delta0_hat = NA_real_
) {
  estimation_result$method <- method
  estimation_result$outcome_rmse <- outcome_rmse
  estimation_result$correction_rmse <- correction_rmse
  estimation_result$mean_delta1_hat <- mean_delta1_hat
  estimation_result$mean_delta0_hat <- mean_delta0_hat
  estimation_result
}

assemble_scenario_results <- function(
    E_data,
    shared_prediction,
    scenario_prediction,
    label_budget,
    replication_id,
    history_best_trees
) {
  oracle_prediction <- data.frame(
    e = E_data$e_true,
    r = E_data$r_true,
    m1 = E_data$m1_true,
    m0 = E_data$m0_true,
    mu1 = E_data$mu1_true,
    mu0 = E_data$mu0_true
  )

  oi_prediction <- data.frame(
    e = shared_prediction$e,
    h = shared_prediction$h,
    q1 = shared_prediction$q1,
    q0 = shared_prediction$q0
  )

  zero_correction <- rep(0, nrow(E_data))

  history_h_rmse <- sqrt(
    mean(
      (shared_prediction$h - E_data$h_observed_true)^2
    )
  )

  results <- dplyr::bind_rows(
    format_method_result(
      estimate_kallus_style(
        E_data,
        oracle_prediction
      ),
      method = "Kallus Oracle",
      outcome_rmse = 0
    ),

    format_method_result(
      estimate_kallus_style(
        E_data,
        scenario_prediction$kallus
      ),
      method = "Kallus-GBM",
      outcome_rmse = calculate_outcome_rmse(
        E_data,
        scenario_prediction$kallus$m1,
        scenario_prediction$kallus$m0
      )
    ),

    format_method_result(
      estimate_oi_style(
        E_data,
        oi_prediction
      ),
      method = "OI-GBM",
      outcome_rmse = calculate_outcome_rmse(
        E_data,
        oi_prediction$h,
        oi_prediction$h
      ),
      correction_rmse = calculate_correction_rmse(
        E_data,
        zero_correction,
        zero_correction
      ),
      mean_delta1_hat = 0,
      mean_delta0_hat = 0
    ),

    format_method_result(
      estimate_kallus_style(
        E_data,
        scenario_prediction$bridge_constant
      ),
      method = "Bridge-Constant",
      outcome_rmse = calculate_outcome_rmse(
        E_data,
        scenario_prediction$bridge_constant$m1,
        scenario_prediction$bridge_constant$m0
      ),
      correction_rmse = calculate_correction_rmse(
        E_data,
        scenario_prediction$bridge_constant$d1,
        scenario_prediction$bridge_constant$d0
      ),
      mean_delta1_hat = mean(
        scenario_prediction$bridge_constant$d1
      ),
      mean_delta0_hat = mean(
        scenario_prediction$bridge_constant$d0
      )
    )
  )

  results$scenario <- E_data$scenario[[1L]]
  results$lambda <- E_data$lambda[[1L]]
  results$label_budget <- label_budget
  results$actual_labeled_N <- sum(E_data$R)
  results$replication <- replication_id
  results$true_ATE <- TRUE_TAU
  results$sample_true_ATE <- mean(E_data$Y1 - E_data$Y0)
  results$history_h_rmse <- history_h_rmse
  results$history_best_trees <- history_best_trees
  results$target_mean_d1_true <- mean(E_data$d1_potential_true)
  results$target_mean_d0_true <- mean(E_data$d0_potential_true)
  results$target_sd_d1_true <- sd(E_data$d1_potential_true)
  results$target_sd_d0_true <- sd(E_data$d0_potential_true)

  results
}

atomic_save_rds <- function(object, final_path) {
  temporary_path <- tempfile(
    pattern = "temporary_checkpoint_",
    tmpdir = dirname(final_path),
    fileext = ".rds"
  )

  saveRDS(
    object,
    temporary_path,
    compress = FALSE
  )

  rename_success <- file.rename(
    temporary_path,
    final_path
  )

  if (!rename_success) {
    if (file.exists(temporary_path)) {
      file.remove(temporary_path)
    }

    stop("Could not create checkpoint: ", final_path)
  }

  invisible(final_path)
}

checkpoint_path <- function(
    scenario,
    label_budget,
    replication_id
) {
  file.path(
    CHECKPOINT_DIRECTORY,
    sprintf(
      "scenario_%s_budget_%04d_rep_%03d.rds",
      scenario,
      label_budget,
      replication_id
    )
  )
}

error_path <- function(label_budget, replication_id) {
  file.path(
    ERROR_DIRECTORY,
    sprintf(
      "budget_%04d_rep_%03d.txt",
      label_budget,
      replication_id
    )
  )
}

run_replication_bundle <- function(
    label_budget,
    replication_id,
    scenarios = SCENARIOS
) {
  pending_scenarios <- scenarios[
    !vapply(
      scenarios,
      function(current_scenario) {
        file.exists(
          checkpoint_path(
            current_scenario,
            label_budget,
            replication_id
          )
        )
      },
      logical(1)
    )
  ]

  if (length(pending_scenarios) == 0L) {
    return(
      data.frame(
        label_budget = label_budget,
        replication = replication_id,
        status = "already_complete",
        seconds = 0
      )
    )
  }

  start_time <- Sys.time()

  current_seed <- (
    BASE_SEED +
      1000L * label_budget +
      replication_id
  )

  message(
    "Starting budget ",
    label_budget,
    " | replication ",
    replication_id,
    " | pending lambdas: ",
    paste(pending_scenarios, collapse = ", ")
  )

  H_data <- simulate_heterogeneous_history(
    n_H = N_H,
    seed = current_seed + 1L
  )

  history_model <- fit_gbm_learner(
    formula = NONLINEAR_FORMULAS$history_outcome,
    data = H_data,
    tree_selection = "fixed",
    max_trees = GBM_MAX_TREES
  )

  history_best_trees <- history_model$best_trees

  rm(H_data)
  invisible(gc(full = TRUE))

  E_data_list <- lapply(
    pending_scenarios,
    function(current_scenario) {
      simulate_heterogeneous_target(
        n_E = N_E,
        label_budget = label_budget,
        scenario = current_scenario,
        seed = current_seed
      )
    }
  )

  names(E_data_list) <- pending_scenarios

  predictions <- crossfit_lambda_bundle(
    E_data_list = E_data_list,
    history_model = history_model,
    K = K_FOLDS,
    seed = current_seed + 2L
  )

  rm(history_model)
  invisible(gc(full = TRUE))

  for (scenario in pending_scenarios) {
    current_result <- assemble_scenario_results(
      E_data = E_data_list[[scenario]],
      shared_prediction = predictions$shared,
      scenario_prediction = predictions$scenarios[[scenario]],
      label_budget = label_budget,
      replication_id = replication_id,
      history_best_trees = history_best_trees
    )

    atomic_save_rds(
      current_result,
      checkpoint_path(
        scenario,
        label_budget,
        replication_id
      )
    )

    rm(current_result)
  }

  rm(E_data_list, predictions)
  invisible(gc(full = TRUE))

  elapsed_seconds <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )

  message(
    "Finished budget ",
    label_budget,
    " | replication ",
    replication_id,
    " | minutes: ",
    round(elapsed_seconds / 60, 2)
  )

  data.frame(
    label_budget = label_budget,
    replication = replication_id,
    status = "completed",
    seconds = elapsed_seconds
  )
}

run_replication_bundle_safe <- function(
    label_budget,
    replication_id,
    scenarios = SCENARIOS
) {
  tryCatch(
    run_replication_bundle(
      label_budget = label_budget,
      replication_id = replication_id,
      scenarios = scenarios
    ),
    error = function(error_condition) {
      error_message <- paste0(
        format(Sys.time()),
        "\n",
        "Budget: ",
        label_budget,
        "\nReplication: ",
        replication_id,
        "\nError: ",
        conditionMessage(error_condition),
        "\n"
      )

      writeLines(
        error_message,
        con = error_path(label_budget, replication_id)
      )

      message(error_message)

      data.frame(
        label_budget = label_budget,
        replication = replication_id,
        status = "error",
        seconds = NA_real_
      )
    }
  )
}

all_scenario_checkpoints_exist <- function(
    label_budget,
    replication_id,
    scenarios = SCENARIOS
) {
  all(
    vapply(
      scenarios,
      function(current_scenario) {
        file.exists(
          checkpoint_path(
            current_scenario,
            label_budget,
            replication_id
          )
        )
      },
      logical(1)
    )
  )
}

build_job_grid <- function(mode) {
  if (mode == "pilot") {
    return(
      data.frame(
        label_budget = 200L,
        replication = 1L
      )
    )
  }

  replications <- if (mode == "fast") {
    5L
  } else {
    50L
  }

  expand.grid(
    label_budget = LABEL_BUDGETS,
    replication = seq_len(replications),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(label_budget, replication)
}

run_parallel_jobs <- function(job_grid, workers) {
  pending <- !vapply(
    seq_len(nrow(job_grid)),
    function(index) {
      all_scenario_checkpoints_exist(
        label_budget = job_grid$label_budget[[index]],
        replication_id = job_grid$replication[[index]]
      )
    },
    logical(1)
  )

  pending_grid <- job_grid[pending, , drop = FALSE]

  if (nrow(pending_grid) == 0L) {
    message("All requested jobs already have checkpoints.")
    return(invisible(data.frame()))
  }

  workers <- max(
    1L,
    min(as.integer(workers), nrow(pending_grid))
  )

  message(
    "Pending budget-replication bundles: ",
    nrow(pending_grid),
    " | workers: ",
    workers
  )

  future::plan(
    future::multisession,
    workers = workers
  )

  on.exit(
    future::plan(future::sequential),
    add = TRUE
  )

  status_list <- future.apply::future_lapply(
    seq_len(nrow(pending_grid)),
    function(index) {
      run_replication_bundle_safe(
        label_budget = pending_grid$label_budget[[index]],
        replication_id = pending_grid$replication[[index]]
      )
    },
    future.seed = TRUE,
    future.scheduling = 2,
    future.packages = c("gbm", "splitTools", "dplyr")
  )

  status <- dplyr::bind_rows(status_list)

  saveRDS(
    status,
    file.path(RESULT_DIRECTORY, "latest_run_status.rds")
  )

  utils::write.csv(
    status,
    file.path(RESULT_DIRECTORY, "latest_run_status.csv"),
    row.names = FALSE
  )

  status
}

collect_and_summarize <- function() {
  checkpoint_files <- list.files(
    CHECKPOINT_DIRECTORY,
    pattern = paste0(
      "^scenario_.*_budget_[0-9]+_rep_[0-9]+\\.rds$"
    ),
    full.names = TRUE
  )

  if (length(checkpoint_files) == 0L) {
    message("No checkpoint files found yet.")
    return(invisible(NULL))
  }

  heterogeneous_results <- checkpoint_files |>
    lapply(readRDS) |>
    dplyr::bind_rows() |>
    dplyr::arrange(
      lambda,
      label_budget,
      replication,
      method
    )

  critical_value <- qnorm(0.975)

  heterogeneous_summary <- heterogeneous_results |>
    dplyr::group_by(
      scenario,
      lambda,
      label_budget,
      method
    ) |>
    dplyr::summarise(
      successful_replications = sum(is.finite(estimate)),
      mean_labeled_N = mean(actual_labeled_N, na.rm = TRUE),
      mean_estimate = mean(estimate, na.rm = TRUE),
      bias = mean(estimate - TRUE_TAU, na.rm = TRUE),
      empirical_sd = sd(estimate, na.rm = TRUE),
      rmse = sqrt(
        mean(
          (estimate - TRUE_TAU)^2,
          na.rm = TRUE
        )
      ),
      mean_standard_error = mean(standard_error, na.rm = TRUE),
      coverage = mean(
        estimate - critical_value * standard_error <= TRUE_TAU &
          estimate + critical_value * standard_error >= TRUE_TAU,
        na.rm = TRUE
      ),
      mean_ci_width = mean(
        2 * critical_value * standard_error,
        na.rm = TRUE
      ),
      mean_outcome_rmse = mean(outcome_rmse, na.rm = TRUE),
      mean_correction_rmse = safe_mean(correction_rmse),
      mean_history_h_rmse = mean(history_h_rmse, na.rm = TRUE),
      mean_delta1_hat = safe_mean(mean_delta1_hat),
      mean_delta0_hat = safe_mean(mean_delta0_hat),
      mean_target_d1_true = mean(target_mean_d1_true, na.rm = TRUE),
      mean_target_d0_true = mean(target_mean_d0_true, na.rm = TRUE),
      mean_target_sd_d1 = mean(target_sd_d1_true, na.rm = TRUE),
      mean_target_sd_d0 = mean(target_sd_d0_true, na.rm = TRUE),
      mean_sample_true_ATE = mean(sample_true_ATE, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      se_to_empirical_sd = mean_standard_error / empirical_sd
    ) |>
    dplyr::arrange(
      lambda,
      label_budget,
      method
    )

  completion_check <- heterogeneous_results |>
    dplyr::group_by(
      scenario,
      lambda,
      label_budget,
      method
    ) |>
    dplyr::summarise(
      completed_replications = dplyr::n_distinct(replication),
      finite_estimates = sum(is.finite(estimate)),
      .groups = "drop"
    ) |>
    dplyr::arrange(lambda, label_budget, method)

  saveRDS(
    heterogeneous_results,
    file.path(RESULT_DIRECTORY, "lambda_transport_results.rds")
  )

  saveRDS(
    heterogeneous_summary,
    file.path(RESULT_DIRECTORY, "lambda_transport_summary.rds")
  )

  utils::write.csv(
    heterogeneous_results,
    file.path(RESULT_DIRECTORY, "lambda_transport_results.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    heterogeneous_summary,
    file.path(RESULT_DIRECTORY, "lambda_transport_summary.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    completion_check,
    file.path(RESULT_DIRECTORY, "completion_check.csv"),
    row.names = FALSE
  )

  message(
    "Saved summaries in: ",
    RESULT_DIRECTORY
  )

  print(
    heterogeneous_summary |>
      dplyr::select(
        lambda,
        label_budget,
        method,
        successful_replications,
        mean_estimate,
        bias,
        empirical_sd,
        rmse,
        mean_standard_error,
        coverage,
        mean_outcome_rmse,
        mean_correction_rmse
      ),
    n = Inf,
    width = Inf
  )

  invisible(
    list(
      results = heterogeneous_results,
      summary = heterogeneous_summary,
      completion = completion_check
    )
  )
}

arguments <- commandArgs(trailingOnly = TRUE)

RUN_MODE <- if (length(arguments) >= 1L) {
  tolower(arguments[[1L]])
} else {
  "pilot"
}

valid_modes <- c(
  "pilot",
  "fast",
  "full",
  "summary"
)

if (!RUN_MODE %in% valid_modes) {
  stop(
    "Mode must be one of: ",
    paste(valid_modes, collapse = ", ")
  )
}

REQUESTED_WORKERS <- if (length(arguments) >= 2L) {
  as.integer(arguments[[2L]])
} else {
  24L
}

if (
  length(REQUESTED_WORKERS) != 1L ||
    is.na(REQUESTED_WORKERS) ||
    REQUESTED_WORKERS < 1L
) {
  stop("The worker count must be a positive integer.")
}

message("Run mode: ", RUN_MODE)
message("Script directory: ", BASE_DIRECTORY)
message("Result directory: ", RESULT_DIRECTORY)

print(SHIFT_CALIBRATION_DIAGNOSTIC)

if (RUN_MODE != "summary") {
  requested_grid <- build_job_grid(RUN_MODE)

  run_status <- run_parallel_jobs(
    job_grid = requested_grid,
    workers = REQUESTED_WORKERS
  )

  if (nrow(run_status) > 0L) {
    print(table(run_status$status, useNA = "ifany"))
  }
}

collect_and_summarize()


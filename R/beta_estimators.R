# Helper function to fit common beta (homogeneous coefficients) model
# via pooled OLS across all units (without lagged outcomes for now)
#' Fit Pooled OLS Model for Homogeneous Covariate Effects
#'
#' @param data A data frame.
#' @param outcome_var Name of the outcome variable.
#' @param covariates A character vector of covariate names.
#' @param treat_time_var Name of the treatment adoption time variable.
#' @param time_var Name of the time variable.
#' @param unit_var Name of the unit identifier variable.
#'
#' @return A vector of coefficients (beta_hat).
fit_common_beta_ols <- function(data,
                                outcome_var,
                                covariate_vars,
                                treat_time_var,
                                time_var,
                                unit_var,
                                degree = 0,
                                pretreatment_window = c("full", "minimal")) {
  pretreatment_window <- match.arg(pretreatment_window)

  # Add time_to_treat
  data <- dplyr::mutate(data, time_to_treat = .data[[time_var]] - .data[[treat_time_var]])

  data <- if (pretreatment_window == "minimal") {
    dplyr::filter(data, time_to_treat <= 0 & time_to_treat >= -degree)
  } else {
    dplyr::filter(data, time_to_treat < 0)
  }

  # Filter to pre-treatment window
  data <- dplyr::filter(data, time_to_treat <= 0 & time_to_treat >= -degree)

  # Build formula
  cov_formula <- as.formula(
    paste(outcome_var, "~", paste(covariate_vars, collapse = " + "))
  )

  # Fit pooled OLS model
  model <- lm(formula = cov_formula, data = data)

  return(stats::coef(model))
}


#' Fit Pooled IV Model with Lagged Outcome (Homogeneous Beta)
#'
#' This function estimates a pooled IV model for covariate adjustment.
#' It instruments the lagged outcome with earlier lags and returns only
#' the beta coefficients for the included covariates (including the intercept).
#'
#' @param data A data.frame with panel structure.
#' @param covariate_vars Character vector of exogenous covariate names.
#' @param outcome_var Name of the outcome variable.
#' @param time_var Name of the time variable.
#' @param unit_var Name of the unit identifier.
#' @param treat_time_var Name of the treatment adoption time variable.
#' @param min_iv_lag Minimum lag to use as instrument (e.g., 2).
#' @param max_iv_lag Maximum lag to use as instrument (e.g., 2 or more).
#' @param degree Polynomial degree for determining pre-treatment window.
#'
#' @return A named vector of beta coefficients for the covariates and intercept.
#' @export
fit_common_beta_iv <- function(data,
                               covariate_vars,
                               outcome_var,
                               time_var,
                               unit_var,
                               treat_time_var,
                               min_iv_lag = 2,
                               max_iv_lag = 2,
                               degree = 0) {
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("AER", quietly = TRUE)

  # Arrange by unit and time
  data <- dplyr::arrange(data, .data[[unit_var]], .data[[time_var]])

  # Create lags of the outcome variable and name them accordingly
  for (l in 1:max(c(1, max_iv_lag))) {
    lag_name <- paste0("lag", l, "_", outcome_var)
    data[[lag_name]] <- dplyr::lag(data[[outcome_var]], n = l)
  }

  # Define lagged outcome (first lag)
  data$lagged_outcome <- data[[paste0("lag1_", outcome_var)]]

  # Determine pretreatment window using polynomial degree
  data <- dplyr::mutate(data, time_to_treat = .data[[time_var]] - .data[[treat_time_var]])
  data <- dplyr::filter(data, time_to_treat <= 0 & time_to_treat >= -degree)

  # Drop incomplete cases
  instruments <- paste0("lag", min_iv_lag:max_iv_lag, "_", outcome_var)
  all_needed <- c("lagged_outcome", instruments, covariate_vars, outcome_var)
  data <- data[stats::complete.cases(data[, all_needed]), ]

  # Build IV formula with explicit intercept
  rhs <- paste(c("1", "lagged_outcome", covariate_vars), collapse = " + ")
  rhs_instr <- paste(c(instruments, covariate_vars), collapse = " + ")
  fml <- stats::as.formula(paste0(outcome_var, " ~ ", rhs, " | ", rhs_instr))

  # Estimate IV model
  iv_model <- AER::ivreg(fml, data = data)

  # Return beta_hat including intercept and covariates (not the endogenous lag)
  beta_full <- stats::coef(iv_model)
  keep <- c("(Intercept)", covariate_vars)
  return(beta_full[keep])
}



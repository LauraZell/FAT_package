#' @title Fit Unit-Specific Polynomial Trend with Unitwise Coefficients
#' @description
#' This function estimates a polynomial trend for a single unit, adjusting for
#' time-varying covariates with **unit-specific** coefficients (βᵢ).
#' For each unit, it fits the model:
#' \deqn{ y_{it} = \sum_{k=0}^{q} c_k (t - t_0)^k + x_{it}' βᵢ + ε_{it} }
#' where the coefficients \eqn{βᵢ} are estimated from pre-treatment data only.
#'
#' This approach allows for heterogeneous covariate effects across units.
#'
#' @param data A `data.frame` containing the panel data.
#' @param unit A single unit identifier (e.g., a state or municipality).
#' @param unit_var Character: name of the unit variable.
#' @param time_var Character: name of the time variable.
#' @param outcome_var Character: name of the outcome variable.
#' @param treat_time_var Character: name of the treatment time variable.
#' @param degree Integer: degree of the polynomial trend to fit.
#' @param covariate_vars Character vector of time-varying covariate names.
#'
#' @return A `tibble` with columns for unit, time, treatment time, actual outcome,
#' and fitted trend values (`preds`).
#'
#' @export

#'
fit_unitwise_trend <- function(data,
                                    unit,
                                    unit_var,
                                    time_var,
                                    outcome_var,
                                    treat_time_var,
                                    degree,
                                    covariate_vars = NULL) {

  unit_data <- data[data[[unit_var]] == unit, ]
  treat_year <- unique(unit_data[[treat_time_var]])
  if (length(treat_year) != 1 || is.na(treat_year)) {
    stop(paste("Invalid treatment year for unit:", unit))
  }

  # Subset to pre-treatment observations
  pre_data <- unit_data[unit_data[[time_var]] <= treat_year, ]

  # Build formula
  rhs_terms <- paste0("poly(", time_var, ", ", degree, ", raw = TRUE)")
  if (!is.null(covariate_vars)) {
    rhs_terms <- paste(rhs_terms, paste(covariate_vars, collapse = " + "), sep = " + ")
  }
  fml <- as.formula(paste(outcome_var, "~", rhs_terms))

  # Estimate the model on pre-treatment data
  mod <- stats::lm(fml, data = pre_data)

  # Predict on full data for that unit
  unit_data$preds <- stats::predict(mod, newdata = unit_data)

  # Return tidy output
  tibble::tibble(
    !!unit_var := unit_data[[unit_var]],
    !!time_var := unit_data[[time_var]],
    !!treat_time_var := unit_data[[treat_time_var]],
    !!outcome_var := unit_data[[outcome_var]],
    preds = unit_data$preds
  )
}




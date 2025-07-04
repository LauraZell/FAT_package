#' Fit Polynomial Regression for One Unit (with optional covariates)
#'
#' This function fits a unit-specific polynomial trend to pre-treatment data.
#' If covariates are provided, it subtracts the expected contribution of those
#' covariates using a common beta (pooled across units) before estimating the trend.
#'
#' @param data The full data frame.
#' @param unit A single unit (e.g., state or municipality) to subset on.
#' @param unit_var Name of the unit variable (string).
#' @param time_var Name of the time variable (string).
#' @param outcome_var Name of the outcome variable (string).
#' @param treat_time_var Name of the treatment year variable (string).
#' @param degree Degree of the polynomial used for fitting (integer).
#' @param covariate_vars Optional character vector of covariate names.
#' @param beta_hat Optional vector of estimated coefficients for covariates.
#'
#' @return A data frame containing all time periods for that unit and a `preds` column.
#' @export
fit_common_trend <- function(data,
                           unit,
                           unit_var,
                           time_var,
                           outcome_var,
                           treat_time_var,
                           degree,
                           covariate_vars = NULL,
                           beta_hat = NULL) {
  # Subset for the specific unit
  unit_data <- data[data[[unit_var]] == unit, ]
  unit_data[[time_var]] <- as.numeric(as.character(unit_data[[time_var]]))
  # Determine the minimum observed year for centering the time trend
  min_year <- min(unit_data[[time_var]], na.rm = TRUE)
  # Extract this unit’s treatment year (should be one fixed value)
  treat_year <- unique(unit_data[[treat_time_var]])
  if (length(treat_year) != 1 || is.na(treat_year)) {
    stop(paste("Invalid treatment year for unit:", unit))
  }

  # Select pre-treatment data (<= treatment year)
  pre_data <- unit_data[unit_data[[time_var]] <= treat_year, ]

  # Adjust outcome if covariates and beta_hat are supplied
  if (!is.null(covariate_vars) && !is.null(beta_hat)) {
    # Create design matrix X for the covariates in pre-treatment data
    X <- model.matrix(
      as.formula(paste("~", paste(covariate_vars, collapse = " + "))),
      data = pre_data
    )

    # Diagnostic checks
    if (nrow(pre_data) != nrow(X)) {
      stop("Mismatch between rows in pre_data and covariate matrix X")
    }
    if (length(beta_hat) != ncol(X)) {
      stop("Length of beta_hat does not match number of covariates (including intercept)")
    }

    # Subtract predicted covariate effects (pooled across units) from the outcome
    # → This gives us "covariate-adjusted" residual variation
    pre_data$outcome_adj <- pre_data[[outcome_var]] - as.vector(X %*% beta_hat)

    # Repeat adjustment for full data (for prediction)
    X_pred <- model.matrix(
      as.formula(paste("~", paste(covariate_vars, collapse = " + "))),
      data = unit_data
    )
    unit_data$outcome_adj <- unit_data[[outcome_var]] - as.vector(X_pred %*% beta_hat)

  } else {
    # If no covariates are supplied, the adjusted outcome is just the raw outcome
    pre_data$outcome_adj <- pre_data[[outcome_var]]
    unit_data$outcome_adj <- unit_data[[outcome_var]]
  }

  # Fit Polynomial Trend to Pre-Treatment (Covariate-Adjusted) Outcome
  if (degree == 0) {
    # A constant pre-treatment average (intercept-only model)
    unit_data$preds <- mean(pre_data$outcome_adj, na.rm = TRUE)
  } else {
    # A unit-specific polynomial regression centered around min_year
    f <- stats::as.formula(
      paste0("outcome_adj ~ poly(I(", time_var, " - ", min_year, "), ", degree, ", raw = TRUE)")
    )
    mod <- stats::lm(f, data = pre_data)
    unit_data$preds <- stats::predict(mod, newdata = unit_data)
  }

  # Return output in tidy format
  tibble::tibble(
    !!unit_var := unit_data[[unit_var]],
    !!time_var := unit_data[[time_var]],
    !!treat_time_var := unit_data[[treat_time_var]],
    !!outcome_var := unit_data[[outcome_var]],
    preds = unit_data$preds
  )
}

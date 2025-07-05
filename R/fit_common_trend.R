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
                             deg,
                             unit_var,
                             time_var,
                             outcome_var,
                             treat_time_var,
                             covariate_vars = NULL,
                             beta_hat = NULL) {

  # Filter for the unit of interest
  unit_data <- data[data[[unit_var]] == unit, ]

  # Fix: Convert time variable to numeric (handles DFAT-related join issues)
  if (is.factor(unit_data[[time_var]])) {
    unit_data[[time_var]] <- as.numeric(as.character(unit_data[[time_var]]))
  }

  # Fix: Convert treatment time variable to numeric if needed
  if (is.factor(unit_data[[treat_time_var]])) {
    unit_data[[treat_time_var]] <- as.numeric(as.character(unit_data[[treat_time_var]]))
  }

  # Remove rows with missing values in the time or treatment time variable
  if (any(is.na(unit_data[[time_var]])) || any(is.na(unit_data[[treat_time_var]]))) {
    stop(paste("Invalid treatment year or time for unit:", unit))
  }

  # Center time variable relative to the treatment
  unit_data$timeToTreat <- unit_data[[time_var]] - unit_data[[treat_time_var]]

  # Construct polynomial time trend terms
  for (d in 1:deg) {
    unit_data[[paste0("ttreat", d)]] <- unit_data$timeToTreat^d
  }

  # Build regression formula
  rhs_terms <- paste0("ttreat", 1:deg)
  if (!is.null(covariate_vars)) {
    rhs_terms <- c(rhs_terms, covariate_vars)
  }
  rhs <- paste(rhs_terms, collapse = " + ")
  formula <- as.formula(paste(outcome_var, "~", rhs))

  # Fit model to pre-treatment data
  pre_data <- subset(unit_data, timeToTreat < 0)
  if (nrow(pre_data) < deg + 1) {
    # Not enough pre-treatment observations
    return(data.frame())
  }

  model <- lm(formula, data = pre_data)

  # Predict for full sample
  unit_data$preds <- predict(model, newdata = unit_data)

  # Return original outcome, predictions, and treatment time
  return(unit_data[, c(unit_var, time_var, outcome_var, "preds", treat_time_var)])
}

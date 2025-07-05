#' Fit Polynomial Regression for One Unit (with optional covariates and pooled beta subtraction)
#'
#' This function fits a unit-specific polynomial trend to pre-treatment data.
#' If covariates and a common beta_hat are provided (from pooled OLS or IV),
#' it subtracts the covariate effect from the outcome before estimating the trend.
#'
#' @param data The full data frame.
#' @param unit A single unit (e.g., state or municipality) to subset on.
#' @param unit_var Name of the unit variable (string).
#' @param time_var Name of the time variable (string).
#' @param outcome_var Name of the outcome variable (string).
#' @param treat_time_var Name of the treatment year variable (string).
#' @param degree Degree of the polynomial used for fitting (integer).
#' @param covariate_vars Optional character vector of covariate names.
#' @param beta_hat Optional named vector of estimated coefficients for covariates.
#'
#' @return A data frame with original outcome, prediction, and treatment time.
#' @export
fit_common_trend <- function(data,
                             unit,
                             degree,
                             unit_var,
                             time_var,
                             outcome_var,
                             treat_time_var,
                             covariate_vars = NULL,
                             beta_hat = NULL) {

  # Subset data for this unit
  unit_data <- data[data[[unit_var]] == unit, ]

  # Ensure time and treatment year are numeric
  unit_data[[time_var]] <- as.numeric(as.character(unit_data[[time_var]]))
  unit_data[[treat_time_var]] <- as.numeric(as.character(unit_data[[treat_time_var]]))

  if (any(is.na(unit_data[[time_var]])) || any(is.na(unit_data[[treat_time_var]]))) {
    stop(paste("Invalid treatment year or time for unit:", unit))
  }

  # Time to treatment
  unit_data$timeToTreat <- unit_data[[time_var]] - unit_data[[treat_time_var]]

  # Construct polynomial terms
  for (d in 1:degree) {
    unit_data[[paste0("ttreat", d)]] <- unit_data$timeToTreat^d
  }

  # Subtract pooled covariate effect if beta_hat is supplied
  if (!is.null(beta_hat) && !is.null(covariate_vars)) {
    # Check for missing variables
    missing_covars <- setdiff(covariate_vars, names(unit_data))
    missing_betas <- setdiff(covariate_vars, names(beta_hat))
    if (length(missing_covars) > 0 || length(missing_betas) > 0) {
      stop("Mismatch in covariate variables or missing beta_hat values.")
    }

    # Subtract pooled covariate contribution
    covariate_effect <- as.matrix(unit_data[, covariate_vars]) %*% beta_hat[covariate_vars]
    unit_data$adjusted_outcome <- unit_data[[outcome_var]] - covariate_effect
  } else {
    unit_data$adjusted_outcome <- unit_data[[outcome_var]]
  }

  # Fit model to pre-treatment data only
  pre_data <- subset(unit_data, timeToTreat < 0)
  if (nrow(pre_data) < degree + 1) return(data.frame())  # Not enough pre-treatment points

  # Build formula using polynomial time trend terms
  rhs_terms <- paste0("ttreat", 1:degree)
  rhs <- paste(rhs_terms, collapse = " + ")
  formula <- as.formula(paste("adjusted_outcome ~", rhs))

  # Fit model on pre-treatment data
  model <- lm(formula, data = pre_data)

  # ==== Predict only for post-treatment years ====
  post_data <- subset(unit_data, timeToTreat >= 0)
  post_data$preds <- predict(model, newdata = post_data)

  # Re-add covariate contribution if previously subtracted
  if (!is.null(beta_hat) && !is.null(covariate_vars)) {
    covariate_matrix <- as.matrix(post_data[, covariate_vars])
    covariate_effect <- as.numeric(covariate_matrix %*% beta_hat)
    post_data$preds <- post_data$preds + covariate_effect
  }

  # ==== Merge predictions back into full unit_data ====
  unit_data$preds <- NA  # default NA
  unit_data[unit_data$timeToTreat >= 0, "preds"] <- post_data$preds

  # Return full dataset with only post-treatment preds filled
  return(unit_data[, c(unit_var, time_var, outcome_var, "preds", treat_time_var)])
}

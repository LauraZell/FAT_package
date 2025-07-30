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
#' @param data The full dataset (data frame).
#' @param unit A single unit (e.g. state or municipality) to subset on.
#' @param degree Degree of the polynomial (e.g. 0, 1, 2).
#' @param unit_var Name of the unit identifier variable (string).
#' @param time_var Name of the time variable (string).
#' @param outcome_var Name of the outcome variable to predict (string).
#' @param treat_time_var Name of the treatment time variable (string).
#' @param covariate_vars Optional character vector of covariate names.
#' @param beta_hat Optional named vector of pooled covariate coefficients.
#'
#' @return A data frame with unit, time, outcome, predicted values (`preds`),
#'         and treatment time.
#' @export
fit_unitwise_trend <- function(data,
                               unit,
                               degree,
                               unit_var,
                               time_var,
                               outcome_var,
                               treat_time_var,
                               covariate_vars = NULL,
                               beta_hat = NULL,
                               forecast_from_treatment_year = FALSE,
                               pretreatment_window = c("full", "minimal")) {
  pretreatment_window <- match.arg(pretreatment_window)
  # Filter data to the unit of interest
  unit_data <- data[data[[unit_var]] == unit, ]

  # Ensure numeric time variables (avoids factor errors)
  unit_data[[time_var]] <- as.numeric(as.character(unit_data[[time_var]]))
  unit_data[[treat_time_var]] <- as.numeric(as.character(unit_data[[treat_time_var]]))

  # Create time-to-treatment variable
  unit_data$timeToTreat <- unit_data[[time_var]] - unit_data[[treat_time_var]]

  # Create polynomial terms for timeToTreat
  for (d in 1:degree) {
    unit_data[[paste0("ttreat", d)]] <- unit_data$timeToTreat^d
  }

  # Subtract covariate contribution if beta_hat is provided
  if (!is.null(beta_hat) && !is.null(covariate_vars)) {
    covariate_matrix <- as.matrix(unit_data[, covariate_vars, drop = FALSE])
    covariate_effect <- as.vector(covariate_matrix %*% beta_hat[covariate_vars])
    unit_data$adjusted_outcome <- unit_data[[outcome_var]] - covariate_effect
  } else {
    unit_data$adjusted_outcome <- unit_data[[outcome_var]]
  }

  # Keep only pre-treatment data for estimation
  if (pretreatment_window == "minimal") {
    pre_data <- subset(unit_data, timeToTreat < 0)
    pre_data <- dplyr::arrange(pre_data, desc(timeToTreat)) %>% head(degree + 1)
  } else {
    pre_data <- subset(unit_data, timeToTreat < 0)
  }

  if (nrow(pre_data) <= degree) {
    warning(paste("Unit", unit, "has too few pre-treatment observations for degree =", degree))
  }


  # Create regression formula (e.g. adjusted_outcome ~ ttreat1 + ttreat2)
  rhs_terms <- paste0("ttreat", 1:degree)
  rhs <- paste(rhs_terms, collapse = " + ")
  formula <- as.formula(paste("adjusted_outcome ~", rhs))

  # Fit model on pre-treatment data
  model <- lm(formula, data = pre_data)
  print(summary(model))

  # Create empty preds
  unit_data$preds <- NA_real_

  # Predict only for post-treatment periods (timeToTreat >= 0)
  post_data <- subset(unit_data, timeToTreat >= if (forecast_from_treatment_year) 0 else 1)
  preds <- predict(model, newdata = post_data)

  # Fill in predictions only for post-treatment years
  unit_data$preds[unit_data$timeToTreat >= 0] <- preds

  # Re-add covariate effects if subtracted
  if (!is.null(beta_hat) && !is.null(covariate_vars)) {
    unit_data$preds[unit_data$timeToTreat >= 0] <-
      unit_data$preds[unit_data$timeToTreat >= 0] +
      covariate_effect[unit_data$timeToTreat >= 0]
  }

  # Return relevant columns
  return(unit_data[, c(unit_var, time_var, outcome_var, "preds", treat_time_var)])
}

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
#' @param forecast_lag Number of periods to wait after treatment before forecasting begins.
#'        Default is 0 (forecast starts in treatment year).
#' @param pretreatment_window Character. Either "full" to use all available pre-treatment observations, or "minimal" to use exactly (degree + 1) most recent ones.
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
                             hh,
                             covariate_vars = NULL,
                             beta_hat = NULL,
                             forecast_lag = 0,
                             pretreatment_window = c("full", "minimal")) {
  pretreatment_window <- match.arg(pretreatment_window)

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
  if (pretreatment_window == "minimal") {
    pre_data <- subset(unit_data, timeToTreat < 0)
    pre_data <- dplyr::arrange(pre_data, desc(timeToTreat)) %>% head(degree + 1)
  } else {
    pre_data <- subset(unit_data, timeToTreat < 0)
  }

  # Include variables to check how many pretreatment years are included
  n_pre_fit <- nrow(pre_data)
  pre_years_used <- paste(sort(unique(pre_data$timeToTreat)), collapse = ",")

  # Need at least degree+1 DISTINCT pre-treatment points
  required_n <- degree + 1
  if (nrow(pre_data) < required_n || dplyr::n_distinct(pre_data$timeToTreat) < required_n) {
    warning(paste("Skipping unit:", unit, "due to insufficient pre-treatment observations."))

    return(
      dplyr::mutate(
        unit_data[, c(unit_var, time_var, outcome_var, treat_time_var)],
        timeToTreat = unit_data[[time_var]] - unit_data[[treat_time_var]],
        preds = NA_real_,
        hh    = 0L,           # will be overwritten later anyway
        deg   = degree,       # just for consistency
        skipped = TRUE,       # optional
        n_pre_fit = n_pre_fit,  # number of pre-treatment rows used
        pre_years_used = pre_years_used  # years used in pre-treatment fit
        )
    )
  }


  # Build formula using polynomial time trend terms and fit on pretreatment data
  if (degree == 0) {
    model <- lm(adjusted_outcome ~ 1, data = pre_data)
  } else {
    rhs_term <- paste0("ttreat", 1:degree)
    rhs <- paste(rhs_term, collapse = " + ")
    formula <- as.formula(paste("adjusted_outcome ~", rhs))
    model <- lm(formula, data = pre_data)
  }

  print(paste("Fitting model for unit:", unit, "with degree:", degree))
  print(summary(model))

  # ==== Predict only for post-treatment years ====
  # New: restrict to horizon hh
  post_data <- subset(unit_data,
                      timeToTreat >= forecast_lag &
                        timeToTreat <= forecast_lag + hh - 1)

  post_data$preds <- predict(model, newdata = post_data)

  print(post_data)
  print(post_data$preds)

  # Re-add covariate contribution if previously subtracted
  if (!is.null(beta_hat) && !is.null(covariate_vars)) {
    covariate_matrix <- as.matrix(post_data[, covariate_vars])
    covariate_effect <- as.numeric(covariate_matrix %*% beta_hat)
    post_data$preds <- post_data$preds + covariate_effect
  }

  # ==== Merge predictions back into full unit_data ====
  unit_data$preds <- NA_real_
  unit_data$preds[unit_data$timeToTreat %in% post_data$timeToTreat] <- post_data$preds


  # Add timeToTreat again (in case subset removed it), and compute hh indicator
  unit_data$timeToTreat <- unit_data[[time_var]] - unit_data[[treat_time_var]]
  unit_data$n_pre_fit <- n_pre_fit
  unit_data$pre_years_used <- pre_years_used

  unit_data$hh <- dplyr::if_else(
    unit_data$timeToTreat >= forecast_lag & unit_data$timeToTreat <= forecast_lag + hh - 1,
    unit_data$timeToTreat - forecast_lag + 1L,
    0L
  )

  # Return all years with preds (only filled for post-treatment forecast horizon)
  return(unit_data[, c(unit_var, time_var, outcome_var, "preds", treat_time_var, "timeToTreat", "hh", "n_pre_fit", "pre_years_used")])

}

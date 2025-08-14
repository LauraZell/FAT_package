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
#' @param forecast_lag Number of periods to wait after treatment before forecasting begins.
#'        Default is 0 (forecast starts in treatment year).
#' @param pretreatment_window Character. Either "full" to use all available pre-treatment observations, or "minimal" to use exactly (degree + 1) most recent ones.
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
                               hh,
                               covariate_vars = NULL,
                               beta_hat = NULL,          # ignored in unitwise; kept for signature parity
                               forecast_lag = 0,
                               pretreatment_window = c("full", "minimal")) {
  pretreatment_window <- match.arg(pretreatment_window)

  # Warn if a pooled beta_hat was (accidentally) passed
  if (!is.null(beta_hat)) {
    warning("fit_unitwise_trend(): beta_hat supplied but will be ignored in unitwise mode.")
  }

  # Subset unit
  unit_data <- data[data[[unit_var]] == unit, , drop = FALSE]

  # Ensure numeric
  unit_data[[time_var]]       <- as.numeric(as.character(unit_data[[time_var]]))
  unit_data[[treat_time_var]] <- as.numeric(as.character(unit_data[[treat_time_var]]))

  # time-to-treatment and polynomial terms
  unit_data$timeToTreat <- unit_data[[time_var]] - unit_data[[treat_time_var]]
  if (degree > 0) {
    for (d in 1:degree) {
      unit_data[[paste0("ttreat", d)]] <- unit_data$timeToTreat^d
    }
  }

  # Choose pre-treatment window
  pre_data <- subset(unit_data, timeToTreat < 0)
  if (pretreatment_window == "minimal") {
    # take the (degree+1) most recent distinct times (and keep all rows at those times)
    # Step 1: get the (degree+1) most recent distinct t's
    distinct_t <- sort(unique(pre_data$timeToTreat), decreasing = TRUE)
    keep_t     <- head(distinct_t, degree + 1)
    # Step 2: keep all rows with those t's (handles duplicates if any)
    pre_data   <- pre_data[pre_data$timeToTreat %in% keep_t, , drop = FALSE]
  }

  # Diagnostics
  n_pre_fit       <- nrow(pre_data)
  pre_years_used  <- paste(sort(unique(pre_data$timeToTreat)), collapse = ",")

  # Identification checks
  p <- if (is.null(covariate_vars)) 0L else length(covariate_vars)
  # number of parameters to estimate = intercept (1) + degree polynomial + p covariates
  required_n <- (degree + 1L) + p

  if (n_pre_fit < required_n || dplyr::n_distinct(pre_data$timeToTreat) < (degree + 1L)) {
    warning(paste("Skipping unit:", unit,
                  "due to insufficient pre-treatment observations (need at least",
                  required_n, "rows and", degree + 1L, "distinct time points)."))
    # return all years w/ NA preds so the caller can keep the pre-period for plotting
    return(
      dplyr::mutate(
        unit_data[, c(unit_var, time_var, outcome_var, treat_time_var), drop = FALSE],
        timeToTreat   = unit_data[[time_var]] - unit_data[[treat_time_var]],
        preds         = NA_real_,
        hh            = 0L,
        deg           = degree,
        skipped       = TRUE,
        n_pre_fit     = n_pre_fit,
        pre_years_used = pre_years_used
      )
    )
  }

  # Build RHS: intercept + polynomial terms + covariates
  poly_terms <- if (degree > 0) paste0("ttreat", 1:degree) else NULL
  rhs_terms  <- c(poly_terms, covariate_vars)
  # If there are no rhs terms (degree==0 & no covariates), keep intercept-only
  fml <- if (length(rhs_terms) > 0) {
    as.formula(paste(outcome_var, "~", paste(rhs_terms, collapse = " + ")))
  } else {
    as.formula(paste(outcome_var, "~ 1"))
  }

  # Fit on pre-treatment
  model <- lm(fml, data = pre_data)
  # Optional: print summary for debugging
  # print(summary(model))

  # Predict ONLY within the forecast horizon [forecast_lag, forecast_lag+hh-1]
  post_data <- subset(unit_data,
                      timeToTreat >= forecast_lag &
                        timeToTreat <= forecast_lag + hh - 1)

  unit_data$preds <- NA_real_
  if (nrow(post_data) > 0) {
    unit_data$preds[unit_data$timeToTreat %in% post_data$timeToTreat] <-
      predict(model, newdata = post_data)
  }

  # hh index for plotting (keep all rows; pre gets 0)
  unit_data$hh <- dplyr::if_else(
    unit_data$timeToTreat >= forecast_lag & unit_data$timeToTreat <= forecast_lag + hh - 1,
    unit_data$timeToTreat - forecast_lag + 1L,
    0L
  )

  # Attach diagnostics
  unit_data$n_pre_fit      <- n_pre_fit
  unit_data$pre_years_used <- pre_years_used

  # Return full time path (pre rows kept; preds only filled on horizon)
  unit_data[, c(unit_var, time_var, outcome_var, "preds",
                treat_time_var, "timeToTreat", "hh", "n_pre_fit", "pre_years_used")]
}

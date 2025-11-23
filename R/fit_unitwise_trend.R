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
#'
#' @param data The full dataset (data frame).
#' @param unit A single unit (e.g. state or municipality) to subset on.
#' @param degree Degree of the polynomial (e.g. 0, 1, 2).
#' @param unit_var Name of the unit identifier variable (string).
#' @param time_var Name of the time variable (string).
#' @param outcome_var Name of the outcome variable to predict (string).
#' @param treat_time_var Name of the treatment time variable (string).
#' @param covariate_vars Optional character vector of covariate names.
#' @param beta_hat Optional pooled coefficients (ignored in unitwise; kept for signature parity).
#' @param forecast_offset Integer; first forecast period is τ + forecast_offset.
#' @param pretreatment_window "full", "minimal", or "manual".
#' @param pretreatment_years If "manual", use the last K pre-treatment years (K ≥ q+1+p).
#' @param hh Horizon length: predict for t ∈ [τ+offset, τ+offset+hh−1].
#'
#' @return Data frame with: unit, time, outcome, preds, treat_time, timeToTreat, hh, n_pre_fit, pre_years_used.
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
                               forecast_offset = 0,
                               pretreatment_window = c("full", "minimal", "manual"),
                               pretreatment_years = NULL) {

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
  if (any(is.na(unit_data[[time_var]])) || any(is.na(unit_data[[treat_time_var]]))) {
    stop(sprintf("Invalid %s or %s for unit: %s", time_var, treat_time_var, unit))
  }

  # time-to-treatment and polynomial terms
  unit_data$timeToTreat <- unit_data[[time_var]] - unit_data[[treat_time_var]]
  if (degree > 0) {
    for (d in 1:degree) {
      unit_data[[paste0("ttreat", d)]] <- unit_data$timeToTreat^d
    }
  }
  # ---------------------------------------------------------------------------
  # Choose pre-treatment window
  pre_all <- unit_data[unit_data$timeToTreat < 0 & !is.na(unit_data[[outcome_var]]), , drop = FALSE]

  if (nrow(pre_all) == 0L) {
    # Return skeleton with required columns (no fit possible)
    out <- unit_data[, c(unit_var, time_var, outcome_var, treat_time_var), drop = FALSE]
    out$timeToTreat     <- unit_data$timeToTreat
    out$preds           <- NA_real_
    out$hh              <- 0L
    out$n_pre_fit       <- 0L
    out$pre_years_used  <- ""
    return(out)
  }

  # Identification checks: Number of unitwise parameters = intercept(1) + q polynomial terms + p covariates
  p_cov <- if (is.null(covariate_vars)) 0L else length(covariate_vars)
  required_n <- (degree + 1L) + p_cov

  if (pretreatment_window == "full") {
    pre_data <- pre_all[order(pre_all[[time_var]]), , drop = FALSE]

  } else if (pretreatment_window == "minimal") {
    # Take the last required_n DISTINCT pre-treatment periods
    distinct_t <- sort(unique(pre_all$timeToTreat), decreasing = TRUE)
    keep_t     <- head(distinct_t, required_n)
    pre_data   <- pre_all[pre_all$timeToTreat %in% keep_t, , drop = FALSE]
    # Keep chronological order
    pre_data   <- pre_data[order(pre_data[[time_var]]), , drop = FALSE]

  } else { # "manual"
    if (is.null(pretreatment_years) || !is.finite(pretreatment_years) || pretreatment_years < 1L) {
      stop("When pretreatment_window = 'manual', provide a valid 'pretreatment_years' (>= 1).")
    }
    if (pretreatment_years < required_n) {
      stop(sprintf("pretreatment_years = %d is too small for unitwise degree q = %d with %d covariate(s) (need at least q+1+p = %d).",
                   pretreatment_years, degree, p_cov, required_n))
    }
    # Take the last K DISTINCT pre periods (closest to treatment)
    distinct_t <- sort(unique(pre_all$timeToTreat), decreasing = TRUE)
    keep_t     <- head(distinct_t, as.integer(pretreatment_years))
    pre_data   <- pre_all[pre_all$timeToTreat %in% keep_t, , drop = FALSE]
    pre_data   <- pre_data[order(pre_data[[time_var]]), , drop = FALSE]
  }
  # ---------------------------------------------------------------------------

  # Diagnostics
  n_pre_fit       <- nrow(pre_data)
  pre_years_used  <- paste(sort(unique(pre_data$timeToTreat)), collapse = ",")


  if (n_pre_fit < required_n || dplyr::n_distinct(pre_data$timeToTreat) < (degree + 1L)) {
    warning(paste("Skipping unit:", unit,
                  "due to insufficient pre-treatment observations (need at least",
                  required_n, "rows and", degree + 1L, "distinct time points)."))
    out <- unit_data[, c(unit_var, time_var, outcome_var, treat_time_var), drop = FALSE]
    out$timeToTreat    = unit_data$timeToTreat
    out$preds          = NA_real_
    out$hh             = 0L
    out$deg            = degree
    out$n_pre_fit      = n_pre_fit
    out$pre_years_used = pre_years_used
    return(out)
  }
  if (n_pre_fit == required_n) {
    warning(sprintf("Saturated fit for unit %s at q = %d with %d covariate(s): usable pre rows = q+1+p = %d. This will result in a perfect fit with zero degrees of freedom, and the forecasts may be unstable.",
                    unit, degree, p_cov, required_n), call. = FALSE)
  }

  # Build RHS: intercept + polynomial terms + covariates
  poly_terms <- if (degree > 0) paste0("ttreat", 1:degree) else character(0)
  rhs_terms  <- c(poly_terms, if (!is.null(covariate_vars)) covariate_vars else character(0))

  # If there are no rhs terms (degree==0 & no covariates), keep intercept-only
  fml <- if (length(rhs_terms) > 0) {
    stats::as.formula(paste(outcome_var, "~", paste(rhs_terms, collapse = " + ")))
  } else {
    stats::as.formula(paste(outcome_var, "~ 1"))
  }

  model <- stats::lm(fml, data = pre_data)

  # Optional: print summary for debugging
  # print(summary(model))

  # Predict only over the requested horizon block [τ+forecast_offset, τ+forecast_offset+hh-1]
  post_data <- unit_data[
    unit_data$timeToTreat >= forecast_offset &
      unit_data$timeToTreat <= forecast_offset + hh - 1L,
    , drop = FALSE
  ]

  unit_data$preds <- NA_real_
  if (nrow(post_data) > 0L) {
    unit_data$preds[unit_data$timeToTreat %in% post_data$timeToTreat] <-
      stats::predict(model, newdata = post_data)
  }

  # hh index for plotting (keep all rows; pre gets 0)
  unit_data$hh <- dplyr::if_else(
    unit_data$timeToTreat >= forecast_offset & unit_data$timeToTreat <= forecast_offset + hh - 1,
    unit_data$timeToTreat - forecast_offset + 1L,
    0L
  )

  # Attach diagnostics
  unit_data$n_pre_fit      <- n_pre_fit
  unit_data$pre_years_used <- pre_years_used

  # Return full time path (pre rows kept; preds only filled on horizon)
  unit_data[, c(unit_var, time_var, outcome_var, "preds",
                treat_time_var, "timeToTreat", "hh", "n_pre_fit", "pre_years_used")]
}

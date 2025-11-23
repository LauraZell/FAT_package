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
#' @param degree Integer polynomial degree q (requires at least q+1 pre points).
#' @param covariate_vars Optional character vector of covariate names.
#' @param beta_hat Optional named vector of estimated coefficients for covariates.
#' @param forecast_offset Number of periods to wait after treatment before forecasting begins.
#'        Default is 0 (forecast starts in treatment year). Integer; first forecast period is tau + forecast_offset.
#' @param pretreatment_window Character. "minimal", "full", or "manual".
#' @param pretreatment_years If "manual", use the last K pre-treatment years (K >= q+1).
#' @param hh Horizon block length (predict for t in [tau+offset, tau+offset+hh-1]).
#' @return Data frame with columns: unit, time, outcome, preds, treat_time, timeToTreat, hh, n_pre_fit, pre_years_used.
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
                             forecast_offset = 0,
                             pretreatment_window = c("minimal", "full", "manual"),
                             pretreatment_years = NULL) {
  pretreatment_window <- match.arg(pretreatment_window)

  # Subset data for this unit
  unit_data <- data[data[[unit_var]] == unit, , drop = FALSE]

  # Ensure time and treatment year are numeric
  unit_data[[time_var]] <- as.numeric(as.character(unit_data[[time_var]]))
  unit_data[[treat_time_var]] <- as.numeric(as.character(unit_data[[treat_time_var]]))
  if (any(is.na(unit_data[[time_var]])) || any(is.na(unit_data[[treat_time_var]]))) {
    stop(paste("Invalid treatment year or time for unit:", unit))
  }

  # Time to treatment
  unit_data$timeToTreat <- unit_data[[time_var]] - unit_data[[treat_time_var]]

  # Construct polynomial terms
  if (degree > 0) {
    for (d in 1:degree) {
      unit_data[[paste0("ttreat", d)]] <- unit_data$timeToTreat^d
    }
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
    cov_mat <- as.matrix(unit_data[, covariate_vars, drop = FALSE])
    unit_data$adjusted_outcome <- as.numeric(unit_data[[outcome_var]] - (cov_mat %*% beta_hat[covariate_vars]))
  } else {
    unit_data$adjusted_outcome <- unit_data[[outcome_var]]
  }


  # ---------- Select pre-treatment rows according to the window rule ----------
  pre_all <- unit_data[unit_data$timeToTreat < 0 & !is.na(unit_data$adjusted_outcome), , drop = FALSE]

  if (nrow(pre_all) == 0L) {
    # Return a "skipped" frame with the required columns
    out <- unit_data[, c(unit_var, time_var, outcome_var, treat_time_var), drop = FALSE]
    out$timeToTreat   <- unit_data$timeToTreat
    out$preds         <- NA_real_
    out$hh            <- 0L
    out$n_pre_fit     <- 0L
    out$pre_years_used<- ""
    return(out)
  }

  if (pretreatment_window == "full") {
    pre_data <- pre_all[order(pre_all[[time_var]]), , drop = FALSE]

  } else if (pretreatment_window == "minimal") {
    need <- degree + 1L
    # take the last 'need' pre points closest to treatment
    ord  <- order(pre_all[[time_var]], decreasing = TRUE)
    take <- seq_len(min(need, length(ord)))
    idx  <- sort(ord[take], decreasing = FALSE)
    pre_data <- pre_all[idx, , drop = FALSE]

  } else { # "manual"
    if (is.null(pretreatment_years) || !is.finite(pretreatment_years) || pretreatment_years < 1L) {
      stop("When pretreatment_window = 'manual', provide a valid 'pretreatment_years' (>= 1).")
    }
    need <- degree + 1L
    if (pretreatment_years < need) {
      stop(sprintf("pretreatment_years = %d is too small for degree q = %d (need at least q+1 = %d).",
                   pretreatment_years, degree, need))
    }
    ord  <- order(pre_all[[time_var]], decreasing = TRUE)
    take <- seq_len(min(as.integer(pretreatment_years), length(ord)))
    idx  <- sort(ord[take], decreasing = FALSE)
    pre_data <- pre_all[idx, , drop = FALSE]
  }

  # ---------------------------------------------------------------------------

  # Identifiability & saturation checks *after* selection and NA drop
  required_n <- degree + 1L
  n_pre_fit  <- nrow(pre_data)
  distinct_t <- dplyr::n_distinct(pre_data$timeToTreat)

  pre_years_used <- paste(sort(unique(pre_data$timeToTreat)), collapse = ",")

  if (n_pre_fit < required_n || distinct_t < required_n) {
    # Not enough distinct pre points to fit degree q
    out <- unit_data[, c(unit_var, time_var, outcome_var, treat_time_var), drop = FALSE]
    out$timeToTreat    <- unit_data$timeToTreat
    out$preds          <- NA_real_
    out$hh             <- 0L
    out$n_pre_fit      <- n_pre_fit      # number of pre-treatment rows used
    out$pre_years_used <- pre_years_used
    return(out)
  }
  if (n_pre_fit == required_n) {
    warning(sprintf("Saturated fit for unit %s at q = %d (usable pre rows = q+1 = %d). Forecasts may be unstable.",
                    unit, degree, required_n), call. = FALSE)
  }


  # Fit polynomial on pre_data
  if (degree == 0) {
    model <- stats::lm(adjusted_outcome ~ 1, data = pre_data)
  } else {
    rhs  <- paste(paste0("ttreat", 1:degree), collapse = " + ")
    form <- stats::as.formula(paste("adjusted_outcome ~", rhs))
    model <- stats::lm(form, data = pre_data)
  }

  # ---------- Predict only for the requested post-treatment horizon block ----------
  post_data <- unit_data[
    unit_data$timeToTreat >= forecast_offset &
      unit_data$timeToTreat <= forecast_offset + hh - 1L,
    , drop = FALSE
  ]

  if (nrow(post_data) > 0L) {
    post_data$preds <- stats::predict(model, newdata = post_data)

    # Re-add pooled covariate effect if we subtracted it
    if (!is.null(beta_hat) && !is.null(covariate_vars)) {
      cov_mat_post <- as.matrix(post_data[, covariate_vars, drop = FALSE])
      post_data$preds <- as.numeric(post_data$preds + (cov_mat_post %*% beta_hat[covariate_vars]))
    }
  } else {
    post_data$preds <- numeric(0)
  }

  # Merge predictions back onto the unit grid
  unit_data$preds <- NA_real_
  if (nrow(post_data) > 0L) {
    # match by timeToTreat to avoid floating-point issues with time
    idx <- match(unit_data$timeToTreat, post_data$timeToTreat)
    unit_data$preds[!is.na(idx)] <- post_data$preds[idx[!is.na(idx)]]
  }

  # Book-keeping
  unit_data$n_pre_fit       <- n_pre_fit
  unit_data$pre_years_used  <- pre_years_used
  unit_data$hh <- dplyr::if_else(
    unit_data$timeToTreat >= forecast_offset & unit_data$timeToTreat <= forecast_offset + hh - 1L,
    unit_data$timeToTreat - forecast_offset + 1L,
    0L
  )

  # Final shape
  unit_data[, c(unit_var, time_var, outcome_var, "preds",
                treat_time_var, "timeToTreat", "hh", "n_pre_fit", "pre_years_used")]
}


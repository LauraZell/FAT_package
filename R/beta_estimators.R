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
                                pretreatment_window = c("full", "minimal", "manual"),
                                pretreatment_years = NULL) {

  pretreatment_window <- match.arg(pretreatment_window)

  # Build a pooled pre-treatment fram using the per-unit window rule
  build_pre_for_unit <- function(df_u) {
    pre_all <- df_u[df_u[[time_var]] - df_u[[treat_time_var]] < 0 &
                      !is.na(df_u[[outcome_var]]), , drop = FALSE]
    if (nrow(pre_all) == 0L) return(pre_all)

    if (pretreatment_window == "full") {
      pre_all[order(pre_all[[time_var]]), , drop = FALSE]

    } else if (pretreatment_window == "minimal") {
      need <- (degree + 1L) + max_iv_lag
      ord  <- order(pre_all[[time_var]], decreasing = TRUE)
      take <- seq_len(min(need, length(ord)))
      idx  <- sort(ord[take], decreasing = FALSE)
      pre_all[idx, , drop = FALSE]

    } else { # "manual"
      if (is.null(pretreatment_years) || !is.finite(pretreatment_years) || pretreatment_years < 1L) {
        stop("When pretreatment_window = 'manual', provide a valid 'pretreatment_years' (>= 1).")
      }
      need <- ((degree + 1L) + max_iv_lag)
      if (pretreatment_years < need) {
        stop(sprintf("pretreatment_years = %d is too small for degree q = %d (need at least q+1+max_iv_lag = %d).",
                     pretreatment_years, degree, need))
      }

      ord  <- order(pre_all[[time_var]], decreasing = TRUE)
      take <- seq_len(min(as.integer(pretreatment_years), length(ord)))
      idx  <- sort(ord[take], decreasing = FALSE)
      pre_all[idx, , drop = FALSE]
    }
  }

  # Split-apply-bind across units
  units <- unique(data[[unit_var]])
  pre_list <- lapply(units, function(u) build_pre_for_unit(data[data[[unit_var]] == u, , drop = FALSE]))
  pre_pool <- do.call(rbind, pre_list)

  if (is.null(pre_pool) || nrow(pre_pool) == 0L) {
    stop("Pooled OLS: no usable pre-treatment observations after windowing.")
  }

  # --- Identifiability check for pooled OLS ---
  k <- length(covariate_vars)
  if (k == 0L) stop("Pooled OLS requires at least one covariate in 'covariate_vars'.")
  # Need at least k+1 total rows to estimate (including intercept)
  if (nrow(pre_pool) <= k) {
    stop(sprintf("Pooled OLS: too few pre-treatment rows (%d) relative to regressors (%d).",
                 nrow(pre_pool), k))
  }

  # Ensure covariates exist and drop NA rows in covariates/outcome
  missing_cols <- setdiff(covariate_vars, names(pre_pool))
  if (length(missing_cols) > 0L) {
    stop("Pooled OLS: missing covariates in data: ", paste(missing_cols, collapse = ", "))
  }
  keep <- stats::complete.cases(pre_pool[, c(outcome_var, covariate_vars), drop = FALSE])
  pre_pool <- pre_pool[keep, , drop = FALSE]
  if (nrow(pre_pool) <= k) {
    stop("Pooled OLS: after removing NAs, not enough rows to estimate the model.")
  }

  # --- Fit pooled OLS on levels ---
  rhs <- paste(covariate_vars, collapse = " + ")
  fml <- stats::as.formula(paste(outcome_var, "~", rhs))
  fit <- stats::lm(fml, data = pre_pool)

  # Return named vector of coefficients for the covariates (exclude intercept)
  beta <- stats::coef(fit)
  beta <- beta[setdiff(names(beta), "(Intercept)")]
  # Keep only covariates explicitly requested (preserve order)
  beta <- beta[covariate_vars]
  return(beta)
}


#' Fit Pooled IV Model with Lagged Outcome (Homogeneous Beta)
#'
#' This function estimates a pooled IV model for covariate adjustment.
#' Pooled IV (AH-style) for common beta (pre-treatment windowed, per unit)
#' It instruments the lagged outcome with earlier lags and returns only
#' the beta coefficients for the included covariates.
#' Treats any regressor named "Y_lag1" as endogenous and instruments it with
#' lags Y_lag{min_iv_lag : max_iv_lag}. Other covariates are assumed exogenous.
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
#' @keywords internal
#'
# IV pooled betas for covariates (return ONLY covariate coefficients)

fit_common_beta_iv <- function(data,
                               outcome_var,
                               covariate_vars,
                               treat_time_var,
                               time_var,
                               unit_var,
                               min_iv_lag = 2,
                               max_iv_lag = 2,
                               degree = 0,
                               pretreatment_window = c("full", "minimal", "manual"),
                               pretreatment_years = NULL) {

  pretreatment_window <- match.arg(pretreatment_window)

  # --- Build a pooled pre-treatment frame using the per-unit window rule ---
  build_pre_for_unit <- function(df_u) {
    # Create requested outcome lags if needed (Y_lagk)
    # If the user included "Y_lag1" in covariate_vars, we also need to construct "Y_lag{2..K}"
    df_u <- df_u[order(df_u[[time_var]]), , drop = FALSE]
    if ("Y_lag1" %in% covariate_vars || max_iv_lag >= 2) {
      y <- df_u[[outcome_var]]
      for (L in 1:max_iv_lag) {
        df_u[[paste0("Y_lag", L)]] <- dplyr::lag(y, n = L)
        # lag is over the unit, since df_u is per-unit already
        y <- y # keep original y pointer
      }
    }

  # --- Define pre-treatment window (consistent with OLS version) ---
  df_u$timeToTreat <- df_u[[time_var]] - df_u[[treat_time_var]]
  pre_all <- df_u[df_u$timeToTreat < 0, , drop = FALSE]

  if (nrow(pre_all) == 0L) return (pre_all)

  if (pretreatment_window == "full") {
    pre_all[order(pre_all[[time_var]]), , drop = FALSE]

  } else if (pretreatment_window == "minimal") {
    need <- degree + 1L
    ord <- order(pre_all[[time_var]], decreasing = TRUE)
    take <- seq_len(min(need, length(ord)))
    idx <- sort(ord[take], decreasing = FALSE)
    pre_all[idx, , drop = FALSE]

  } else { # manual
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
    pre_all[idx, , drop = FALSE]
  }
  }

  units <- unique(data[[unit_var]])
  pre_list <- lapply(units, function(u) build_pre_for_unit(data[data[[unit_var]] == u, , drop = FALSE]))
  pre_pool <- do.call(rbind, pre_list)

  if (is.null(pre_pool) || nrow(pre_pool) == 0L) {
    stop("Pooled IV: no usable pre-treatment observations after windowing.")
  }

  # --- Define regressors and instruments ---
  if (length(covariate_vars) == 0L) {
    stop("Pooled IV requires at least one regressor in 'covariate_vars'.")
  }

  # Identify endogenous lag regressor(s) and construct instrument set
  endog_vars <- intersect(covariate_vars, "Y_lag1")  # currently only treat Y_lag1 as endogenous
  exog_vars  <- setdiff(covariate_vars, endog_vars)


  # Ensure instrument lags exist if needed
  inst_lags <- seq.int(min_iv_lag, max_iv_lag)
  inst_vars <- character(0L)
  if (length(endog_vars) > 0L) {
    # Need Y_lag{min_iv_lag:max_iv_lag} in the data
    need_cols <- paste0("Y_lag", inst_lags)
    missing_iv_cols <- setdiff(need_cols, names(pre_pool))
    if (length(missing_iv_cols) > 0L) {
      stop("Pooled IV: required instrument columns missing: ", paste(missing_iv_cols, collapse = ", "),
           ". Did you include the lagged outcome in 'covariate_vars' as 'Y_lag1'?")
    }
    inst_vars <- c(inst_vars, need_cols)
  }

  # Build model frames; drop rows with NAs in outcome, regressors, or instruments
  rhs_vars <- c(endog_vars, exog_vars)
  keep_cols <- unique(c(outcome_var, rhs_vars, inst_vars))
  keep <- stats::complete.cases(pre_pool[, keep_cols, drop = FALSE])
  pre_pool <- pre_pool[keep, , drop = FALSE]

  # Identifiability after NA-drop from lags: if too few rows remain, fall back to OLS
  k <- length(rhs_vars)

  # Build RHS and instrument strings (used in both branches below)
  rhs1 <- paste(rhs_vars, collapse = " + ")
  rhs2 <- paste(c(inst_vars, exog_vars), collapse = " + ")

  # Hard error if we have an endogenous regressor but literally no instruments constructed
  if (length(inst_vars) == 0L && length(endog_vars) > 0L) {
    stop("Pooled IV: endogenous regressor detected but no instruments were constructed.")
  }

  # If there are no endogenous vars at all, this path is just OLS
  if (length(endog_vars) == 0L) {
    # Ensure enough rows for OLS (intercept handled by lm)
    if (nrow(pre_pool) <= k) {
      stop(sprintf("Pooled OLS fallback: too few usable rows (%d) for %d regressors after NA removal.", nrow(pre_pool), k))
    }
    fml_ols <- stats::as.formula(paste(outcome_var, "~", rhs1))
    fit <- stats::lm(fml_ols, data = pre_pool)

  } else {
    # We do have endogenous vars and instruments available
    if (nrow(pre_pool) <= k) {
      # Too tight after lag-NA removal → graceful fallback to OLS
      warning(sprintf(
        "Pooled IV: too few usable rows (%d) for %d regressors after lag-NA removal; falling back to pooled OLS.",
        nrow(pre_pool), k
      ), call. = FALSE)
      fml_ols <- stats::as.formula(paste(outcome_var, "~", rhs1))
      fit <- stats::lm(fml_ols, data = pre_pool)
    } else {
      # Proceed with IV
      if (!requireNamespace("AER", quietly = TRUE)) {
        stop("Package 'AER' is required for ivreg(). Please install it.")
      }
      fml_iv <- stats::as.formula(paste0(outcome_var, " ~ ", rhs1, " | ", rhs2))
      fit <- AER::ivreg(fml_iv, data = pre_pool)
    }
  }

  # Extract β for the covariates only (drop intercept)
  beta <- stats::coef(fit)
  beta <- beta[setdiff(names(beta), "(Intercept)")]
  beta <- beta[covariate_vars]

  if (any(!is.finite(beta))) {
    stop("Pooled IV/OLS: non-finite coefficient(s). Widen pre window, add instruments, or reduce multicollinearity.")
  }


  return(beta)
}

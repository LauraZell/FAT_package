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
                                pretreatment_window = c("full", "minimal")) {
  pretreatment_window <- match.arg(pretreatment_window)

  # 1) Compute time_to_treat
  df <- dplyr::mutate(
    data,
    time_to_treat = .data[[time_var]] - .data[[treat_time_var]]
  )

  # 2) Keep strictly pre-treatment observations (< 0)
  if (pretreatment_window == "full") {
    df_pre <- dplyr::filter(df, time_to_treat < 0)
  } else { # "minimal": take exactly (degree + 1) most-recent pre rows per unit
    df_pre <- df |>
      dplyr::filter(time_to_treat < 0) |>
      dplyr::group_by(.data[[unit_var]]) |>
      dplyr::arrange(dplyr::desc(time_to_treat), .by_group = TRUE) |>
      dplyr::slice_head(n = degree + 1) |>
      dplyr::ungroup()
  }

  # 3) Basic checks
  if (!all(covariate_vars %in% names(df_pre))) {
    missing_covars <- setdiff(covariate_vars, names(df_pre))
    stop("Missing covariates in data: ", paste(missing_covars, collapse = ", "))
  }

  # Drop rows with NA in outcome or covariates
  df_pre <- df_pre |>
    dplyr::filter(!is.na(.data[[outcome_var]])) |>
    tidyr::drop_na(dplyr::all_of(covariate_vars))

  if (nrow(df_pre) == 0) {
    stop("No pre-treatment observations available after filtering (time_to_treat < 0).")
  }

  if (pretreatment_window == "minimal") {
    # Warn if some units contributed fewer than (degree+1) rows
    counts <- df_pre |>
      dplyr::count(.data[[unit_var]], name = "n_pre")
    if (any(counts$n_pre < (degree + 1))) {
      warning("Some units have fewer than (degree + 1) pre-treatment rows in 'minimal' mode.")
    }
  }

  # 4) Fit pooled OLS with only covariates
  cov_formula <- stats::as.formula(
    paste(outcome_var, "~", paste(covariate_vars, collapse = " + "))
  )
  model <- stats::lm(formula = cov_formula, data = df_pre)

  # 5) Return only the covariate coefficients, in the order requested
  bhat <- stats::coef(model)
  # Ensure we return exactly the covariates (intercept dropped even if present)
  bhat_covars <- bhat[covariate_vars]

  # If any are completely un-identified, they'll be NA; warn rather than error
  if (any(is.na(bhat_covars))) {
    na_names <- names(bhat_covars)[is.na(bhat_covars)]
    warning("Some pooled OLS covariate coefficients are NA (collinearity?): ",
            paste(na_names, collapse = ", "),
            ". They will propagate as NA into adjusted_outcome.")
  }

  return(bhat_covars)
}



#' Fit Pooled IV Model with Lagged Outcome (Homogeneous Beta)
#'
#' This function estimates a pooled IV model for covariate adjustment.
#' It instruments the lagged outcome with earlier lags and returns only
#' the beta coefficients for the included covariates (including the intercept).
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
                               pretreatment_window = c("full", "minimal")) {

  pretreatment_window <- match.arg(pretreatment_window)

  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("AER", quietly = TRUE)

  # Work on a copy
  df <- data

  # --- Create within-unit lags of the outcome (by unit, ordered by time) ---
  df <- dplyr::group_by(df, .data[[unit_var]]) |>
    dplyr::arrange(.data[[time_var]], .by_group = TRUE)

  # Make sure at least lag 1 exists for the endogenous regressor
  max_needed_lag <- max(1, max_iv_lag)
  for (l in 1:max_needed_lag) {
    lag_name <- paste0("lag", l, "_", outcome_var)
    df[[lag_name]] <- dplyr::lag(df[[outcome_var]], n = l)
  }
  df <- dplyr::ungroup(df)

  # Endogenous regressor: first lag of outcome
  df$lagged_outcome <- df[[paste0("lag1_", outcome_var)]]

  # --- Define pre-treatment window (consistent with OLS version) ---
  df <- dplyr::mutate(df, time_to_treat = .data[[time_var]] - .data[[treat_time_var]])

  if (pretreatment_window == "minimal") {
    # use exactly the last `degree` pre-treatment years:  time_to_treat in [-degree, -1]
    # if degree == 0, this yields an empty set; in that case, fall back to full (<0)
    if (degree > 0) {
      df <- dplyr::filter(df, time_to_treat < 0, time_to_treat >= -degree)
    } else {
      df <- dplyr::filter(df, time_to_treat < 0)
    }
  } else {
    # "full" = all pre-treatment years strictly before adoption
    df <- dplyr::filter(df, time_to_treat < 0)
  }

  # --- Build instrument set: lag(min_iv_lag) ... lag(max_iv_lag) of outcome
  instruments <- paste0("lag", min_iv_lag:max_iv_lag, "_", outcome_var)

  # Variables required to be non-missing
  needed <- c("lagged_outcome", instruments, covariate_vars, outcome_var)

  # Keep complete cases
  if (!all(needed %in% names(df))) {
    missing_cols <- setdiff(needed, names(df))
    stop("IV beta: missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  df <- df[stats::complete.cases(df[, needed]), ]

  # Guardrails: enough rows?
  if (nrow(df) == 0) {
    stop("IV beta: 0 (non-NA) cases after filtering to pre-treatment window and complete cases.")
  }

  # --- Build IV formula ---
  # Structural RHS includes intercept + lagged_outcome + covariates
  rhs_str     <- paste(c("1", "lagged_outcome", covariate_vars), collapse = " + ")
  rhs_instr   <- paste(c(instruments, covariate_vars), collapse = " + ")
  iv_formula  <- stats::as.formula(paste0(outcome_var, " ~ ", rhs_str, " | ", rhs_instr))

  # --- Estimate IV ---
  iv_model <- AER::ivreg(iv_formula, data = df)

  # Return ONLY covariate coefficients (no intercept, no endogenous regressor)
  bhat <- stats::coef(iv_model)
  keep <- covariate_vars
  bhat_cov <- bhat[keep]

  # Safety check: names + length
  if (any(is.na(bhat_cov)) || length(bhat_cov) != length(covariate_vars)) {
    stop("IV beta: could not recover covariate coefficients cleanly. ",
         "Check instrument strength / sample size.")
  }

  return(bhat_cov)
}

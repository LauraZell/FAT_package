#' Estimate Forecasted Average Treatment Effects (FAT / DFAT)
#'
#' @description
#' Estimates forecasted average treatment effects (FAT) or difference-in-FAT
#' (DFAT) across a range of polynomial degrees and forecast horizons. Optionally
#' returns unit-level prediction trajectories for diagnostic and plotting use.
#'
#' @param data A long-format dataframe with panel structure.
#' @param unit_var Name of the column for the unit ID (e.g., "state", "mun_id").
#' @param time_var Name of the column for the time variable (e.g., "year").
#' @param outcome_var Name of the outcome variable to be forecasted. (e.g., "ln_age_mort_rate").
#' @param treat_time_var Name of the column indicating treatment year per unit. (e.g., "adopt_year").
#' @param units_to_include Optional vector of unit names to restrict analysis to (e.g., "state1", "state2").
#' @param degrees Vector of polynomial degrees for trend fitting (e.g., 0:2).
#' @param horizons Vector of forecast horizons to estimate (e.g., 1:5).
#' @param se_method Method for standard errors: "analytic", "bootstrap", "clustered", "unitwise".
#' Note: when dynamic parameters are estimated via pooled IV (e.g., beta_estimator = "iv"
#'   with lagged outcome among covariates), analytic SEs are not available and will be
#'   automatically switched to "bootstrap" with a warning.
#' @param n_bootstrap Number of bootstrap repetitions (only if `se_method` is "bootstrap").
#' @param covariate_vars Optional vector of column names to include as covariates.
#' @param beta_estimator Pooled beta estimation method: "none", "ols", "iv", or "unitwise".
#' @param min_iv_lag Minimum lag for IV estimation (if used).
#' @param max_iv_lag Maximum lag for IV estimation (if used).
#' @param control_group_value Optional. If set (e.g., control_group_value = FALSE), DFAT mode is activated.
#'                            Treated units are expected to be marked in a "treated" column (TRUE/FALSE).
#' @param forecast_offset Number of periods to wait after treatment before forecasting begins.
#'        Default is 0 (forecast starts in treatment year).
#' @param pretreatment_window Character. Defines how much pre-treatment data to use for trend fitting:
#'                            "full" uses all pre-treatment years (`timeToTreat < 0`);
#'                            "minimal" uses only the minimum number of years required by the polynomial degree.
#'                            Default is "full".
#'
#' @return A list with:
#' \describe{
#'   \item{results}{Data frame with summary estimates: deg, hh, FAT, sdFAT.}
#'   \item{predictions}{Data frame with unit-level observed and predicted values across time, for plotting.}
#' }
#' @export


estimate_fat <- function(data,
                         unit_var,
                         time_var,
                         outcome_var,
                         treat_time_var,
                         units_to_include = NULL,
                         degrees = 0:2,
                         horizons = 1:5,
                         se_method = "analytic",
                         n_bootstrap = 1000,
                         covariate_vars = NULL,
                         beta_estimator = c("none", "ols", "iv", "unitwise"),
                         min_iv_lag = 2,
                         max_iv_lag = 2,
                         control_group_value = NULL,
                         forecast_offset = 1,
                         pretreatment_window = c("full", "minimal", "manual"),
                         pretreatment_years = NULL) {
  pretreatment_window <- match.arg(pretreatment_window)
  beta_estimator <- match.arg(beta_estimator)

  # SE guardrails for model-based / IV paths
  # If pooled IV is used and the lagged outcome is among covariates (e.g., "Y_lag1"),
  # we treat this as a "dynamic parameter estimated" case (rho-like). In this case,
  # analytic SE for FAT is not available; we force bootstrap with a warning.
  is_dynamic_iv <- (beta_estimator == "iv") &&
    !is.null(covariate_vars) &&
    any(covariate_vars %in% c("Y_lag1"))

  if (is_dynamic_iv && identical(se_method, "analytic")) {
    warning(
      paste0(
        "se_method='analytic' is not available when beta_estimator='iv' ",
        "with lagged outcome among covariates (dynamic parameter estimated). ",
        "Switching to se_method='bootstrap'."
      ),
      call. = FALSE
    )
    se_method <- "bootstrap"
    # bump bootstrap reps if user left the default very low
    if (is.null(n_bootstrap) || !is.finite(n_bootstrap) || n_bootstrap < 199L) {
      n_bootstrap <- max(499L, as.integer(n_bootstrap %||% 0L))
    }
  }

  # Also prevent "unitwise" + analytic in DFAT with tiny N (optional, conservative):
  if (!is.null(control_group_value) && beta_estimator == "unitwise" && identical(se_method, "analytic")) {
    # keep analytic allowed, but warn if too few units (no robust CLT)
    n_units <- dplyr::n_distinct(data[[unit_var]])
    if (n_units < 25L) {
      warning("DFAT + unitwise + se_method='analytic' with <25 units may be anti-conservative. Consider 'bootstrap'.", call. = FALSE)
    }
  }


  # Validate manual mode
  if (pretreatment_window == "manual") {
    if (is.null(pretreatment_years) || length(pretreatment_years) != 1L || !is.finite(pretreatment_years)) {
      stop("When pretreatment_window = 'manual', provide a single numeric 'pretreatment_years' value.")
    }
    if (pretreatment_years < 1L) {
      stop("'pretreatment_years' must be at least 1.")
    }
  }

  # Argument validation for covariates vs estimator
  # Ensure covariate_vars is a proper character vector when needed
  is_missing_covars <- is.null(covariate_vars) || length(covariate_vars) == 0

  if (beta_estimator == "unitwise") {
    if (is_missing_covars) {
      stop(
        "beta_estimator = 'unitwise' requires 'covariate_vars' to be provided.\n",
        "Hint: pass a character vector of column names, e.g., covariate_vars = c('x1','x2')."
      )
    }
    # check existence
    missing_cols <- setdiff(covariate_vars, names(data))
    if (length(missing_cols) > 0) {
      stop(
        "The following 'covariate_vars' are not in 'data': ",
        paste(missing_cols, collapse = ", "), "\n",
        "Please ensure all covariates exist prior to calling estimate_fat()."
      )
    }
  }

  # (Recommended) enforce covariates for pooled estimators too
  if (beta_estimator %in% c("ols", "iv")) {
    if (is_missing_covars) {
      stop(
        "beta_estimator = '", beta_estimator,
        "' requires 'covariate_vars'. Provide e.g. covariate_vars = c('x1','x2')."
      )
    }
    missing_cols <- setdiff(covariate_vars, names(data))
    if (length(missing_cols) > 0) {
      stop(
        "The following 'covariate_vars' are not in 'data': ",
        paste(missing_cols, collapse = ", "), "\n",
        "Please ensure all covariates exist prior to calling estimate_fat()."
      )
    }
  }


  # Optional subsetting to user-specified units
  if (!is.null(units_to_include)) {
    data <- data[data[[unit_var]] %in% units_to_include, ]
  }

  # Flag if we are in DFAT (difference-in-FAT) mode
  dfat_mode <- !is.null(control_group_value)

  if (dfat_mode && !"treated" %in% names(data)) {
    stop("DFAT mode requires a column named 'treated' in the dataset.")
  }

  # Set common treatment time for DFAT (median among treated) or use individual treat_time for regular FAT
  if (dfat_mode) {
    # Use median treatment year among treated units as a common reference
    treated_years <- data[data$treated != control_group_value & !is.na(data[[treat_time_var]]), treat_time_var]
    if (length(treated_years) == 0) stop("No treated units with valid treatment year.")
    reference_treat_time <- stats::median(treated_years, na.rm = TRUE)
    data$treat_time_for_fit <- reference_treat_time
  } else {
    data$treat_time_for_fit <- data[[treat_time_var]]
  }
  # Create time-to-treatment variable (centered time)
  data <- dplyr::mutate(data, timeToTreat = as.numeric(.data[[time_var]]) - .data$treat_time_for_fit)


  # Helper: choose which pre-treatment rows to use for a given unit & degree
  .select_pre_rows <- function(df_unit, deg, window, k_years) {
    pre <- df_unit[df_unit$timeToTreat < 0 & !is.na(df_unit[[outcome_var]]), , drop = FALSE]
    if (nrow(pre) == 0L) return(pre)

    if (window == "full") {
      return(pre[order(pre[[time_var]]), , drop = FALSE])
    }
    if (window == "minimal") {
      need <- deg + 1L
      idx <- order(pre[[time_var]], decreasing = TRUE)[seq_len(min(need, nrow(pre)))]
      return(pre[sort(idx, decreasing = FALSE), , drop = FALSE])
    }
    # manual
    need <- deg + 1L
    if (k_years < need) {
      stop(sprintf("pretreatment_years = %d is too small for degree q = %d (need at least q+1 = %d).",
                   k_years, deg, need))
    }
    take <- min(k_years, nrow(pre))
    idx <- order(pre[[time_var]], decreasing = TRUE)[seq_len(take)]
    pre[sort(idx, decreasing = FALSE), , drop = FALSE]
  }

  # Internal function to estimate FAT/DFAT with specified degree and horizon
  fat_for_combo <- function(deg, hh) {

    # Saturation/identifiability check per unit, using chosen pretreatment window
    required_params <- deg + 1L

    pre_counts <- data %>%
      dplyr::group_by(.data[[unit_var]]) %>%
      dplyr::group_modify(~{
        pre_use <- .select_pre_rows(.x, deg, pretreatment_window, pretreatment_years)
        tibble::tibble(n_pre = nrow(pre_use))
      }) %>%
      dplyr::ungroup() %>%
      dplyr::rename(.unit = !!rlang::sym(unit_var))

    if (nrow(pre_counts) == 0L) {
      stop("No units present after initial filtering.")
    }

    too_few <- pre_counts %>% dplyr::filter(n_pre < required_params)
    if (nrow(too_few) > 0L) {
      ex <- paste(utils::head(too_few$.unit, 5L), collapese = ", ")
      stop(sprintf(
        "Cannot fit q = %d: some units have fewer than q+1 = %d usable pretreatment rows under '%s' window. Examples: %s",
        deg, required_params, pretreatment_window, ex))
    }

    exactly_equal <- pre_counts %>% dplyr::filter(n_pre == required_params)
    if (nrow(exactly_equal) > 0L) {
      ex <- paste(utils::head(exactly_equal$.unit, 5L), collapse = ", ")
      warning(sprintf(
        "Saturated at q = %d for some units: usable pretreatment rows equal q + 1 = %d under '%s'. Forecasts may be unstable. Examples: %s",
        deg, required_params, pretreatment_window, ex), call. = FALSE)
    }

    # If covariates + pooled estimation selected, compute beta
    beta_hat <- NULL
    if (!is.null(covariate_vars) && beta_estimator != "none") {
      if (beta_estimator == "ols") {
        beta_hat <- fit_common_beta_ols(data = data,
                                        outcome_var = outcome_var,
                                        covariate_vars = covariate_vars,
                                        treat_time_var = "treat_time_for_fit",
                                        time_var = time_var,
                                        unit_var = unit_var,
                                        degree = deg,
                                        pretreatment_window = pretreatment_window,
                                        pretreatment_years = pretreatment_years)

      } else if (beta_estimator == "iv") {
        beta_hat <- fit_common_beta_iv(data          = data,
                                       outcome_var   = outcome_var,
                                       covariate_vars= covariate_vars,
                                       treat_time_var= "treat_time_for_fit",
                                       time_var      = time_var,
                                       unit_var      = unit_var,
                                       degree        = deg,
                                       min_iv_lag    = min_iv_lag,
                                       max_iv_lag    = max_iv_lag,
                                       pretreatment_window = pretreatment_window,
                                       pretreatment_years = pretreatment_years)
      }
    }

    # Unitwise trend estimation requires at least a linear polynomial
    if (beta_estimator == "unitwise" && deg < 1) {
      return(list(summary = data.frame(deg = deg, hh = hh, FAT = NA_real_, sdFAT = NA_real_),
                  preds = NULL))
    }

    # Forecast outcome for each unit
    all_preds <- unique(data[[unit_var]]) %>%
      purrr::map_df(~ {
        if (beta_estimator == "unitwise") {
          fit_unitwise_trend(
            data = data,
            unit = .x,
            degree = deg,
            unit_var = unit_var,
            time_var = time_var,
            outcome_var = outcome_var,
            treat_time_var = "treat_time_for_fit",
            covariate_vars = covariate_vars,
            beta_hat = NULL,
            forecast_offset = forecast_offset,
            pretreatment_window = pretreatment_window,
            pretreatment_years = pretreatment_years,
            hh = hh
          )
        } else {
          fit_common_trend(
            data = data,
            unit = .x,
            degree = deg,
            unit_var = unit_var,
            time_var = time_var,
            outcome_var = outcome_var,
            treat_time_var = "treat_time_for_fit",
            covariate_vars = covariate_vars,
            beta_hat = beta_hat,
            forecast_offset = forecast_offset,
            pretreatment_window = pretreatment_window,
            pretreatment_years = pretreatment_years,
            hh = hh
          )
        }
      })

    # Set up to enforce uniqueness of predictions per (unit, time, deg, hh)
    key_unit <- rlang::sym(unit_var)
    key_time <- rlang::sym(time_var)


    # Attach deg and hh to savely use them in grouping/filters
    all_preds <- all_preds %>%
      dplyr::mutate(deg = deg, hh = hh) %>%
      dplyr::distinct(!!key_unit, !!key_time, deg, hh, .keep_all = TRUE)

    # For DFAT: Merge treatment indicator back (used only in post-treatment comparison) (new: include guardrails against duplicates)
          # if (dfat_mode) {
          #   all_preds <- dplyr::left_join(
          #     all_preds,
          #     dplyr::select(data, .data[[unit_var]], .data[[time_var]], treated),
          #     by = c(unit_var, time_var)
          #   )
          # }


    if (dfat_mode) {
      treated_lookup <- data %>%
        dplyr::select(!!key_unit, !!key_time, treated) %>%
        dplyr::distinct(!!key_unit, !!key_time, .keep_all = TRUE)

      all_preds <- dplyr::left_join(
        all_preds,
        treated_lookup,
        by = c(unit_var, time_var)
      )
    }

    # Keep only needed columns:
    sel_cols <- c(unit_var, time_var, outcome_var,
                  "preds", "timeToTreat", "hh", "deg", "n_pre_fit", "pre_years_used")
    if (dfat_mode) sel_cols <- c(sel_cols, "treated")

    all_preds <- all_preds %>%
      dplyr::select(dplyr::all_of(sel_cols)) %>%
      # belt & suspenders: re‑enforce uniqueness of the key
      dplyr::distinct(!!key_unit, !!key_time, deg, hh, .keep_all = TRUE)

    # Remove later: one-line diagnostic
    dup_check <- all_preds %>%
      dplyr::count(!!key_unit, !!key_time, deg, hh, name = "n") %>%
      dplyr::filter(n > 1)

    if (nrow(dup_check) > 0) {
      warning("Internal duplicate keys survived to all_preds; collapsing to first row per key.\n",
              "Examples: ",
              paste(utils::capture.output(print(utils::head(dup_check, 5))), collapse = "\n"))
    }

    # Target rows for FAT at the specified horizon only
    # Pick the single forecast step 'hh' (i.e., timeToTreat == forecast_offset + hh - 1)
    .h <- hh  # avoid NSE name collision
    target_data <- all_preds %>%
      dplyr::filter(.data$hh == .h) %>%
      dplyr::mutate(diff = .data[[outcome_var]] - preds)


    # ===================== DFAT logic =====================
    if (dfat_mode) {
      treat_diff <- dplyr::filter(target_data, treated != control_group_value)$diff
      control_diff <- dplyr::filter(target_data, treated == control_group_value)$diff

      FAT <- mean(treat_diff, na.rm = TRUE) - mean(control_diff, na.rm = TRUE)

      if (se_method == "analytic") {
        sdFAT <- sqrt(stats::var(treat_diff, na.rm = TRUE) / length(treat_diff) +
                        stats::var(control_diff, na.rm = TRUE) / length(control_diff))
      } else if (se_method == "bootstrap") {
        FAT_boot <- replicate(n_bootstrap, {
          sampled_treat <- sample(treat_diff, replace = TRUE)
          sampled_control <- sample(control_diff, replace = TRUE)
          mean(sampled_treat, na.rm = TRUE) - mean(sampled_control, na.rm = TRUE)
        })
        sdFAT <- stats::sd(FAT_boot, na.rm = TRUE)
      } else if (se_method == "clustered") {
        target_data$treated_dummy <- as.integer(target_data$treated != control_group_value)
        mod <- stats::lm(diff ~ treated_dummy, data = target_data)
        cluster <- target_data[[unit_var]]
        sdFAT <- sqrt(sandwich::vcovCL(mod, cluster = cluster)["treated_dummy", "treated_dummy"])
      } else {
        stop("Unsupported SE method in DFAT mode.")
      }

    } else {
      # ===================== Regular FAT =====================
      FAT <- mean(target_data$diff, na.rm = TRUE)

      if (se_method == "analytic") {
        sdFAT <- (1 / sqrt(nrow(target_data))) * sqrt(mean((target_data$diff - FAT)^2, na.rm = TRUE))
      } else if (se_method == "bootstrap") {
        units <- unique(target_data[[unit_var]])
        FAT_boot <- replicate(n_bootstrap, {
          sampled_units <- sample(units, replace = TRUE)
          sampled_data <- target_data[target_data[[unit_var]] %in% sampled_units, ]
          mean(sampled_data$diff, na.rm = TRUE)
        })
        sdFAT <- stats::sd(FAT_boot, na.rm = TRUE)
      } else if (se_method == "clustered") {
        mod <- stats::lm(diff ~ 1, data = target_data)
        cluster <- target_data[[unit_var]]
        sdFAT <- sqrt(sandwich::vcovCL(mod, cluster = cluster)[1, 1])
      } else if (se_method == "unitwise") {
        unit_fats <- target_data |>
          dplyr::group_by(.data[[unit_var]]) |>
          dplyr::summarise(diff = mean(diff, na.rm = TRUE), .groups = "drop")
        sdFAT <- stats::sd(unit_fats$diff, na.rm = TRUE)
      } else {
        stop("Invalid SE method for regular FAT.")
      }
    }

    # Return both summary and per-unit prediction data
    return(list(
      summary = data.frame(deg = deg, hh = hh, FAT = FAT, sdFAT = sdFAT),
      predictions = all_preds
    ))

  }

  # Generate all (degree, horizon) combinations
  combos <- base::expand.grid(deg = degrees, hh = horizons)
  results_list <- purrr::pmap(combos, fat_for_combo)
  summary_df <- purrr::map_dfr(results_list, "summary")
  preds_df <- purrr::map_dfr(results_list, "predictions")

  # belt-and-suspenders: keep a single row per (unit, time, deg, hh)
  preds_df <- preds_df %>%
    dplyr::distinct(!!rlang::sym(unit_var), !!rlang::sym(time_var), deg, hh, .keep_all = TRUE)


  return(list(results = summary_df, predictions = preds_df))
}



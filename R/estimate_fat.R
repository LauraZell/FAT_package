# =============================================================================
# estimate_fat: Forecasted Average Treatment (FAT) and Difference-in-FAT (DFAT)
#               estimator with optional standard errors and return of forecasts
#
# INPUTS:
# - data: a long-format dataframe with panel structure
# - unit_var: column name for unit ID (e.g., state, municipality)
# - time_var: column name for time (e.g., year)
# - outcome_var: column name for observed outcome
# - treat_time_var: column name for treatment timing
# - units_to_include: optional vector of units to subset on
# - degrees: polynomial degrees to use for trend fitting
# - horizons: forecast horizons to estimate (e.g. 1 to 5 years after treatment)
# - se_method: standard error method ("analytic", "bootstrap", "clustered", etc.)
# - n_bootstrap: number of bootstrap reps if bootstrap SEs used
# - covariate_vars: optional covariates to include in trend fitting
# - beta_estimator: pooled beta estimation method ("none", "ols", "iv", "unitwise")
# - min_iv_lag, max_iv_lag: used if beta_estimator = "iv"
# - control_group_value: if not NULL, activates DFAT mode using "treated" column
#
# RETURNS:
# - A list with:
#     $results      => summary results by (deg, hh)
#     $predictions  => unit-level forecasts and outcomes for plotting
# =============================================================================

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
                         control_group_value = NULL) {

  # Match estimator option
  beta_estimator <- match.arg(beta_estimator)

  # Optional subsetting of units
  if (!is.null(units_to_include)) {
    data <- data[data[[unit_var]] %in% units_to_include, ]
  }

  # Flag if we are in DFAT (difference-in-FAT) mode
  dfat_mode <- !is.null(control_group_value)

  if (dfat_mode && !"treated" %in% names(data)) {
    stop("DFAT mode requires a column named 'treated' in the dataset.")
  }

  # Define the treatment time to use for fitting
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

  # Internal function to estimate FAT/DFAT for given (degree, horizon)
  fat_for_combo <- function(deg, hh) {
    beta_hat <- NULL

    # If covariates + pooled estimation selected, compute beta
    if (!is.null(covariate_vars) && beta_estimator != "none") {
      if (beta_estimator == "ols") {
        beta_hat <- fit_common_beta_ols(data, outcome_var, covariate_vars, time_var,
                                        "treat_time_for_fit", unit_var, deg)
      } else if (beta_estimator == "iv") {
        beta_hat <- fit_common_beta_iv(data, outcome_var, covariate_vars, time_var,
                                       "treat_time_for_fit", unit_var, deg, min_iv_lag, max_iv_lag)
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
          fit_unitwise_trend(data, .x, deg, unit_var, time_var, outcome_var,
                             "treat_time_for_fit", covariate_vars)
        } else {
          fit_common_trend(data, .x, deg, unit_var, time_var, outcome_var,
                           "treat_time_for_fit", covariate_vars, beta_hat)
        }
      })

    # Add treated status in DFAT mode
    if (dfat_mode) {
      all_preds <- dplyr::left_join(
        all_preds,
        dplyr::select(data, .data[[unit_var]], .data[[time_var]], treated),
        by = c(unit_var, time_var)
      )
    }

    # Add meta info for faceting/plotting
    all_preds$deg <- deg
    all_preds$hh <- hh

    # Compute outcome difference at the treatment + hh year
    target_data <- dplyr::filter(all_preds, .data[[time_var]] == (.data$treat_time_for_fit + hh)) %>%
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
    list(summary = data.frame(deg = deg, hh = hh, FAT = FAT, sdFAT = sdFAT),
         preds = all_preds)
  }

  # Generate all (degree, horizon) combinations
  combos <- base::expand.grid(deg = degrees, hh = horizons)

  # Run main function over all combos
  results_list <- purrr::pmap(combos, ~ fat_for_combo(..1, ..2))

  # Separate summary and predictions
  summary_df <- purrr::map_dfr(results_list, "summary")
  preds_list <- purrr::map(results_list, "preds")
  all_predictions <- dplyr::bind_rows(preds_list, .id = "combo_id")

  return(list(results = summary_df, predictions = all_predictions))
}

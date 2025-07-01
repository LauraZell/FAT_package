#' Estimate Forecasted Average Treatment (FAT) or Differential FAT (DFAT)
#'
#' This function estimates forecasted average treatment effects (FAT) or differential
#' forecasted effects (DFAT) over a grid of polynomial degrees and post-treatment horizons.
#'
#' @param data Panel dataset with unit, time, outcome, and treatment variables.
#' @param unit_var Name of the unit variable.
#' @param time_var Name of the time variable.
#' @param outcome_var Name of the outcome variable.
#' @param treat_time_var Name of the treatment timing variable.
#' @param units_to_include Optional vector of units to include in estimation.
#' @param degrees Vector of polynomial degrees for fitting trends.
#' @param horizons Vector of post-treatment horizons to estimate FATs.
#' @param se_method Method for computing standard errors: "analytic", "bootstrap", "clustered", "unitwise".
#' @param n_bootstrap Number of bootstrap samples for "bootstrap" SEs.
#' @param covariate_vars Optional covariates to partial out.
#' @param beta_estimator One of "none", "ols", "iv", "unitwise".
#' @param min_iv_lag Minimum lag to use as instrument (for "iv" model).
#' @param max_iv_lag Maximum lag to use as instrument (for "iv" model).
#' @param control_group_value Optional value (e.g., FALSE) that identifies the control group in a column called "treated".
#'
#' @return A data.frame with deg, hh, FAT or DFAT, sdFAT
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
                         control_group_value = NULL) {

  beta_estimator <- match.arg(beta_estimator)

  # Optionally restrict to selected units
  if (!is.null(units_to_include)) {
    data <- data[data[[unit_var]] %in% units_to_include, ]
  }

  # Check if user wants to run DFAT mode by specifying a control group value
  dfat_mode <- !is.null(control_group_value)

  # In DFAT mode, a column called 'treated' must exist
  if (dfat_mode && !"treated" %in% names(data)) {
    stop("If 'control_group_value' is specified, a column named 'treated' must exist.")
  }

  # In DFAT mode, define a common pseudo-treatment year for all units
  if (dfat_mode) {
    # This step computes the reference treatment year as the median adoption year across treated units
    treated_years <- data[data$treated != control_group_value & !is.na(data[[treat_time_var]]), treat_time_var]
    if (length(treated_years) == 0) stop("No treated units with valid treatment year.")
    reference_treat_time <- stats::median(treated_years, na.rm = TRUE)

    # Replace original treatment time for forecasting
    data$treat_time_for_fit <- reference_treat_time
  } else {
    # In regular FAT mode, use unit-specific treatment times
    data$treat_time_for_fit <- data[[treat_time_var]]
  }

  # Compute relative time to treatment for each observation
  data <- dplyr::mutate(data, timeToTreat = .data[[time_var]] - .data$treat_time_for_fit)

  # This is the core estimation function used in the grid expansion (deg, horizon)
  fat_for_combo <- function(deg, hh) {
    beta_hat <- NULL

    # Compute beta_hat if a model is specified
    if (!is.null(covariate_vars) && beta_estimator != "none") {
      if (beta_estimator == "ols") {
        beta_hat <- fit_common_beta_ols(
          data = data,
          outcome_var = outcome_var,
          covariate_vars = covariate_vars,
          time_var = time_var,
          treat_time_var = "treat_time_for_fit",
          unit_var = unit_var,
          degree = deg
        )
      } else if (beta_estimator == "iv") {
        beta_hat <- fit_common_beta_iv(
          data = data,
          outcome_var = outcome_var,
          covariate_vars = covariate_vars,
          time_var = time_var,
          treat_time_var = "treat_time_for_fit",
          unit_var = unit_var,
          degree = deg,
          min_iv_lag = min_iv_lag,
          max_iv_lag = max_iv_lag
        )
      }
    }

    # Skip trend fitting in unitwise model if degree is zero
    if (beta_estimator == "unitwise" && deg < 1) {
      return(data.frame(deg = deg, hh = hh, FAT = NA_real_, sdFAT = NA_real_))
    }

    # Estimate trends and predictions for each unit
    all_preds <- unique(data[[unit_var]]) %>%
      purrr::map_df(~ {
        if (beta_estimator == "unitwise") {
          fit_unitwise_trend(
            data = data,
            unit = .x,
            deg = deg,
            unit_var = unit_var,
            time_var = time_var,
            outcome_var = outcome_var,
            treat_time_var = "treat_time_for_fit",
            covariate_vars = covariate_vars
          )
        } else {
          fit_common_trend(
            data = data,
            unit = .x,
            deg = deg,
            unit_var = unit_var,
            time_var = time_var,
            outcome_var = outcome_var,
            treat_time_var = "treat_time_for_fit",
            covariate_vars = covariate_vars,
            beta_hat = beta_hat
          )
        }
      })

    # --- Fix: join "treated" status into all_preds manually ---
    # This fixes a common source of bugs when the treated column is not preserved in prediction data
    all_preds <- dplyr::left_join(
      all_preds,
      dplyr::select(data, .data[[unit_var]], .data[[time_var]], treated),
      by = c(unit_var, time_var)
    )
    # Future improvement: consider ensuring treated is always part of fit_* return if needed

    # Filter to post-treatment target horizon
    target_data <- dplyr::filter(all_preds, .data[[time_var]] == (.data$treat_time_for_fit + hh)) %>%
      dplyr::mutate(diff = .data[[outcome_var]] - preds)

    # ---------- DFAT mode ----------
    if (dfat_mode) {
      # In DFAT mode, compute difference in mean forecast error between treated and control units
      treat_diff <- dplyr::filter(target_data, treated != control_group_value)$diff
      control_diff <- dplyr::filter(target_data, treated == control_group_value)$diff

      FAT <- mean(treat_diff, na.rm = TRUE) - mean(control_diff, na.rm = TRUE)

      if (se_method == "analytic") {
        sdFAT <- sqrt(stats::var(treat_diff, na.rm = TRUE) / length(treat_diff) +
                        stats::var(control_diff, na.rm = TRUE) / length(control_diff))
      } else {
        stop("Only analytic SE supported for DFAT at the moment.")
      }

    } else {
      # ---------- Regular FAT mode ----------
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
        cluster_se <- sqrt(sandwich::vcovCL(mod, cluster = cluster)[1, 1])
        sdFAT <- cluster_se
      } else if (se_method == "unitwise") {
        unit_fats <- target_data |>
          dplyr::group_by(.data[[unit_var]]) |>
          dplyr::summarise(diff = mean(diff, na.rm = TRUE), .groups = "drop")
        sdFAT <- stats::sd(unit_fats$diff, na.rm = TRUE)
      } else {
        stop("Invalid se_method.")
      }
    }

    # Return summary result
    data.frame(deg = deg, hh = hh, FAT = FAT, sdFAT = sdFAT)
  }

  # Loop over all degree-horizon combos
  combos <- base::expand.grid(deg = degrees, hh = horizons)
  purrr::pmap_dfr(combos, ~ fat_for_combo(..1, ..2))
}

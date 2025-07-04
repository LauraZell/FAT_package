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

  # Subset data if user provides a unit subset
  if (!is.null(units_to_include)) {
    data <- data[data[[unit_var]] %in% units_to_include, ]
  }

  # ========== DFAT MODE: Logic for difference-in-forecast-errors ==========
  dfat_mode <- !is.null(control_group_value)

  # Design decision: We do not allow the user to specify a treated column name.
  # The logic expects "treated" column to exist if DFAT mode is activated.
  if (dfat_mode && !"treated" %in% names(data)) {
    stop("DFAT mode requires a column named 'treated' in the dataset.")
  }

  # Compute treatment reference time
  if (dfat_mode) {
    # Reference treatment time for fitting: median treat year among treated units
    treated_years <- data[data$treated != control_group_value & !is.na(data[[treat_time_var]]), treat_time_var]
    if (length(treated_years) == 0) stop("No treated units with valid treatment year.")
    reference_treat_time <- stats::median(treated_years, na.rm = TRUE)

    # Use reference treatment time for fitting all units
    data$treat_time_for_fit <- reference_treat_time
  } else {
    # Regular FAT mode: fit per-unit adoption years
    data$treat_time_for_fit <- data[[treat_time_var]]
  }

  # Create time-to-treatment variable
  data <- dplyr::mutate(data, timeToTreat = as.numeric(.data[[time_var]]) - .data$treat_time_for_fit)

  # ========== Define internal combo function ==========

  fat_for_combo <- function(deg, hh) {
    beta_hat <- NULL

    # Optional pooled regression step if covariates included
    if (!is.null(covariate_vars) && beta_estimator != "none") {
      if (beta_estimator == "ols") {
        beta_hat <- fit_common_beta_ols(data, outcome_var, covariate_vars, time_var,
                                        "treat_time_for_fit", unit_var, deg)
      } else if (beta_estimator == "iv") {
        beta_hat <- fit_common_beta_iv(data, outcome_var, covariate_vars, time_var,
                                       "treat_time_for_fit", unit_var, deg, min_iv_lag, max_iv_lag)
      }
    }

    # Defensive: unitwise needs a polynomial trend
    if (beta_estimator == "unitwise" && deg < 1) {
      return(data.frame(deg = deg, hh = hh, FAT = NA_real_, sdFAT = NA_real_))
    }

    # Fit trend and predict for each unit
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

    # FIX: In DFAT mode, merge back "treated" column from original dataset
    if (dfat_mode) {
      all_preds <- dplyr::left_join(
        all_preds,
        dplyr::select(data, .data[[unit_var]], .data[[time_var]], treated),
        by = c(unit_var, time_var)
      )
    }

    # Evaluate prediction horizon: 1, 2, ... years after adoption
    target_data <- dplyr::filter(all_preds, .data[[time_var]] == (.data$treat_time_for_fit + hh)) %>%
      dplyr::mutate(diff = .data[[outcome_var]] - preds)

    # ==================== DFAT LOGIC ====================
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
        # Clustered SEs: use diff ~ treatment dummy
        target_data$treated_dummy <- as.integer(target_data$treated != control_group_value)
        mod <- stats::lm(diff ~ treated_dummy, data = target_data)
        cluster <- target_data[[unit_var]]
        cluster_se <- sqrt(sandwich::vcovCL(mod, cluster = cluster)["treated_dummy", "treated_dummy"])
        sdFAT <- cluster_se
      } else {
        stop("Unsupported SE method in DFAT mode.")
      }

    } else {
      # ==================== Regular FAT ====================
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

    data.frame(deg = deg, hh = hh, FAT = FAT, sdFAT = sdFAT)
  }

  # Generate all combinations of degrees × horizons
  combos <- base::expand.grid(deg = degrees, hh = horizons)

  # Run estimation
  purrr::pmap_dfr(combos, ~ fat_for_combo(..1, ..2))
}

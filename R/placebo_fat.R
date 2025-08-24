#' Placebo FAT / DFAT (wrapper around estimate_fat)
#'
#' @description
#' Runs FAT or DFAT after pretending treatment happened earlier by each placebo lag.
#' Fully reuses \code{estimate_fat()} so all current package options/logic apply.
#'
#' @param data A long-format dataframe with panel structure.
#' @param unit_var Name of the column for the unit ID (e.g., "state", "mun_id").
#' @param time_var Name of the column for the time variable (e.g., "year").
#' @param outcome_var Name of the outcome variable to be forecasted.
#' @param treat_time_var Name of the column indicating treatment year per unit.
#' @param units_to_include Optional vector of unit names to restrict analysis.
#' @param degrees Vector of polynomial degrees for trend fitting.
#' @param horizons Vector of forecast horizons to estimate.
#' @param lags Integer vector of placebo lags (years moved earlier), e.g. 0:2 or 1:3.
#' @param se_method Method for standard errors: "analytic", "bootstrap", "clustered", "unitwise".
#' @param n_bootstrap Number of bootstrap repetitions (only if `se_method` == "bootstrap").
#' @param covariate_vars Optional vector of column names to include as covariates.
#' @param beta_estimator Pooled beta estimation method: "none", "ols", "iv", or "unitwise".
#' @param min_iv_lag Minimum lag for IV estimation (if used).
#' @param max_iv_lag Maximum lag for IV estimation (if used).
#' @param control_group_value Optional. If set (e.g., FALSE), DFAT mode is activated (requires `treated`).
#' @param forecast_lag Number of periods to wait after treatment before forecasting begins (0 = in treatment year).
#' @param pretreatment_window "full" or "minimal".
#'
#' @return A list with:
#' \describe{
#'   \item{results}{Data frame with summary estimates across lags:
#'         \code{deg, hh, llag, FAT, sdFAT}.}
#'   \item{predictions}{Data frame with unit-level observed and predicted values, stacked across lags, with \code{llag}.}
#' }
#' @export
estimate_placebo_fat <- function(data,
                                 unit_var,
                                 time_var,
                                 outcome_var,
                                 treat_time_var,
                                 units_to_include = NULL,
                                 degrees = 0:2,
                                 horizons = 1:5,
                                 lags = 0:2,
                                 se_method = "analytic",
                                 n_bootstrap = 1000,
                                 covariate_vars = NULL,
                                 beta_estimator = c("none", "ols", "iv", "unitwise"),
                                 min_iv_lag = 2,
                                 max_iv_lag = 2,
                                 control_group_value = NULL,
                                 forecast_lag = 0,
                                 pretreatment_window = c("full", "minimal")) {

  # Match args to keep parity with estimate_fat()
  pretreatment_window <- match.arg(pretreatment_window)
  beta_estimator <- match.arg(beta_estimator)

  # Optional subsetting to user-specified units (same as estimate_fat)
  if (!is.null(units_to_include)) {
    data <- data[data[[unit_var]] %in% units_to_include, ]
  }

  # Storage
  res_list <- list()
  pred_list <- list()

  # Loop over placebo lags
  for (llag in lags) {
    # Build a pseudo-treatment column name to avoid mutating the original column
    pseudo_col <- paste0(treat_time_var, "_placebo_llag", llag)

    df_placebo <- data
    df_placebo[[pseudo_col]] <- df_placebo[[treat_time_var]] - llag

    # Run your existing estimator with the pseudo treatment column
    out_llag <- estimate_fat(
      data               = df_placebo,
      unit_var           = unit_var,
      time_var           = time_var,
      outcome_var        = outcome_var,
      treat_time_var     = pseudo_col,          # <— key handoff
      units_to_include   = NULL,                # already applied above if provided
      degrees            = degrees,
      horizons           = horizons,
      se_method          = se_method,
      n_bootstrap        = n_bootstrap,
      covariate_vars     = covariate_vars,
      beta_estimator     = beta_estimator,
      min_iv_lag         = min_iv_lag,
      max_iv_lag         = max_iv_lag,
      control_group_value= control_group_value,
      forecast_lag       = forecast_lag,
      pretreatment_window= pretreatment_window
    )

    # Tag with llag and collect
    res_list[[length(res_list) + 1]] <- transform(out_llag$results, llag = llag)
    pred_list[[length(pred_list) + 1]] <- transform(out_llag$predictions, llag = llag)
  }

  # Bind + enforce uniqueness on predictions just like estimate_fat does
  results_df <- dplyr::bind_rows(res_list) |>
    dplyr::select(deg, hh, llag, FAT, sdFAT) |>
    dplyr::arrange(deg, hh, llag)

  preds_df <- dplyr::bind_rows(pred_list) |>
    dplyr::distinct(!!rlang::sym(unit_var), !!rlang::sym(time_var), deg, hh, llag, .keep_all = TRUE)

  list(results = results_df, predictions = preds_df)
}

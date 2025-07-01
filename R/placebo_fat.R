#' Compute Placebo FAT estimates
#'
#' @param data Data.frame with required variables.
#' @param unit_var Unit variable name.
#' @param time_var Time variable name.
#' @param outcome_var Outcome variable name.
#' @param treat_year_var Adoption year variable.
#' @param degrees Vector of polynomial degrees.
#' @param lags Vector of placebo lags to apply.
#'
#' @return Data.frame of placebo FAT results.
#' @export

estimate_placebo_fat <- function(data, unit_var, time_var, outcome_var, treat_year_var,
                                 degrees = 1:2, lags = 0:2) {
  results <- list()

  for (d in degrees) {
    for (l in lags) {
      df <- data %>%
        group_by(.data[[unit_var]]) %>%
        mutate(adopt_year = .data[[treat_year_var]] - l,
               timeToTreat = .data[[time_var]] - adopt_year) %>%
        ungroup()

      df_preds <- df %>%
        group_by(.data[[unit_var]]) %>%
        group_modify(~ state_poly_fit(.x, time_var, outcome_var, unique(.x$adopt_year), d)) %>%
        ungroup()

      fat_val <- df_preds %>%
        filter(.data[[time_var]] == (adopt_year + 1)) %>%
        summarise(
          placeboFAT = mean(.data[[outcome_var]] - preds, na.rm = TRUE),
          sdFAT = sd(.data[[outcome_var]] - preds, na.rm = TRUE) / sqrt(n())
        ) %>%
        mutate(deg = d, llag = l)

      results[[length(results) + 1]] <- fat_val
    }
  }

  return(bind_rows(results))
}

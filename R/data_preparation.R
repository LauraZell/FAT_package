#' Prepare data for FAT estimation
#'
#' @param data A data.frame including outcome, unit, time, and treatment variables.
#' @param unit_var Name of the unit identifier variable (e.g., "state").
#' @param time_var Name of the time variable (e.g., "Year").
#' @param treat_var Name of the treatment variable (e.g., "Medical_Cannabis_Law").
#'
#' @return A data.frame with added variables: `adopt_year` and `timeToTreat`.
#' @export
prepare_fat_data <- function(data, unit_var, time_var, treat_var) {
  data <- data %>%
    group_by(.data[[unit_var]]) %>%
    mutate(adopt = ifelse(.data[[treat_var]] > 0 & lag(.data[[treat_var]]) == 0, 1, 0)) %>%
    mutate(adopt_year = ifelse(adopt == 1, .data[[time_var]], NA)) %>%
    fill(adopt_year, .direction = "downup") %>%
    mutate(timeToTreat = .data[[time_var]] - adopt_year)

  return(data)
}

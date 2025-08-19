# R/sim_data_helpers.R
# -------------------------------------------------------------------
# Two tiny helpers that reproduce your current simulation datasets
# without changing column names, formulas, or distributions.
# -------------------------------------------------------------------

#' Simulate FAT test data (no explicit control group)
#'
#' Reproduces the dataset created by dev/simulate_fat_testdata.R (first part).
#' Columns (unchanged): state, Year, adopt_year, timeToTreat, covariate1,
#' covariate2, instrument, ln_age_mort_rate.
#'
#' @param seed RNG seed. Default 123 to reproduce your script exactly.
#' @param n_units Number of units. Default 20.
#' @param n_years Number of years per unit. Default 10.
#' @param min_year First year. Default 2000.
#' @return A tibble with the simulated panel.
#' @export
simulate_fat_testdata <- function(seed = 123,
                                  n_units = 20,
                                  n_years = 10,
                                  min_year = 2000) {
  set.seed(seed)

  unit_names <- paste0("state", 1:n_units)

  data_list <- lapply(unit_names, function(unit) {
    treat_year <- sample(2004:2007, 1)  # Treatment years between 2004 and 2007
    years <- min_year:(min_year + n_years - 1)

    tibble::tibble(
      state = unit,
      Year = years,
      adopt_year = treat_year,
      timeToTreat = Year - treat_year,
      covariate1 = stats::rnorm(n_years, mean = 0, sd = 1),
      covariate2 = stats::rnorm(n_years, mean = 0, sd = 1),
      instrument = dplyr::lag(covariate1, 2),
      ln_age_mort_rate = 3 +
        0.05 * (years - min_year) -
        0.2 * (years >= treat_year) +
        0.3 * stats::rnorm(n_years) +
        0.5 * stats::rnorm(1) # base trend + noise + unit FE
    )
  })

  sim_data <- dplyr::bind_rows(data_list)

  # Clean up instrument column (remove NA due to lag) exactly as in your script
  sim_data <- sim_data %>%
    dplyr::group_by(state) %>%
    dplyr::mutate(instrument = dplyr::lag(covariate1, 2)) %>%
    dplyr::ungroup()

  sim_data
}


#' Simulate FAT test data with treated vs. control units
#'
#' Reproduces the dataset created by dev/simulate_fat_testdata.R (second part).
#' Columns (unchanged): state, Year, treated, adopt_year, time_to_treat,
#' covariate1, covariate2, y0, ln_age_mort_rate.
#'
#' @param seed RNG seed. Default 42 to reproduce your script exactly.
#' @param n_units Number of units. Default 40.
#' @param n_years Number of years per unit. Default 10.
#' @param start_year First year. Default 2000.
#' @return A data.frame with the simulated panel (column types match your code).
#' @export
simulate_dfat_testdata <- function(seed = 42,
                                   n_units = 40,
                                   n_years = 10,
                                   start_year = 2000) {
  set.seed(seed)

  units <- paste0("state", 1:n_units)
  years <- start_year + 0:(n_years - 1)

  sim_data_with_controls <- expand.grid(state = units, Year = years)
  sim_data_with_controls <- sim_data_with_controls[
    order(sim_data_with_controls$state, sim_data_with_controls$Year), ]

  # Half treated, half control
  treated_units <- sample(units, n_units / 2)
  sim_data_with_controls$treated <- sim_data_with_controls$state %in% treated_units

  # Assign adoption year only to treated units (keep your exact logic)
  sim_data_with_controls$adopt_year <- ifelse(
    sim_data_with_controls$treated,
    sample(2003:2006, length(treated_units), replace = TRUE)[
      match(sim_data_with_controls$state[sim_data_with_controls$treated], treated_units)
    ],
    NA
  )

  # Time to treatment (note: underscore name kept exactly)
  sim_data_with_controls$time_to_treat <-
    sim_data_with_controls$Year - sim_data_with_controls$adopt_year

  # Covariates
  sim_data_with_controls$covariate1 <- stats::rnorm(nrow(sim_data_with_controls))
  sim_data_with_controls$covariate2 <- stats::rnorm(nrow(sim_data_with_controls))

  # True untreated outcome (y0), unchanged
  sim_data_with_controls$y0 <- with(
    sim_data_with_controls,
    0.5 * covariate1 - 0.3 * covariate2 + 0.05 * (Year - start_year) +
      stats::rnorm(nrow(sim_data_with_controls), sd = 0.2)
  )

  # Add treatment effect after adoption for treated units, unchanged
  sim_data_with_controls$ln_age_mort_rate <- sim_data_with_controls$y0 +
    ifelse(sim_data_with_controls$treated &
             !is.na(sim_data_with_controls$adopt_year) &
             sim_data_with_controls$Year > sim_data_with_controls$adopt_year,
           0.25, 0)

  sim_data_with_controls
}

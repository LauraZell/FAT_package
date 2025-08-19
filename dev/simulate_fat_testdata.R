# dev/simulate_fat_testdata.R
# ------------------------------------------------------------
# This script generates simulated panel data for testing the
# estimate_fat() function under different specifications:
#   - Baseline (no covariates)
#   - Pooled OLS-adjusted
#   - IV-adjusted
#
# It is meant for development and internal testing purposes.
# ------------------------------------------------------------

set.seed(123)

# Parameters
n_units <- 20
n_years <- 10
min_year <- 2000
unit_names <- paste0("state", 1:n_units)

data_list <- lapply(unit_names, function(unit) {
  treat_year <- sample(2004:2007, 1)  # Treatment years between 2004 and 2007
  years <- min_year:(min_year + n_years - 1)

  tibble::tibble(
    state = unit,
    Year = years,
    adopt_year = treat_year,
    timeToTreat = Year - treat_year,
    covariate1 = rnorm(n_years, mean = 0, sd = 1),
    covariate2 = rnorm(n_years, mean = 0, sd = 1),
    instrument = lag(covariate1, 2),
    ln_age_mort_rate = 3 + 0.05 * (years - min_year) - 0.2 * (years >= treat_year) +
      0.3 * rnorm(n_years) + 0.5 * rnorm(1)  # base trend + noise + unit FE
  )
})

# Combine into panel
sim_data <- dplyr::bind_rows(data_list)

# Clean up instrument column (remove NA due to lag)
sim_data <- dplyr::group_by(sim_data, state) %>%
  dplyr::mutate(instrument = dplyr::lag(covariate1, 2)) %>%
  dplyr::ungroup()

# Save to global environment for immediate testing
assign("sim_data", sim_data, envir = .GlobalEnv)

# Optional: write to file for persistent testing or examples
# readr::write_csv(sim_data, "dev/simulated_fat_data.csv")


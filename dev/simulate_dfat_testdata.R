set.seed(42)

# Parameters
n_units <- 40
n_years <- 10
units <- paste0("state", 1:n_units)
years <- 2000 + 0:(n_years - 1)

# Simulate unit-level data
sim_data_with_controls <- expand.grid(state = units, Year = years)
sim_data_with_controls <- sim_data_with_controls[order(sim_data_with_controls$state,
                                                       sim_data_with_controls$Year), ]

# Half treated, half control
treated_units <- sample(units, n_units / 2)
sim_data_with_controls$treated <- sim_data_with_controls$state %in% treated_units

# Assign adoption year only to treated units
sim_data_with_controls$adopt_year <- ifelse(sim_data_with_controls$treated,
                                            sample(2003:2006, length(treated_units), replace = TRUE)[
                                              match(sim_data_with_controls$state[sim_data_with_controls$treated], treated_units)
                                            ],
                                            NA)

# Time to treatment
sim_data_with_controls$time_to_treat <- sim_data_with_controls$Year - sim_data_with_controls$adopt_year

# Simulate covariates and outcome
sim_data_with_controls$covariate1 <- rnorm(nrow(sim_data_with_controls))
sim_data_with_controls$covariate2 <- rnorm(nrow(sim_data_with_controls))

# Simulate true untreated outcome
sim_data_with_controls$y0 <- with(sim_data_with_controls,
  0.5 * covariate1 - 0.3 * covariate2 + 0.05 * (Year - 2000) +
    rnorm(nrow(sim_data_with_controls), sd = 0.2)
)

# Add treatment effect after adoption for treated units
sim_data_with_controls$ln_age_mort_rate <- sim_data_with_controls$y0 +
  ifelse(sim_data_with_controls$treated &
           !is.na(sim_data_with_controls$adopt_year) &
           sim_data_with_controls$Year > sim_data_with_controls$adopt_year,
         0.25, 0)

# Inspect
head(sim_data_with_controls)

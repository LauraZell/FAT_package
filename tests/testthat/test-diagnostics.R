test_that("diagnostics pass on baseline FAT output", {
  skip_on_cran()

  sim_data <- make_sim_data(n_units = 6, years = 2000:2008)

  res <- estimate_fat(
    data = sim_data,
    unit_var = "state",
    time_var = "Year",
    outcome_var = "ln_age_mort_rate",
    treat_time_var = "adopt_year",
    degrees = 0:2,
    horizons = 1:3,
    beta_estimator = "none",
    forecast_lag = 0,
    pretreatment_window = "full"
  )

  chk <- fat_validate(
    predictions = res$predictions,
    results     = res$results,
    data        = dplyr::mutate(sim_data, treat_time_for_fit = adopt_year),  # mirror your estimator
    unit_var    = "state",
    time_var    = "Year",
    outcome_var = "ln_age_mort_rate",
    forecast_lag = 0,
    dfat_mode    = FALSE
  )

  expect_true(all(chk$passed), info = paste(chk$check[!chk$passed], collapse = ", "))
})

test_that("diagnostics pass with 'minimal' window", {
  skip_on_cran()

  sim_data <- make_sim_data(n_units = 6, years = 2000:2008)

  res <- estimate_fat(
    data = sim_data,
    unit_var = "state",
    time_var = "Year",
    outcome_var = "ln_age_mort_rate",
    treat_time_var = "adopt_year",
    degrees = 1:2,           # deg=0 + minimal might be fragile
    horizons = 1:3,
    beta_estimator = "none",
    forecast_lag = 0,
    pretreatment_window = "minimal"
  )

  chk <- fat_validate(
    predictions = res$predictions,
    results     = res$results,
    data        = dplyr::mutate(sim_data, treat_time_for_fit = adopt_year),
    unit_var    = "state",
    time_var    = "Year",
    outcome_var = "ln_age_mort_rate",
    forecast_lag = 0,
    dfat_mode    = FALSE
  )

  expect_true(all(chk$passed), info = paste(chk$check[!chk$passed], collapse = ", "))
})

test_that("DFAT diagnostics (treated vs control) recompute", {
  skip_on_cran()

  sim_data <- make_sim_data(n_units = 8, years = 2000:2008)
  # fabricate treated flag (half treated)
  treated_units <- unique(sim_data$state)[1:4]
  sim_data$treated <- sim_data$state %in% treated_units

  res <- estimate_fat(
    data = sim_data,
    unit_var = "state",
    time_var = "Year",
    outcome_var = "ln_age_mort_rate",
    treat_time_var = "adopt_year",
    degrees = 1:2,
    horizons = 1:2,
    beta_estimator = "none",
    forecast_lag = 0,
    pretreatment_window = "full",
    control_group_value = FALSE
  )

  chk <- fat_validate(
    predictions = res$predictions,
    results     = res$results,
    data        = dplyr::mutate(sim_data, treat_time_for_fit = median(sim_data$adopt_year[sim_data$treated])),
    unit_var    = "state",
    time_var    = "Year",
    outcome_var = "ln_age_mort_rate",
    forecast_lag = 0,
    dfat_mode    = TRUE,
    control_group_value = FALSE
  )

  expect_true(all(chk$passed), info = paste(chk$check[!chk$passed], collapse = ", "))
})

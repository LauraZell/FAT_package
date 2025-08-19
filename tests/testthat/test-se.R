# tests/testthat/test-se.R
#testthat::context("SE methods for estimate_fat()")

data <- simulate_fat_testdata()

# Common quick args for speed in CI
deg_vec <- 0:1
hh_vec  <- 1:2

# ---- FAT + analytic SE ----------------------------------------------------
test_that("estimate_fat() returns valid results with analytic SE (FAT)", {
  data <- simulate_fat_testdata()

  out <- fatEstimation::estimate_fat(
    data = data,
    unit_var = "state",
    time_var = "Year",
    outcome_var = "ln_age_mort_rate",
    treat_time_var = "adopt_year",
    degrees = deg_vec,
    horizons = hh_vec,
    se_method = "analytic",
    beta_estimator = "none",
    pretreatment_window = "full"
  )

  expect_type(out, "list")
  expect_true(all(c("results", "predictions") %in% names(out)))

  res <- out$results
  expect_s3_class(res, "data.frame")
  expect_true(all(c("deg", "hh", "FAT", "sdFAT") %in% names(res)))
  expect_equal(nrow(res), length(deg_vec) * length(hh_vec))
  expect_true(all(is.finite(res$FAT)))
  expect_true(all(is.finite(res$sdFAT)))
  expect_true(all(res$sdFAT >= 0))
})

# ---- FAT + bootstrap SE ---------------------------------------------------
test_that("estimate_fat() returns valid results with bootstrap SE (FAT)", {
  set.seed(1)  # reproducible bootstraps
  data <- simulate_fat_testdata()

  out <- fatEstimation::estimate_fat(
    data = data,
    unit_var = "state",
    time_var = "Year",
    outcome_var = "ln_age_mort_rate",
    treat_time_var = "adopt_year",
    degrees = deg_vec,
    horizons = hh_vec,
    se_method = "bootstrap",
    n_bootstrap = 20,         # keep small for CI speed
    beta_estimator = "none",
    pretreatment_window = "full"
  )

  res <- out$results
  expect_s3_class(res, "data.frame")
  expect_true(all(c("deg", "hh", "FAT", "sdFAT") %in% names(res)))
  expect_equal(nrow(res), length(deg_vec) * length(hh_vec))
  expect_true(all(is.finite(res$FAT)))
  expect_true(all(is.finite(res$sdFAT)))
  expect_true(all(res$sdFAT >= 0))
})

# ---- FAT + clustered SE ---------------------------------------------------
test_that("estimate_fat() returns valid results with clustered SE (FAT)", {
  data <- simulate_fat_testdata()

  out <- fatEstimation::estimate_fat(
    data = data,
    unit_var = "state",
    time_var = "Year",
    outcome_var = "ln_age_mort_rate",
    treat_time_var = "adopt_year",
    degrees = deg_vec,
    horizons = hh_vec,
    se_method = "clustered",
    beta_estimator = "none",
    pretreatment_window = "full"
  )

  res <- out$results
  expect_s3_class(res, "data.frame")
  expect_true(all(c("deg", "hh", "FAT", "sdFAT") %in% names(res)))
  expect_equal(nrow(res), length(deg_vec) * length(hh_vec))
  expect_true(all(is.finite(res$FAT)))
  expect_true(all(is.finite(res$sdFAT)))
  expect_true(all(res$sdFAT >= 0))
})

# ---- DFAT (analytic SE) ---------------------------------------------------
test_that("estimate_fat() returns valid results with analytic SE (DFAT)", {
  data <- simulate_dfat_testdata()

  out <- fatEstimation::estimate_fat(
    data = data,
    unit_var = "state",
    time_var = "Year",
    outcome_var = "ln_age_mort_rate",
    treat_time_var = "adopt_year",
    degrees = deg_vec,
    horizons = hh_vec,
    se_method = "analytic",
    control_group_value = FALSE,   # control is treated == FALSE
    beta_estimator = "none",
    pretreatment_window = "full"
  )

  res <- out$results
  preds <- out$predictions

  expect_s3_class(res, "data.frame")
  expect_true(all(c("deg", "hh", "FAT", "sdFAT") %in% names(res)))
  expect_true(all(is.finite(res$FAT)))
  expect_true(all(is.finite(res$sdFAT)))

  # Predictions should carry treated and be unique per (unit,time,deg,hh)
  expect_true("treated" %in% names(preds))
  keys_n <- nrow(preds)
  keys_unique <- preds %>%
    distinct(.data$state, .data$Year, .data$deg, .data$hh) %>%
    nrow()
  expect_equal(keys_n, keys_unique)
})

# ---- Minimal vs Full window sanity check (FAT) ----------------------------
test_that("Minimal vs Full pretreatment windows both run and return expected shapes", {
  data <- simulate_fat_testdata()

  out_full <- fatEstimation::estimate_fat(
    data = data,
    unit_var = "state",
    time_var = "Year",
    outcome_var = "ln_age_mort_rate",
    treat_time_var = "adopt_year",
    degrees = deg_vec,
    horizons = hh_vec,
    se_method = "analytic",
    beta_estimator = "none",
    pretreatment_window = "full"
  )

  out_min <- fatEstimation::estimate_fat(
    data = data,
    unit_var = "state",
    time_var = "Year",
    outcome_var = "ln_age_mort_rate",
    treat_time_var = "adopt_year",
    degrees = deg_vec,
    horizons = hh_vec,
    se_method = "analytic",
    beta_estimator = "none",
    pretreatment_window = "minimal"
  )

  expect_equal(nrow(out_full$results), nrow(out_min$results))
  expect_setequal(names(out_full$results), names(out_min$results))
  expect_true(all(is.finite(out_full$results$FAT)))
  expect_true(all(is.finite(out_min$results$FAT)))
})

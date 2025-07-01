library(testthat)
library(fatEstimation)
library(dplyr)

# Create test dataset
test_data <- data.frame(
  state = rep(c("A", "B", "C"), each = 10),
  Year = rep(2001:2010, times = 3),
  ln_age_mort_rate = rnorm(30, mean = 2, sd = 0.5),
  Medical_Cannabis_Law = c(rep(0, 4), rep(1, 6), rep(0, 5), rep(1, 5), rep(0, 6), rep(1, 4))
) %>%
  group_by(state) %>%
  mutate(adopt = ifelse(Medical_Cannabis_Law > 0 & lag(Medical_Cannabis_Law) == 0, 1, 0),
         adopt_year = min(Year[adopt == 1], na.rm = TRUE)) %>%
  ungroup()

test_that("estimate_fat() returns valid results with analytic SE", {
  result <- estimate_fat(
    data = test_data,
    unit_var = "state",
    time_var = "Year",
    outcome_var = "ln_age_mort_rate",
    treat_time_var = "adopt_year",
    degrees = 0:1,
    horizons = 1:2,
    se_method = "analytic"
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("deg", "hh", "FAT", "sdFAT") %in% names(result)))
  expect_false(any(is.na(result$FAT)))
  expect_false(any(is.na(result$sdFAT)))
})

test_that("estimate_fat() returns valid results with bootstrap SE", {
  result <- estimate_fat(
    data = test_data,
    unit_var = "state",
    time_var = "Year",
    outcome_var = "ln_age_mort_rate",
    treat_time_var = "adopt_year",
    degrees = 0:1,
    horizons = 1:2,
    se_method = "bootstrap",
    n_bootstrap = 5
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("deg", "hh", "FAT", "sdFAT") %in% names(result)))
  expect_false(any(is.na(result$FAT)))
  expect_false(any(is.na(result$sdFAT)))
})

test_that("estimate_fat() returns valid results with clustered SE", {
  result <- estimate_fat(
    data = test_data,
    unit_var = "state",
    time_var = "Year",
    outcome_var = "ln_age_mort_rate",
    treat_time_var = "adopt_year",
    degrees = 0:1,
    horizons = 1:2,
    se_method = "clustered"
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("deg", "hh", "FAT", "sdFAT") %in% names(result)))
  expect_false(any(is.na(result$FAT)))
  expect_false(any(is.na(result$sdFAT)))
})

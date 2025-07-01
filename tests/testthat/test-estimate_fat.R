library(testthat)
library(dplyr)       # ← Add this
library(fatEstimation)


test_that("estimate_fat() runs and returns expected structure", {
  # Create test dataset
  test_data <- data.frame(
    state = rep(c("A", "B"), each = 10),
    Year = rep(2001:2010, times = 2),
    ln_age_mort_rate = rnorm(20, mean = 2, sd = 0.5),
    Medical_Cannabis_Law = c(rep(0, 4), rep(1, 6), rep(0, 6), rep(1, 4))
  ) %>%
    group_by(state) %>%
    mutate(adopt = ifelse(Medical_Cannabis_Law > 0 & lag(Medical_Cannabis_Law) == 0, 1, 0),
           adopt_year = min(Year[adopt == 1], na.rm = TRUE)) %>%
    ungroup()

  # Estimate FAT with bootstrapping (small n for test)
  fat_result <- estimate_fat(
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

  # Structure tests
  expect_s3_class(fat_result, "data.frame")
  expect_true(all(c("deg", "hh", "FAT", "sdFAT") %in% names(fat_result)))

  # Value tests
  expect_true(nrow(fat_result) > 0)
  expect_false(any(is.na(fat_result$FAT)))
})

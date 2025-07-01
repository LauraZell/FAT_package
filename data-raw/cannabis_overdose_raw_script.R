# data-raw/cannabis_overdose.R
# Script to prepare the example dataset for the fatEstimation package

library(dplyr)
library(readr)
library(usethis)

# Load the raw data (you must adjust this path to your local file)
raw_data <- read_csv("data-raw/pnas.csv")  # Replace with actual path if needed

# Prepare the data
cannabis_overdose <- raw_data %>%
  select(state, Year, Medical_Cannabis_Law, ln_age_mort_rate) %>%
  group_by(state) %>%
  mutate(adopt = ifelse(Medical_Cannabis_Law > 0 & lag(Medical_Cannabis_Law) == 0, 1, 0)) %>%
  filter(!is.na(Medical_Cannabis_Law)) %>%
  mutate(adopt_year = ifelse(adopt == 1, Year, NA)) %>%
  fill(adopt_year, .direction = "downup") %>%
  ungroup() %>%
  filter(!is.na(adopt_year)) %>%
  mutate(timeToTreat = Year - adopt_year)

# Save to data/
usethis::use_data(cannabis_overdose, overwrite = TRUE)

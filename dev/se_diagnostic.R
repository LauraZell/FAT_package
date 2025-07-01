# dev/se_diagnostic.R

# Load your package
devtools::load_all()

# Load required packages
library(dplyr)
library(ggplot2)
library(purrr)

# Load the real example data
data("cannabis_overdose")

# Run FAT estimation for all three SE methods
results_all_se <- purrr::map_dfr(
  c("analytic", "bootstrap", "clustered"),
  function(method) {
    estimate_fat(
      data = cannabis_overdose,
      unit_var = "state",
      time_var = "Year",
      outcome_var = "ln_age_mort_rate",
      treat_time_var = "adopt_year",
      se_method = method,
      degrees = 0:3,
      horizons = 1:5,
      n_bootstrap = 500
    ) %>%
      dplyr::mutate(se_method = method)
  }
)

# Plot results
ggplot(results_all_se, aes(x = hh, y = FAT, color = se_method)) +
  geom_line() +
  facet_wrap(~ deg, labeller = label_both) +
  geom_errorbar(aes(ymin = FAT - sdFAT, ymax = FAT + sdFAT), width = 0.2) +
  theme_minimal() +
  labs(
    title = "FAT Estimates by Standard Error Method",
    subtitle = "With 95% confidence intervals",
    x = "Post-Treatment Horizon (hh)",
    y = "FAT Estimate",
    color = "SE Method"
  )

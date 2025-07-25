
# fatEstimation

The `fatEstimation` R package implements the **Forecasted Average
Treatment (FAT)** methodology for estimating causal effects when a clean
control group is unavailable or unconvincing based on the research paper
by Botosaru, Giacomini, and Weidner 2023 (link?).

It allows researchers and policy analysts to estimate treatment effects
using **forecasted counterfactuals** based on pre-treatment trends.

For a full example and methodological overview, see the [package
vignette](https://your.package.website/articles/fat_vignette.html).

------------------------------------------------------------------------

## 🔍 Motivation

Traditional causal inference methods like DiD rely on finding comparable
control units. The FAT approach bypasses this need by:

- Fitting pre-treatment outcome trends for treated units individually
- Forecasting what would have happened without treatment
- Comparing the actual post-treatment outcomes to these forecasts

------------------------------------------------------------------------

## 📦 Installation

You can install the development version of the package directly from
GitHub using:

``` r
# install.packages("remotes")
remotes::install_github("laurazell/fatEstimation")

# From local source (if you're developing it)
devtools::load_all()
```

## Key features

- Forecasts counterfactual outcomes using unit-specific pre-trends.
- Supports multiple polynomial degrees and forecast horizons.
- Provides analytic, clustered, and bootstrap standard errors.
- Built-in diagnostic and placebo test functions.

## Example using the cannabis_overdose dataset

The package includes a cleaned dataset cannabis_overdose, derived from
the study by Shover et al. It tracks opioid overdose mortality rates
across U.S. states and the adoption of medical cannabis laws. (I assume
we might not have the right to use this, when making the package public,
but for now, it’s the data set for which we have FAT results for
comparison.)

``` r
library(fatEstimation)

data(cannabis_overdose)

# Run FAT estimation for deg = 0:2, hh = 1:3
results <- estimate_fat(
  data = cannabis_overdose,
  unit_var = "state",
  time_var = "Year",
  outcome_var = "ln_age_mort_rate",
  treat_time_var = "adopt_year",
  units_to_include = unique(cannabis_overdose$state),
  degrees = 0:2,
  horizons = 1:3,
  se_method = "analytic"
)

print(results)

# Plot results
library(ggplot2)
ggplot(results, aes(x = hh, y = FAT, color = as.factor(deg))) +
  geom_point() +
  geom_errorbar(aes(ymin = FAT - 1.96 * sdFAT, ymax = FAT + 1.96 * sdFAT)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Forecasted Average Treatment Effects",
       x = "Horizon (Years after adoption)",
       color = "Polynomial Degree")
       
# Diagnostic plot
plot_fat_diagnostics(results)
```

## License

MIT + file LICENSE

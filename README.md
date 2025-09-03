
# fatEstimator

The **`fatEstimator`** R package implements the *Forecasted Average
Treatment Effects (FAT)* methodology from [Botosaru, Giacomini, and
Weidner (2024)](https://arxiv.org/abs/2401.02121), which allows for
causal inference **without requiring a traditional control group**.

This method is especially useful in cases of:

- Universal or staggered treatment adoption,
- Endogenous policy timing,
- No clear or valid comparison units.

Instead of comparing treated units to a control group, the package
forecasts **unit-specific counterfactual outcomes** using pre-treatment
trends, and estimates the treatment effect as the deviation from these
forecasts.

For a full overview, see the
[**vignette**](./vignettes/fat_vignette.html) or run:

``` r
vignette("fat_vignette", package = "fatEstimator")
```

------------------------------------------------------------------------

## Motivation

Standard treatment effect estimators (e.g. Difference-in-Differences)
depend on finding valid comparison units. This can be difficult when:

- All units are eventually treated,
- Treatment timing is endogenous,
- Comparison groups differ structurally from treated units.

The **FAT estimator** bypasses this by: - Fitting a unit-specific model
to the *pre-treatment outcome history*, - Forecasting the counterfactual
outcome path, - Estimating treatment effects as the difference between
observed and forecasted outcomes.

It is formally grounded in recent econometric theory and allows for
flexible polynomial trends, covariates, placebo testing, and extensions
with control groups (DFAT).

------------------------------------------------------------------------

## Installation

Install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("laurazell/fatEstimator")

# From local source (if you're developing it)
devtools::load_all()
```

## Features

- Forecasted counterfactuals via unit-level polynomial trend regression
- Supports multiple forecast horizons and degrees
- Allows for model extensions with:
  - Homogeneous and heterogeneous covariates
  - Instrumental variables (IV) for endogeneity
  - Difference-in-FAT (DFAT) when a noisy control group is available
- Fast and flexible:
  - Computes analytic, clustered, or bootstrap standard errors
  - Diagnostics, placebo checks, and plotting functions built-in

## Citation

If you use this package, please cite:

Botosaru, I., Giacomini, R., & Weidner, M. (2024). Forecasted Average
Treatment Effects in the Absence of a Control Group. Working Paper.
arXiv:2401.02121

## License

MIT + file LICENSE

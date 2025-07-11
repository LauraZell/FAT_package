% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_fat.R
\name{estimate_fat}
\alias{estimate_fat}
\title{Estimate Forecasted Average Treatment Effects (FAT / DFAT)}
\usage{
estimate_fat(
  data,
  unit_var,
  time_var,
  outcome_var,
  treat_time_var,
  units_to_include = NULL,
  degrees = 0:2,
  horizons = 1:5,
  se_method = "analytic",
  n_bootstrap = 1000,
  covariate_vars = NULL,
  beta_estimator = c("none", "ols", "iv", "unitwise"),
  min_iv_lag = 2,
  max_iv_lag = 2,
  control_group_value = NULL
)
}
\arguments{
\item{data}{A long-format dataframe with panel structure.}

\item{unit_var}{Name of the column for the unit ID (e.g., "state", "mun_id").}

\item{time_var}{Name of the column for the time variable (e.g., "year").}

\item{outcome_var}{Name of the outcome variable to be forecasted. (e.g., "ln_age_mort_rate").}

\item{treat_time_var}{Name of the column indicating treatment year per unit. (e.g., "adopt_year").}

\item{units_to_include}{Optional vector of unit names to restrict analysis to (e.g., "state1", "state2").}

\item{degrees}{Vector of polynomial degrees for trend fitting (e.g., 0:2).}

\item{horizons}{Vector of forecast horizons to estimate (e.g., 1:5).}

\item{se_method}{Method for standard errors: "analytic", "bootstrap", "clustered", "unitwise".}

\item{n_bootstrap}{Number of bootstrap repetitions (only if \code{se_method} is "bootstrap").}

\item{covariate_vars}{Optional vector of column names to include as covariates.}

\item{beta_estimator}{Pooled beta estimation method: "none", "ols", "iv", or "unitwise".}

\item{min_iv_lag}{Minimum lag for IV estimation (if used).}

\item{max_iv_lag}{Maximum lag for IV estimation (if used).}

\item{control_group_value}{Optional. If set (e.g., control_group_value = FALSE), DFAT mode is activated.
Treated units are expected to be marked in a "treated" column (TRUE/FALSE).}
}
\value{
A list with:
\describe{
\item{results}{Data frame with summary estimates: deg, hh, FAT, sdFAT.}
\item{predictions}{Data frame with unit-level observed and predicted values across time, for plotting.}
}
}
\description{
Estimates forecasted average treatment effects (FAT) or difference-in-FAT
(DFAT) across a range of polynomial degrees and forecast horizons. Optionally
returns unit-level prediction trajectories for diagnostic and plotting use.
}

#' Plot Forecasted Average Treatment (FAT) Estimates
#'
#' Creates a general plot of FAT estimates with options to include error bars,
#' facet by polynomial degree, and color by standard error method.
#'
#' @param fat_df A data.frame returned by `estimate_fat()`, or a combined result with multiple `se_method` values.
#' @param include_ci Logical. Whether to include error bars based on `sdFAT` (default: TRUE).
#' @param facet_by_degree Logical. Whether to facet the plot by polynomial degree (default: TRUE).
#' @param color_by_se_method Logical. Whether to color lines by standard error method (default: TRUE).
#'
#' @return A ggplot object.
#' @export
plot_fat_old <- function(fat_df,
                     include_ci = TRUE,
                     facet_by_degree = TRUE,
                     color_by_se_method = TRUE) {
  p <- ggplot(fat_df, aes(x = hh, y = FAT))

  if (color_by_se_method && "se_method" %in% colnames(fat_df)) {
    p <- p + aes(color = se_method)
  }

  p <- p + geom_line()

  if (include_ci && "sdFAT" %in% colnames(fat_df)) {
    p <- p + geom_errorbar(aes(ymin = FAT - sdFAT, ymax = FAT + sdFAT),
                           width = 0.2, alpha = 0.6)
  }

  if (facet_by_degree && "deg" %in% colnames(fat_df)) {
    p <- p + facet_wrap(~ deg)
  }

  p <- p + theme_minimal() +
    labs(
      title = "Forecasted Average Treatment (FAT)",
      x = "Forecast Horizon (hh)",
      y = "FAT Estimate",
      color = if (color_by_se_method) "SE Method" else NULL
    )

  return(p)
}

#' Diagnostic Plot: FAT Estimates for One SE Method
#'
#' Focused visualization of FAT estimates with 95% confidence intervals for one SE method.
#'
#' @param fat_df A data.frame from `estimate_fat()` with column `se_method`.
#' @param se_filter Character: which standard error method to plot (e.g., "analytic").
#' @param x_var Variable on the x-axis (default: "hh").
#' @param color_var Variable used for grouping (default: "deg").
#'
#' @return A ggplot object.
#' @export
plot_fat_diagnostics <- function(fat_df,
                                 se_filter = "analytic",
                                 x_var = "hh",
                                 color_var = "deg") {
  if (!"se_method" %in% colnames(fat_df)) {
    stop("`se_method` column not found in `fat_df`. Please add it or set `se_method` in `estimate_fat()`.")
  }

  filtered <- dplyr::filter(fat_df, se_method == se_filter)

  ggplot(filtered, aes_string(x = x_var, y = "FAT", color = color_var)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(aes(ymin = FAT - 1.96 * sdFAT,
                      ymax = FAT + 1.96 * sdFAT),
                  width = 0.2) +
    theme_minimal() +
    labs(
      title = paste("FAT estimates with", se_filter, "standard errors"),
      x = x_var,
      y = "FAT",
      color = color_var
    )
}


#' Plot FAT Estimates with Error Bars
#'
#' @param fat_df A data.frame with FAT results. Must include columns: 'hh', 'FAT', 'sdFAT', 'deg'.
#' @param show_ci Logical. If TRUE, plots error bars using sdFAT. Default is TRUE.
#' @param facet_by_degree Logical. If TRUE, creates one plot per polynomial degree. Default is TRUE.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
#' @export
plot_fat <- function(fat_df, show_ci = TRUE, facet_by_degree = TRUE, title = NULL) {
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  # Basic plot setup
  p <- ggplot2::ggplot(fat_df, ggplot2::aes(x = hh, y = FAT)) +
    ggplot2::geom_line(ggplot2::aes(group = factor(deg), color = factor(deg))) +
    ggplot2::geom_point(ggplot2::aes(color = factor(deg)))

  # Optionally add standard error bars
  if (show_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = FAT - 1.96 * sdFAT,
        ymax = FAT + 1.96 * sdFAT,
        color = factor(deg)
      ),
      width = 0.2
    )
  }

  # Facet by degree if needed
  if (facet_by_degree) {
    p <- p + ggplot2::facet_wrap(~deg)
  }

  # Customize plot
  p <- p +
    ggplot2::labs(
      x = "Horizon (hh)",
      y = "FAT Estimate",
      color = "Degree",
      title = title
    ) +
    ggplot2::theme_minimal()

  return(p)
}



plot_fat_dfat_avg_trajectory <- function(predictions_df, mode = c("fat", "dfat"),
                                         unit_var = "state", time_var = "Year",
                                         outcome_var = "ln_age_mort_rate", pred_var = "preds") {

  mode <- match.arg(mode)
  library(dplyr)
  library(ggplot2)

  # Label treatment group
  if (mode == "dfat") {
    if (!"treated" %in% names(predictions_df)) stop("DFAT mode requires a 'treated' column.")
    predictions_df <- predictions_df %>%
      mutate(group = ifelse(treated == 1, "Treated", "Control"))
  } else {
    predictions_df$group <- "Treated"
  }

  # Aggregate: mean over units, by Year × group × deg
  agg_df <- predictions_df %>%
    group_by(Year, group, deg) %>%
    summarise(
      obs = mean(.data[[outcome_var]], na.rm = TRUE),
      pred = mean(.data[[pred_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(cols = c("obs", "pred"),
                        names_to = "type",
                        values_to = "value")

  # Plot: Average observed and predicted for each group
  ggplot(agg_df, aes(x = Year, y = value, color = type, linetype = type)) +
    geom_line(size = 1.1) +
    facet_wrap(~ group + deg, ncol = 2, labeller = label_both) +
    scale_color_manual(values = c("obs" = "black", "pred" = "blue"),
                       labels = c("Observed", "Forecasted")) +
    scale_linetype_manual(values = c("obs" = "solid", "pred" = "dashed")) +
    scale_x_continuous(breaks = scales::pretty_breaks(), labels = scales::number_format(accuracy = 1)) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom") +
    labs(
      title = paste0("Average Outcome Trajectories (", toupper(mode), " Mode)"),
      subtitle = "Solid = Observed, Dashed = Forecasted | Facets: Group × Degree",
      x = "Year",
      y = outcome_var,
      color = "Line Type",
      linetype = "Line Type"
    )
}


#' Plot FAT or DFAT trajectories for observed vs. forecasted values
#'
#' @param predictions_df A dataframe with observed and forecasted outcomes, e.g. `results_fat$predictions`
#' @param mode Character, either "fat" or "dfat"
#' @param unit_var Character, name of the unit identifier (e.g. "state")
#' @param time_var Character, name of the time variable (e.g. "Year")
#' @param outcome_var Character, observed outcome variable (e.g. "ln_age_mort_rate")
#' @param pred_var Character, forecasted value variable (e.g. "preds")
#'
#' @return A ggplot object showing average observed and forecasted trajectories
#' @export
plot_fat_dfat_trajectory <- function(predictions_df,
                                     mode = c("fat", "dfat"),
                                     unit_var = "state",
                                     time_var = "Year",
                                     outcome_var = "ln_age_mort_rate",
                                     pred_var = "preds") {

  mode <- match.arg(mode)

  library(ggplot2)
  library(dplyr)
  library(tidyr)

  # === Grouping logic ===
  if (mode == "dfat") {
    if (!"treated" %in% names(predictions_df)) stop("DFAT mode requires a 'treated' column.")
    predictions_df <- predictions_df %>%
      mutate(group = ifelse(treated == 1, "Treated", "Control"))
  } else {
    predictions_df$group <- "Treated"
  }

  # === Reshape into long format ===
  df_long <- predictions_df %>%
    select(deg, !!sym(time_var), !!sym(unit_var), group,
           obs = !!sym(outcome_var),
           pred = !!sym(pred_var)) %>%
    pivot_longer(cols = c("obs", "pred"), names_to = "type", values_to = "value") %>%
    filter(!is.na(value))  # Drop NAs (especially forecasts before treatment year)

  # === Average over groups per year ===
  df_avg <- df_long %>%
    group_by(deg, group, !!sym(time_var), type) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  # === Final plot ===
  ggplot(df_avg, aes(x = .data[[time_var]], y = value, color = type, linetype = type)) +
    geom_line(size = 1) +
    facet_wrap(~ deg, ncol = 1, labeller = label_both) +
    scale_color_manual(values = c("obs" = "black", "pred" = "blue"),
                       labels = c("Observed", "Forecasted")) +
    scale_linetype_manual(values = c("obs" = "solid", "pred" = "dashed")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom") +
    labs(
      title = paste0("Group-Averaged ", toupper(mode), " Trajectories"),
      subtitle = "Observed vs Forecasted Outcomes by Group and Polynomial Degree",
      x = "Year",
      y = "Outcome Value",
      color = "Type",
      linetype = "Type"
    )
}

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



plot_dfat_trajectory <- function(data, unit_var, time_var, outcome_var, treat_time_var,
                                 predictions, treated_value = TRUE, control_value = FALSE) {
  library(ggplot2)
  library(dplyr)

  # Get a treated unit and a control unit
  treated_unit <- data %>%
    filter(treated == treated_value) %>%
    pull(.data[[unit_var]]) %>%
    unique() %>%
    .[1]

  control_unit <- data %>%
    filter(treated == control_value) %>%
    pull(.data[[unit_var]]) %>%
    unique() %>%
    .[1]

  treat_year <- data %>%
    filter(.data[[unit_var]] == treated_unit) %>%
    pull(.data[[treat_time_var]]) %>%
    unique()

  # Merge predictions with data
  df <- data %>%
    filter(.data[[unit_var]] %in% c(treated_unit, control_unit)) %>%
    left_join(predictions, by = c(unit_var, time_var)) %>%
    mutate(group = ifelse(.data[[unit_var]] == treated_unit, "Treated", "Control"),
           observed = .data[[outcome_var]],
           forecast = preds,
           segment = ifelse(.data[[time_var]] < treat_year, "Observed", "Forecasted"))

  ggplot(df, aes(x = .data[[time_var]], group = interaction(group, segment))) +
    geom_line(aes(y = observed, color = group, linetype = "Observed"), linewidth = 1.1) +
    geom_line(aes(y = forecast, color = group, linetype = "Forecasted"), linewidth = 1.1) +
    geom_vline(xintercept = treat_year, linetype = "dotted") +
    scale_linetype_manual(values = c("Observed" = "solid", "Forecasted" = "dashed")) +
    labs(title = "DFAT trajectories: Treated vs. Control",
         y = outcome_var, x = time_var, color = "Group", linetype = "Data type") +
    theme_minimal()
}




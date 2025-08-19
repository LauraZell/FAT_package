#' Validate internal consistency of FAT/DFAT outputs
#'
#' Runs a suite of invariants on the outputs of `estimate_fat()` and returns
#' a tibble of checks with pass/fail and short details.
#'
#' @param predictions The `results$predictions` data.frame from `estimate_fat()`.
#' @param results      The `results$results` data.frame from `estimate_fat()`.
#' @param data         The original long panel used as input to `estimate_fat()`.
#' @param unit_var     Unit id column (string).
#' @param time_var     Time column (string).
#' @param outcome_var  Outcome column (string).
#' @param forecast_lag Integer; same lag you passed to `estimate_fat()`.
#' @param dfat_mode    Logical; TRUE if you ran DFAT (i.e., you had a `treated` column and set `control_group_value`).
#' @param control_group_value Value used to mark the control group in `treated` (e.g., FALSE or 0). Required if `dfat_mode=TRUE`.
#'
#' @return tibble with columns: check, passed (logical), details (character)
#' @export
fat_validate <- function(predictions,
                         results,
                         data,
                         unit_var,
                         time_var,
                         outcome_var,
                         forecast_lag = 0,
                         dfat_mode = FALSE,
                         control_group_value = NULL) {
  stopifnot(is.data.frame(predictions), is.data.frame(results), is.data.frame(data))
  if (dfat_mode && is.null(control_group_value)) {
    stop("In DFAT mode, please supply `control_group_value`.")
  }

  # symbols for tidy-eval
  key_unit <- rlang::sym(unit_var)
  key_time <- rlang::sym(time_var)
  ysym     <- rlang::sym(outcome_var)

  # ---- 1) Uniqueness: at most one row per (unit,time,deg,hh)
  dup_df <- predictions %>%
    dplyr::count(!!key_unit, !!key_time, .data$deg, .data$hh, name = "n") %>%
    dplyr::filter(.data$n > 1)

  chk1 <- tibble::tibble(
    check   = "unique_keys",
    passed  = nrow(dup_df) == 0,
    details = if (nrow(dup_df) == 0) "OK"
    else paste0("Found ", nrow(dup_df), " duplicate (unit,time,deg,hh) keys; e.g. ",
                utils::capture.output(print(utils::head(dup_df, 5))) |> paste(collapse = " "))
  )

  # ---- 2) hh labeling rule: outside forecast window => hh==0; inside => hh==(timeToTreat-forecast_lag+1)
  # Only evaluate where timeToTreat is present
  if (!("timeToTreat" %in% names(predictions))) {
    chk2 <- tibble::tibble(
      check = "hh_labeling",
      passed = FALSE,
      details = "Column `timeToTreat` missing from predictions."
    )
  } else {
    lab <- predictions %>%
      dplyr::mutate(hh_expected = dplyr::if_else(
        .data$timeToTreat >= forecast_lag,
        .data$timeToTreat - forecast_lag + 1L,
        0L
      )) %>%
      dplyr::mutate(mismatch = .data$hh != .data$hh_expected)

    n_mismatch <- sum(lab$mismatch %in% TRUE, na.rm = TRUE)

    chk2 <- tibble::tibble(
      check   = "hh_labeling",
      passed  = n_mismatch == 0,
      details = if (n_mismatch == 0) "OK"
      else paste0("Found ", n_mismatch, " rows with hh != expected.")
    )
  }

  # ---- 3) Horizon bounds for predictions:
  # we only check rows where preds is non-NA; they must satisfy
  # forecast_lag <= timeToTreat <= forecast_lag + (hh - 1) for the *current* hh.
  if (!("preds" %in% names(predictions))) {
    chk3 <- tibble::tibble(
      check = "horizon_bounds",
      passed = FALSE,
      details = "Column `preds` missing from predictions."
    )
  } else {
    # rows with a filled prediction
    pr_non_na <- predictions %>%
      dplyr::filter(!is.na(.data$preds))

    # For those rows, hh should equal timeToTreat-forecast_lag+1 and be >=1
    bad_bounds <- pr_non_na %>%
      dplyr::mutate(hh_expected = .data$timeToTreat - forecast_lag + 1L) %>%
      dplyr::filter(.data$hh != .data$hh_expected | .data$hh < 1L)

    chk3 <- tibble::tibble(
      check   = "horizon_bounds",
      passed  = nrow(bad_bounds) == 0,
      details = if (nrow(bad_bounds) == 0) "OK"
      else paste0("Rows with preds but out-of-bounds hh/timeToTreat: ", nrow(bad_bounds))
    )
  }

  # ---- 4) Pretreatment presence in predictions: all pre-treatment rows should be present with preds==NA
  # (This checks you kept all observed outcomes pre-treatment.)
  pre_in_data <- data %>%
    dplyr::mutate(timeToTreat = !!key_time - data$treat_time_for_fit) %>%
    dplyr::filter(.data$timeToTreat < 0) %>%
    dplyr::select(!!key_unit, !!key_time) %>%
    dplyr::distinct()

  pre_in_pred <- predictions %>%
    dplyr::filter(.data$timeToTreat < 0) %>%
    dplyr::select(!!key_unit, !!key_time, .data$preds)

  pre_join <- pre_in_data %>%
    dplyr::left_join(pre_in_pred, by = c(unit_var, time_var))

  missing_pre_rows <- sum(is.na(pre_join$preds)) # NA preds are expected; we want to see if any pre rows are missing entirely
  # If left_join dropped nothing, nrow(pre_join) == nrow(pre_in_data)
  pre_rows_lost <- nrow(pre_in_data) - nrow(pre_join)

  chk4 <- tibble::tibble(
    check   = "pretreatment_rows_present",
    passed  = pre_rows_lost == 0,
    details = if (pre_rows_lost == 0) "OK" else paste0("Lost ", pre_rows_lost, " pre-treatment rows in predictions.")
  )

  # ---- 5) Recompute FAT (or DFAT) from predictions and compare to results
  # For each (deg, hh), target rows are those where predictions correspond to that step:
  # i.e. hh (column) equals that horizon. diff = y - preds.
  if (!all(c("deg", "hh") %in% names(predictions)) || !("FAT" %in% names(results))) {
    chk5 <- tibble::tibble(
      check = "fat_recompute",
      passed = FALSE,
      details = "Missing columns `deg`, `hh` in predictions or `FAT` in results."
    )
  } else {
    # build target rows for each (deg,hh)
    targ <- predictions %>%
      dplyr::filter(.data$hh >= 1L) %>%          # only forecasted steps
      dplyr::mutate(diff = !!ysym - .data$preds)

    if (!dfat_mode) {
      recomputed <- targ %>%
        dplyr::group_by(.data$deg, .data$hh) %>%
        dplyr::summarise(FAT_hat = mean(.data$diff, na.rm = TRUE), .groups = "drop")
    } else {
      # need treated flag from data
      treated_df <- data %>%
        dplyr::select(!!key_unit, !!key_time, .data$treated)

      targ2 <- targ %>%
        dplyr::left_join(treated_df, by = c(unit_var, time_var))

      recomputed <- targ2 %>%
        dplyr::group_by(.data$deg, .data$hh) %>%
        dplyr::summarise(
          FAT_hat = mean(.data$diff[.data$treated != control_group_value], na.rm = TRUE) -
            mean(.data$diff[.data$treated == control_group_value], na.rm = TRUE),
          .groups = "drop"
        )
    }

    comp <- results %>%
      dplyr::select(.data$deg, .data$hh, .data$FAT) %>%
      dplyr::left_join(recomputed, by = c("deg", "hh")) %>%
      dplyr::mutate(err = abs(.data$FAT - .data$FAT_hat))

    tol <- 1e-8
    n_bad <- sum(is.na(comp$FAT_hat) | comp$err > tol)

    chk5 <- tibble::tibble(
      check   = ifelse(dfat_mode, "dfat_recompute", "fat_recompute"),
      passed  = n_bad == 0,
      details = if (n_bad == 0) "OK"
      else paste0("Found ", n_bad, " (deg,hh) with mismatch > ", tol,
                  ". Examples: ",
                  utils::capture.output(print(utils::head(comp[is.na(comp$FAT_hat) | comp$err > tol, ], 5))) |>
                    paste(collapse = " "))
    )
  }

  dplyr::bind_rows(chk1, chk2, chk3, chk4, chk5)
}

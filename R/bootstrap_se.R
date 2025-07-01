#' Compute bootstrap standard error for FAT
#'
#' @param df Data.frame with observed and predicted outcomes.
#' @param outcome_var Outcome variable name.
#' @param pred_var Predicted outcome variable name.
#' @param B Number of bootstrap replications.
#' @param seed Optional seed for reproducibility.
#'
#' @return Bootstrap standard error (numeric).
#' @export
bootstrap_fat_se <- function(df, outcome_var, pred_var, B = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(df)
  fat_vals <- replicate(B, {
    sample_idx <- sample(1:n, size = n, replace = TRUE)
    mean(df[[outcome_var]][sample_idx] - df[[pred_var]][sample_idx])
  })

  return(sd(fat_vals))
}

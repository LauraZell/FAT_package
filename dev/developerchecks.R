# dev/developerchecks.R
# ------------------------------------------------------------
# Developer-friendly wrapper to run and summarize package tests.
# Provides shortcuts for:
#   - All tests
#   - Subsets (SE tests, diagnostics tests, etc.)
#
# Usage: source("dev/developerchecks.R")
# ------------------------------------------------------------

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Please install the 'devtools' package to run developer checks.")
}

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Please install the 'testthat' package to run developer checks.")
}

cat("=== Developer Checks for fatEstimation ===\n")

# Run all tests
cat("\n--- Running ALL tests ---\n")
res_all <- devtools::test()
print(res_all)

# Run only SE tests
cat("\n--- Running SE tests only ---\n")
res_se <- devtools::test(filter = "se")
print(res_se)

# Run only diagnostics tests
cat("\n--- Running diagnostics tests only ---\n")
res_diag <- devtools::test(filter = "diagnostics")
print(res_diag)

# Run only plotting tests (if you add them later)
if (file.exists("tests/testthat/test-plot.R")) {
  cat("\n--- Running plotting tests only ---\n")
  res_plot <- devtools::test(filter = "plot")
  print(res_plot)
}

cat("\n=== Developer Checks Complete ===\n")

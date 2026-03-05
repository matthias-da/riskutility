# Tests for drisk (Disclosure Risk for Continuous Variables)

# --- Setup: shared test data ---

make_numeric_data <- function(n = 100, p = 3, seed = 123) {
  set.seed(seed)
  data.frame(
    age = rnorm(n, 40, 10),
    income = rnorm(n, 50000, 15000),
    hours = rnorm(n, 40, 8)
  )[, seq_len(p), drop = FALSE]
}

make_independent_synth <- function(n = 100, p = 3, seed = 456) {
  set.seed(seed)
  data.frame(
    age = rnorm(n, 40, 10),
    income = rnorm(n, 50000, 15000),
    hours = rnorm(n, 40, 8)
  )[, seq_len(p), drop = FALSE]
}

# --- Class structure and fields ---

test_that("drisk returns correct S3 class structure", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "both")

  expect_s3_class(result, "drisk")
  expect_true("drisk_interval" %in% names(result))
  expect_true("drisk_rmd" %in% names(result))
  expect_true("at_risk_interval" %in% names(result))
  expect_true("at_risk_rmd" %in% names(result))
  expect_true("min_rmd" %in% names(result))
  expect_true("interval_widths" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("privacy_pass" %in% names(result))
  expect_true("n_original" %in% names(result))
  expect_true("n_synthetic" %in% names(result))
  expect_true("vars" %in% names(result))
  expect_true("outlier_par" %in% names(result))
  expect_true("alpha" %in% names(result))
})

test_that("drisk privacy_pass is logical", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y)

  expect_type(result$privacy_pass, "logical")
  expect_length(result$privacy_pass, 1)
})

test_that("drisk stores correct dataset sizes", {
  X <- make_numeric_data(n = 80)
  Y <- make_independent_synth(n = 60)

  result <- drisk(X, Y)

  expect_equal(result$n_original, 80)
  expect_equal(result$n_synthetic, 60)
})

test_that("drisk stores correct variables", {
  X <- make_numeric_data(p = 3)
  Y <- make_independent_synth(p = 3)

  result <- drisk(X, Y)

  expect_equal(sort(result$vars), c("age", "hours", "income"))
})

# --- Method selection ---

test_that("drisk method='interval' only computes interval", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "interval")

  expect_false(is.na(result$drisk_interval))
  expect_true(is.na(result$drisk_rmd))
  expect_true(!is.null(result$at_risk_interval))
  expect_null(result$at_risk_rmd)
  expect_null(result$min_rmd)
  expect_true(!is.null(result$interval_widths))
})

test_that("drisk method='rmd' only computes RMD", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "rmd")

  expect_true(is.na(result$drisk_interval))
  expect_false(is.na(result$drisk_rmd))
  expect_null(result$at_risk_interval)
  expect_true(!is.null(result$at_risk_rmd))
  expect_true(!is.null(result$min_rmd))
  expect_null(result$interval_widths)
})

test_that("drisk method='both' computes both methods", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "both")

  expect_false(is.na(result$drisk_interval))
  expect_false(is.na(result$drisk_rmd))
  expect_true(!is.null(result$at_risk_interval))
  expect_true(!is.null(result$at_risk_rmd))
  expect_true(!is.null(result$min_rmd))
  expect_true(!is.null(result$interval_widths))
})

# --- Interval method (dRisk) ---

test_that("drisk_interval is between 0 and 1", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "interval")

  expect_true(result$drisk_interval >= 0 && result$drisk_interval <= 1)
})

test_that("at_risk_interval has correct length", {
  X <- make_numeric_data(n = 80)
  Y <- make_independent_synth(n = 60)

  result <- drisk(X, Y, method = "interval")

  expect_length(result$at_risk_interval, 80)
  expect_type(result$at_risk_interval, "logical")
})

test_that("drisk_interval equals fraction of at-risk records", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "interval")

  expect_equal(result$drisk_interval, mean(result$at_risk_interval))
})

test_that("interval_widths are named and positive", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "interval")

  expect_true(all(result$interval_widths >= 0))
  expect_equal(names(result$interval_widths), result$vars)
})

test_that("interval method detects exact copies", {
  set.seed(123)
  X <- make_numeric_data(n = 50)
  # Y is near-exact copy with tiny noise
  Y <- X + matrix(rnorm(50 * 3, sd = 0.001), ncol = 3)
  names(Y) <- names(X)

  result <- drisk(X, Y, method = "interval")

  # Most records should be at risk when Y is nearly identical to X
  expect_true(result$drisk_interval > 0.5)
})

test_that("interval method shows low risk for independent data", {
  set.seed(42)
  X <- data.frame(
    x1 = rnorm(200),
    x2 = rnorm(200),
    x3 = rnorm(200)
  )
  Y <- data.frame(
    x1 = rnorm(200),
    x2 = rnorm(200),
    x3 = rnorm(200)
  )

  result <- drisk(X, Y, method = "interval")

  # With default outlier_par=0.01, independent data should have low risk
  expect_true(result$drisk_interval < 0.5)
})

test_that("larger outlier_par increases interval risk", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result_small <- drisk(X, Y, method = "interval", outlier_par = 0.001)
  result_large <- drisk(X, Y, method = "interval", outlier_par = 0.1)

  expect_true(result_large$drisk_interval >= result_small$drisk_interval)
})

test_that("interval widths scale with outlier_par", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  r1 <- drisk(X, Y, method = "interval", outlier_par = 0.01)
  r2 <- drisk(X, Y, method = "interval", outlier_par = 0.02)

  # Widths should be exactly 2x with 2x outlier_par
  ratio <- r2$interval_widths / r1$interval_widths
  expect_equal(unname(ratio),
               rep(2, length(r1$interval_widths)),
               tolerance = 1e-10)
})

# --- RMD method (dRiskRMD) ---

test_that("drisk_rmd is between 0 and 1", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "rmd")

  expect_true(result$drisk_rmd >= 0 && result$drisk_rmd <= 1)
})

test_that("at_risk_rmd has correct length", {
  X <- make_numeric_data(n = 80)
  Y <- make_independent_synth(n = 60)

  result <- drisk(X, Y, method = "rmd")

  expect_length(result$at_risk_rmd, 80)
  expect_type(result$at_risk_rmd, "logical")
})

test_that("drisk_rmd equals fraction of at-risk records", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "rmd")

  expect_equal(result$drisk_rmd, mean(result$at_risk_rmd))
})

test_that("min_rmd has correct length and is non-negative", {
  X <- make_numeric_data(n = 80)
  Y <- make_independent_synth(n = 60)

  result <- drisk(X, Y, method = "rmd")

  expect_length(result$min_rmd, 80)
  expect_true(all(result$min_rmd >= 0))
})

test_that("RMD method detects near-copies", {
  set.seed(123)
  X <- make_numeric_data(n = 50)
  # Y is near-exact copy with tiny noise
  Y <- X + matrix(rnorm(50 * 3, sd = 0.01), ncol = 3)
  names(Y) <- names(X)

  result <- drisk(X, Y, method = "rmd")

  # Near-copies should have high risk
  expect_true(result$drisk_rmd > 0.5)
})

test_that("smaller alpha means fewer flagged records", {
  X <- make_numeric_data(n = 100)
  Y <- make_independent_synth(n = 100)

  r_lenient <- drisk(X, Y, method = "rmd", alpha = 0.10)
  r_strict <- drisk(X, Y, method = "rmd", alpha = 0.01)

  expect_true(r_strict$drisk_rmd <= r_lenient$drisk_rmd)
})

test_that("RMD method requires at least 2 variables", {
  X <- data.frame(x = rnorm(50))
  Y <- data.frame(x = rnorm(50))

  expect_error(drisk(X, Y, method = "rmd"), "at least 2")
})

# --- Variable selection ---

test_that("drisk auto-selects common numeric variables", {
  set.seed(123)
  X <- data.frame(
    x1 = rnorm(50), x2 = rnorm(50),
    cat = sample(c("a", "b"), 50, replace = TRUE)
  )
  Y <- data.frame(
    x1 = rnorm(50), x2 = rnorm(50),
    cat = sample(c("a", "b"), 50, replace = TRUE)
  )

  result <- drisk(X, Y, method = "interval")

  # Should auto-select only numeric variables
  expect_equal(sort(result$vars), c("x1", "x2"))
})

test_that("drisk uses specified vars", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, vars = c("age", "income"), method = "interval")

  expect_equal(result$vars, c("age", "income"))
})

test_that("drisk errors on non-numeric specified vars", {
  X <- data.frame(x = rnorm(50), cat = letters[1:50])
  Y <- data.frame(x = rnorm(50), cat = letters[1:50])

  expect_error(drisk(X, Y, vars = c("x", "cat")), "numeric")
})

test_that("drisk errors when no common numeric vars found", {
  X <- data.frame(a = c("x", "y", "z"))
  Y <- data.frame(a = c("x", "y", "z"))

  expect_error(drisk(X, Y), "No common numeric variables")
})

# --- Input validation ---

test_that("drisk input validation works", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  expect_error(drisk(1:10, Y), "X must be a data frame")
  expect_error(drisk(X, 1:10), "Y must be a data frame")
})

test_that("drisk errors with invalid outlier_par", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  expect_error(drisk(X, Y, outlier_par = -1), "positive number")
  expect_error(drisk(X, Y, outlier_par = 0), "positive number")
  expect_error(drisk(X, Y, outlier_par = "abc"), "positive number")
})

test_that("drisk errors with invalid alpha", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  expect_error(drisk(X, Y, method = "rmd", alpha = 0), "between 0 and 1")
  expect_error(drisk(X, Y, method = "rmd", alpha = 1), "between 0 and 1")
  expect_error(drisk(X, Y, method = "rmd", alpha = -0.5), "between 0 and 1")
})

test_that("drisk errors when no common variables", {
  X <- data.frame(a = rnorm(10))
  Y <- data.frame(b = rnorm(10))

  expect_error(drisk(X, Y), "No common variables")
})

test_that("drisk errors with missing variables", {
  X <- data.frame(a = rnorm(10))
  Y <- data.frame(a = rnorm(10), b = rnorm(10))

  expect_error(drisk(X, Y, vars = c("a", "b")), "Variables missing in X")
})

test_that("drisk handles NA values", {
  set.seed(123)
  X <- data.frame(x1 = c(rnorm(48), NA, NA), x2 = rnorm(50))
  Y <- data.frame(x1 = c(rnorm(29), NA), x2 = rnorm(30))

  result <- drisk(X, Y, method = "interval", na.rm = TRUE)

  expect_s3_class(result, "drisk")
  expect_equal(result$n_original, 48)
  expect_equal(result$n_synthetic, 29)
})

test_that("drisk errors with no complete cases", {
  X <- data.frame(x1 = rep(NA_real_, 5), x2 = rep(NA_real_, 5))
  Y <- data.frame(x1 = as.numeric(1:5), x2 = as.numeric(1:5))

  expect_error(drisk(X, Y, na.rm = TRUE), "No complete cases")
})

# --- synth_pair method ---

test_that("drisk.synth_pair works", {
  set.seed(123)
  original <- data.frame(
    x1 = rnorm(60),
    x2 = rnorm(60),
    cat = sample(c("a", "b"), 60, replace = TRUE)
  )
  synthetic <- data.frame(
    x1 = rnorm(40),
    x2 = rnorm(40),
    cat = sample(c("a", "b"), 40, replace = TRUE)
  )

  pair <- synth_pair(original = original, synthetic = synthetic)
  result <- drisk(pair, method = "interval")

  expect_s3_class(result, "drisk")
  # Should only use numeric vars from synth_pair
  expect_equal(sort(result$vars), c("x1", "x2"))
})

# --- Privacy pass ---

test_that("privacy_pass is TRUE for low risk", {
  set.seed(42)
  X <- data.frame(x1 = rnorm(200), x2 = rnorm(200), x3 = rnorm(200))
  Y <- data.frame(x1 = rnorm(200), x2 = rnorm(200), x3 = rnorm(200))

  result <- drisk(X, Y, method = "interval")

  # Independent data with small outlier_par should pass
  if (result$drisk_interval <= 0.1) {
    expect_true(result$privacy_pass)
  }
})

test_that("privacy_pass is FALSE for high risk", {
  set.seed(123)
  X <- make_numeric_data(n = 50)
  Y <- X + matrix(rnorm(50 * 3, sd = 0.001), ncol = 3)
  names(Y) <- names(X)

  result <- drisk(X, Y, method = "interval")

  # Near-copies should fail
  expect_false(result$privacy_pass)
})

# --- print, summary, plot ---

test_that("drisk print method works", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y)

  expect_output(print(result), "Disclosure Risk")
  expect_output(print(result), "Interval Method")
  expect_output(print(result), "RMD Method")
  expect_output(print(result), "Privacy Assessment")
})

test_that("drisk print method works for interval only", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "interval")

  expect_output(print(result), "Interval Method")
})

test_that("drisk print method works for rmd only", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "rmd")

  expect_output(print(result), "RMD Method")
})

test_that("drisk summary returns correct class", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y)
  s <- summary(result)

  expect_s3_class(s, "summary.drisk")
  expect_true("drisk_interval" %in% names(s))
  expect_true("drisk_rmd" %in% names(s))
  expect_true("privacy_pass" %in% names(s))
  expect_true("vars" %in% names(s))
})

test_that("drisk summary with RMD includes quantiles", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "rmd")
  s <- summary(result)

  expect_true("min_rmd_quantiles" %in% names(s))
  expect_true("mean_min_rmd" %in% names(s))
  expect_true("sd_min_rmd" %in% names(s))
  expect_true("threshold_rmd" %in% names(s))
})

test_that("drisk print.summary works", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y)
  s <- summary(result)

  expect_output(print(s), "Summary: Disclosure Risk")
  expect_output(print(s), "Method:")
  expect_output(print(s), "Privacy:")
})

test_that("drisk plot which=1 works for both method", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y)

  expect_no_error(plot(result, which = 1))
})

test_that("drisk plot which=1 works for interval only", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "interval")

  expect_no_error(plot(result, which = 1))
})

test_that("drisk plot which=2 works for rmd method", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "rmd")

  expect_no_error(plot(result, which = 2))
})

test_that("drisk plot which=1:2 works", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y)

  expect_no_error(plot(result, which = 1:2))
})

test_that("drisk plot which=2 shows message when no RMD data", {
  X <- make_numeric_data()
  Y <- make_independent_synth()

  result <- drisk(X, Y, method = "interval")

  expect_message(plot(result, which = 2), "requires method")
})

# --- Two-variable case ---

test_that("drisk works with exactly 2 variables", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  Y <- data.frame(x1 = rnorm(50), x2 = rnorm(50))

  result <- drisk(X, Y, method = "both")

  expect_s3_class(result, "drisk")
  expect_false(is.na(result$drisk_interval))
  expect_false(is.na(result$drisk_rmd))
})

# --- Edge case: constant variable ---

test_that("interval method handles constant variable (IQR=0)", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(50), x2 = rep(5, 50))
  Y <- data.frame(x1 = rnorm(50), x2 = rep(5, 50))

  result <- drisk(X, Y, method = "interval")

  # Constant variable has IQR=0, so interval width = 0
  # Only exact matches on that variable count
  expect_s3_class(result, "drisk")
  expect_equal(result$interval_widths["x2"], c(x2 = 0))
})

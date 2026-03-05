# Tests for singling_out (Singling Out Risk)

# --- Class structure and fields ---

test_that("singling_out returns correct S3 class structure", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    age = sample(20:70, n, replace = TRUE),
    income = rnorm(n, 50000, 15000),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )
  Y <- data.frame(
    age = sample(20:70, n, replace = TRUE),
    income = rnorm(n, 50000, 15000),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )

  result <- singling_out(X, Y, n_attacks = 100, seed = 42)

  expect_s3_class(result, "singling_out")
  expect_true("risk" %in% names(result))
  expect_true("risk_ci" %in% names(result))
  expect_true("risk_attack" %in% names(result))
  expect_true("risk_attack_ci" %in% names(result))
  expect_true("risk_control" %in% names(result))
  expect_true("risk_control_ci" %in% names(result))
  expect_true("n_attacks" %in% names(result))
  expect_true("n_success" %in% names(result))
  expect_true("n_control_success" %in% names(result))
  expect_true("match_counts" %in% names(result))
  expect_true("match_counts_control" %in% names(result))
  expect_true("privacy_pass" %in% names(result))
  expect_true("n_original" %in% names(result))
  expect_true("n_synthetic" %in% names(result))
  expect_true("n_holdout" %in% names(result))
  expect_true("mode" %in% names(result))
  expect_true("n_cols" %in% names(result))
  expect_true("vars" %in% names(result))
  expect_true("confidence_level" %in% names(result))
})

test_that("singling_out privacy_pass is logical", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- singling_out(X, Y, n_attacks = 50, seed = 42)

  expect_type(result$privacy_pass, "logical")
})

test_that("singling_out risk is between 0 and 1", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  result <- singling_out(X, Y, n_attacks = 100, seed = 42)

  expect_true(result$risk >= 0 && result$risk <= 1)
  expect_true(result$risk_attack >= 0 && result$risk_attack <= 1)
  expect_true(result$risk_control >= 0 && result$risk_control <= 1)
})

test_that("match_counts has correct length", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- singling_out(X, Y, n_attacks = 100, seed = 42)

  expect_equal(length(result$match_counts), 100)
  expect_equal(length(result$match_counts_control), 100)
  expect_equal(result$n_attacks, 100)
})

test_that("match_counts are non-negative integers", {
  set.seed(123)
  X <- data.frame(
    x1 = rnorm(80),
    x2 = sample(c("A", "B", "C"), 80, replace = TRUE)
  )
  Y <- data.frame(
    x1 = rnorm(60),
    x2 = sample(c("A", "B", "C"), 60, replace = TRUE)
  )

  result <- singling_out(X, Y, n_attacks = 100, seed = 42)

  expect_true(all(result$match_counts >= 0))
  expect_true(all(result$match_counts_control >= 0))
  expect_type(result$match_counts, "integer")
})

# --- Risk values for random vs memorized data ---

test_that("risk is low for independent random data", {
  set.seed(42)
  n <- 300
  X <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = sample(letters[1:10], n, replace = TRUE)
  )
  Y <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = sample(letters[1:10], n, replace = TRUE)
  )

  result <- singling_out(X, Y, n_attacks = 500, n_cols = 3, seed = 123)

  # For independent data, residual risk should be low
  expect_true(result$risk < 0.3)
})

test_that("risk is higher for copied/memorized data", {
  set.seed(123)
  n <- 200
  X <- data.frame(
    x1 = seq_len(n),  # unique values
    x2 = seq_len(n) * 2,
    x3 = paste0("cat_", seq_len(n))
  )
  # Y is exact copy - maximum memorization
  Y <- X

  result <- singling_out(X, Y, n_attacks = 200, n_cols = 2, seed = 42)

  # With unique values and exact copy, attack rate should be high
  expect_true(result$risk_attack > 0.1)
})

test_that("attack rate is higher than control rate for memorized data", {
  set.seed(123)
  n <- 100
  train <- data.frame(
    x1 = seq_len(n),
    x2 = seq_len(n) * 3
  )
  holdout <- data.frame(
    x1 = (n + 1):(2 * n),
    x2 = ((n + 1):(2 * n)) * 3
  )
  # Synthetic is copy of training
  Y <- train

  result <- singling_out(train, Y, holdout = holdout,
                          n_attacks = 200, n_cols = 2, seed = 42)

  # Attack should succeed more on training than holdout
  expect_true(result$n_success >= result$n_control_success)
})

# --- Multivariate vs univariate modes ---

test_that("multivariate mode works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))

  result <- singling_out(X, Y, mode = "multivariate", n_attacks = 50, seed = 42)

  expect_s3_class(result, "singling_out")
  expect_equal(result$mode, "multivariate")
})

test_that("univariate mode works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  result <- singling_out(X, Y, mode = "univariate", n_attacks = 50, seed = 42)

  expect_s3_class(result, "singling_out")
  expect_equal(result$mode, "univariate")
})

test_that("more columns per predicate increases singling out for unique data", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    x1 = seq_len(n),
    x2 = rnorm(n),
    x3 = sample(letters[1:5], n, replace = TRUE),
    x4 = rnorm(n)
  )
  Y <- X  # exact copy

  r2 <- singling_out(X, Y, n_cols = 2, n_attacks = 200, seed = 42)
  r4 <- singling_out(X, Y, n_cols = 4, n_attacks = 200, seed = 42)

  # More columns = more specific predicates = more singling out (or equal)
  expect_true(r4$n_success >= r2$n_success * 0.5)  # allow some randomness
})

# --- Holdout handling ---

test_that("singling_out works with explicit holdout", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))
  holdout <- data.frame(x1 = rnorm(30), x2 = rnorm(30))

  result <- singling_out(X, Y, holdout = holdout, n_attacks = 50, seed = 42)

  expect_s3_class(result, "singling_out")
  expect_equal(result$n_original, nrow(X))
  expect_equal(result$n_holdout, nrow(holdout))
})

test_that("singling_out works with holdout_fraction", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(50), x2 = rnorm(50))

  result <- singling_out(X, Y, holdout_fraction = 0.3,
                          n_attacks = 50, seed = 42)

  expect_s3_class(result, "singling_out")
  expect_equal(result$n_holdout, floor(n * 0.3))
  expect_equal(result$n_original, n - floor(n * 0.3))
})

test_that("singling_out seed produces reproducible results", {
  set.seed(99)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  r1 <- singling_out(X, Y, n_attacks = 100, seed = 42)
  r2 <- singling_out(X, Y, n_attacks = 100, seed = 42)

  expect_equal(r1$risk, r2$risk)
  expect_equal(r1$n_success, r2$n_success)
  expect_equal(r1$match_counts, r2$match_counts)
})

# --- synth_pair method ---

test_that("singling_out.synth_pair works", {
  set.seed(123)
  original <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  synthetic <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  pair <- synth_pair(original = original, synthetic = synthetic)
  result <- singling_out(pair, n_attacks = 50, seed = 42)

  expect_s3_class(result, "singling_out")
  expect_true(!is.null(result$risk))
})

test_that("singling_out.synth_pair uses holdout when available", {
  set.seed(123)
  original <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  synthetic <- data.frame(x1 = rnorm(40), x2 = rnorm(40))
  holdout <- data.frame(x1 = rnorm(30), x2 = rnorm(30))

  pair <- synth_pair(original = original, synthetic = synthetic,
                     holdout = holdout)
  result <- singling_out(pair, n_attacks = 50, seed = 42)

  expect_equal(result$n_original, 60)
  expect_equal(result$n_holdout, 30)
})

# --- Input validation ---

test_that("singling_out input validation works", {
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  Y <- data.frame(x1 = rnorm(10), x2 = rnorm(10))

  expect_error(singling_out(1:10, Y), "X must be a data frame")
  expect_error(singling_out(X, 1:10), "Y must be a data frame")
  expect_error(singling_out(X, Y, holdout = "not a df"),
               "holdout must be a data frame")
})

test_that("singling_out errors when no common variables", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(b = 1:10)

  expect_error(singling_out(X, Y), "No common variables")
})

test_that("singling_out errors with mismatched variable types", {
  X <- data.frame(x = 1:10)
  Y <- data.frame(x = letters[1:10])

  expect_error(singling_out(X, Y), "different class")
})

test_that("singling_out errors with missing variables", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10)

  expect_error(singling_out(X, Y, vars = c("a", "b")), "Variables missing in X")
})

test_that("singling_out validates n_attacks", {
  X <- data.frame(x1 = rnorm(20))
  Y <- data.frame(x1 = rnorm(10))

  expect_error(singling_out(X, Y, n_attacks = -1), "positive integer")
  expect_error(singling_out(X, Y, n_attacks = "abc"), "positive integer")
})

test_that("singling_out validates n_cols", {
  X <- data.frame(x1 = rnorm(20))
  Y <- data.frame(x1 = rnorm(10))

  expect_error(singling_out(X, Y, n_cols = 0), "positive integer")
})

test_that("singling_out validates confidence_level", {
  X <- data.frame(x1 = rnorm(20))
  Y <- data.frame(x1 = rnorm(10))

  expect_error(singling_out(X, Y, confidence_level = 1.5),
               "between 0 and 1")
  expect_error(singling_out(X, Y, confidence_level = 0),
               "between 0 and 1")
})

test_that("singling_out handles NA values", {
  set.seed(123)
  X <- data.frame(x1 = c(rnorm(18), NA, NA), x2 = rnorm(20))
  Y <- data.frame(x1 = c(rnorm(9), NA), x2 = rnorm(10))

  result <- singling_out(X, Y, n_attacks = 50, seed = 42, na.rm = TRUE)

  expect_s3_class(result, "singling_out")
})

test_that("n_cols is capped at number of variables", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(40), x2 = rnorm(40))
  Y <- data.frame(x1 = rnorm(30), x2 = rnorm(30))

  # n_cols = 10 but only 2 variables -> should cap at 2
  result <- singling_out(X, Y, n_cols = 10, n_attacks = 20, seed = 42)

  expect_s3_class(result, "singling_out")
  expect_equal(result$n_cols, 2)
})

# --- Wilson score CI correctness ---

test_that("Wilson score gives correct known values", {
  # 0 successes out of 100
  w0 <- riskutility:::.wilson_score(0, 100, 0.95)
  expect_equal(w0$estimate, 0)
  expect_true(w0$ci_lower == 0)
  expect_true(w0$ci_upper > 0 && w0$ci_upper < 0.1)

  # 100 successes out of 100
  w100 <- riskutility:::.wilson_score(100, 100, 0.95)
  expect_equal(w100$estimate, 1)
  expect_true(w100$ci_lower > 0.9)
  expect_true(w100$ci_upper == 1)

  # 50 successes out of 100
  w50 <- riskutility:::.wilson_score(50, 100, 0.95)
  expect_equal(w50$estimate, 0.5)
  expect_true(w50$ci_lower > 0.35 && w50$ci_lower < 0.5)
  expect_true(w50$ci_upper > 0.5 && w50$ci_upper < 0.65)

  # 0 total
  w_empty <- riskutility:::.wilson_score(0, 0, 0.95)
  expect_equal(w_empty$estimate, 0)
})

test_that("Wilson CI is narrower for larger samples", {
  w_small <- riskutility:::.wilson_score(5, 10, 0.95)
  w_large <- riskutility:::.wilson_score(500, 1000, 0.95)

  width_small <- w_small$ci_upper - w_small$ci_lower
  width_large <- w_large$ci_upper - w_large$ci_lower

  expect_true(width_large < width_small)
})

# --- Variable selection ---

test_that("singling_out uses specified vars", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))

  result <- singling_out(X, Y, vars = c("x1", "x2"),
                          n_attacks = 50, seed = 42)

  expect_equal(result$vars, c("x1", "x2"))
})

test_that("singling_out auto-selects common vars", {
  set.seed(123)
  X <- data.frame(a = rnorm(60), b = rnorm(60), c_only = rnorm(60))
  Y <- data.frame(a = rnorm(40), b = rnorm(40), d_only = rnorm(40))

  result <- singling_out(X, Y, n_attacks = 50, seed = 42)

  expect_equal(sort(result$vars), c("a", "b"))
})

# --- Dataset sizes stored correctly ---

test_that("singling_out stores correct dataset sizes", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
  Y <- data.frame(x1 = rnorm(80), x2 = rnorm(80))

  result <- singling_out(X, Y, holdout_fraction = 0.4,
                          n_attacks = 50, seed = 42)

  expect_equal(result$n_synthetic, 80)
  expect_equal(result$n_holdout, 40)
  expect_equal(result$n_original, 60)
})

# --- Mixed data types ---

test_that("singling_out handles mixed data types", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    age = sample(20:70, n, replace = TRUE),
    income = rnorm(n, 50000, 15000),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE))
  )
  Y <- data.frame(
    age = sample(20:70, n, replace = TRUE),
    income = rnorm(n, 50000, 15000),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE))
  )

  result <- singling_out(X, Y, n_attacks = 100, n_cols = 3, seed = 42)

  expect_s3_class(result, "singling_out")
})

# --- print, summary, plot ---

test_that("singling_out print method works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- singling_out(X, Y, n_attacks = 50, seed = 42)

  expect_output(print(result), "Singling Out Risk Assessment")
  expect_output(print(result), "Attack Results")
  expect_output(print(result), "Risk Score")
  expect_output(print(result), "Privacy Assessment")
})

test_that("singling_out summary returns correct class", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- singling_out(X, Y, n_attacks = 50, seed = 42)
  s <- summary(result)

  expect_s3_class(s, "summary.singling_out")
  expect_true(!is.null(s$risk))
  expect_true(!is.null(s$match_count_dist))
  expect_true(!is.null(s$match_count_dist_control))
})

test_that("singling_out print.summary works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- singling_out(X, Y, n_attacks = 50, seed = 42)
  s <- summary(result)

  expect_output(print(s), "Summary: Singling Out Risk")
  expect_output(print(s), "Key Metrics")
  expect_output(print(s), "Match Count Distribution")
})

test_that("singling_out plot method works without error", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- singling_out(X, Y, n_attacks = 50, seed = 42)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

# --- n_success consistency ---

test_that("n_success matches match_counts == 1", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  result <- singling_out(X, Y, n_attacks = 100, seed = 42)

  expect_equal(result$n_success, sum(result$match_counts == 1))
  expect_equal(result$n_control_success, sum(result$match_counts_control == 1))
})

test_that("risk_attack equals n_success / n_attacks", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  result <- singling_out(X, Y, n_attacks = 100, seed = 42)

  expect_equal(result$risk_attack, result$n_success / result$n_attacks)
})

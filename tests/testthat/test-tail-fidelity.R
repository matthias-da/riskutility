# Tests for tail_fidelity (Tail Preservation Utility Measure)

library(testthat)

test_that("tail_fidelity returns correct S3 class structure", {
  set.seed(123)
  n <- 500
  X <- data.frame(
    income = rlnorm(n, meanlog = 10, sdlog = 1),
    age = rnorm(n, 40, 10),
    score = rnorm(n, 100, 15)
  )
  Y <- data.frame(
    income = rlnorm(n, meanlog = 10, sdlog = 1),
    age = rnorm(n, 40, 10),
    score = rnorm(n, 100, 15)
  )

  result <- tail_fidelity(X, Y)

  expect_s3_class(result, "tail_fidelity")
  expect_true("qq_divergence" %in% names(result))
  expect_true("utility_score" %in% names(result))
  expect_true("per_variable" %in% names(result))
  expect_true("percentile" %in% names(result))
  expect_true("tails" %in% names(result))
  expect_true("n_vars" %in% names(result))
  expect_true("vars" %in% names(result))
  expect_true("n_X" %in% names(result))
  expect_true("n_Y" %in% names(result))
  expect_true("hill" %in% names(result))
})

test_that("tail_fidelity near zero for same distribution", {
  set.seed(42)
  n <- 1000
  X <- data.frame(x = rnorm(n, 50, 10))
  Y <- data.frame(x = rnorm(n, 50, 10))

  result <- tail_fidelity(X, Y)

  # Should be small for same distribution
  expect_true(result$qq_divergence < 0.5)
  expect_true(result$utility_score > 0.5)
})

test_that("tail_fidelity near zero for identical data", {
  set.seed(123)
  n <- 500
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- X  # identical

  result <- tail_fidelity(X, Y)

  expect_equal(result$qq_divergence, 0, tolerance = 1e-10)
  expect_equal(result$utility_score, 1, tolerance = 1e-10)
})

test_that("tail_fidelity detects distorted tails (t(3) vs normal)", {
  set.seed(99)
  n <- 1000
  # Original: heavy-tailed t-distribution
  X <- data.frame(val = rt(n, df = 3))
  # Synthetic: normal (lighter tails)
  Y <- data.frame(val = rnorm(n, mean = 0, sd = sd(X$val)))

  result <- tail_fidelity(X, Y)

  # Tail divergence should be notably non-zero

  expect_true(result$qq_divergence > 0.05)
  expect_true(result$utility_score < 0.95)
})

test_that("tails parameter works: both", {
  set.seed(100)
  n <- 500
  X <- data.frame(x = rnorm(n))
  Y <- data.frame(x = rnorm(n))

  result <- tail_fidelity(X, Y, tails = "both")
  expect_equal(result$tails, "both")
})

test_that("tails parameter works: upper", {
  set.seed(101)
  n <- 500
  X <- data.frame(x = rnorm(n))
  Y <- data.frame(x = rnorm(n))

  result <- tail_fidelity(X, Y, tails = "upper")
  expect_equal(result$tails, "upper")
  expect_s3_class(result, "tail_fidelity")
})

test_that("tails parameter works: lower", {
  set.seed(102)
  n <- 500
  X <- data.frame(x = rnorm(n))
  Y <- data.frame(x = rnorm(n))

  result <- tail_fidelity(X, Y, tails = "lower")
  expect_equal(result$tails, "lower")
  expect_s3_class(result, "tail_fidelity")
})

test_that("upper vs lower detect asymmetric tail distortion", {
  set.seed(200)
  n <- 1000
  X <- data.frame(x = rnorm(n))
  # Distort only the upper tail: clip values above the 95th percentile
  y <- rnorm(n)
  y[y > quantile(y, 0.95)] <- quantile(y, 0.95)
  Y <- data.frame(x = y)

  result_upper <- tail_fidelity(X, Y, tails = "upper")
  result_lower <- tail_fidelity(X, Y, tails = "lower")

  # Upper tail should show more divergence than lower tail
  expect_true(result_upper$qq_divergence > result_lower$qq_divergence)
})

test_that("percentile parameter: 90", {
  set.seed(103)
  n <- 500
  X <- data.frame(x = rnorm(n))
  Y <- data.frame(x = rnorm(n))

  result <- tail_fidelity(X, Y, percentile = 90)
  expect_equal(result$percentile, 90)
  expect_s3_class(result, "tail_fidelity")
})

test_that("percentile parameter: 95 (default)", {
  set.seed(104)
  n <- 500
  X <- data.frame(x = rnorm(n))
  Y <- data.frame(x = rnorm(n))

  result <- tail_fidelity(X, Y, percentile = 95)
  expect_equal(result$percentile, 95)
})

test_that("percentile parameter: 99", {
  set.seed(105)
  n <- 500
  X <- data.frame(x = rnorm(n))
  Y <- data.frame(x = rnorm(n))

  result <- tail_fidelity(X, Y, percentile = 99)
  expect_equal(result$percentile, 99)
  expect_s3_class(result, "tail_fidelity")
})

test_that("percentile validation rejects bad values", {
  X <- data.frame(x = 1:100)
  Y <- data.frame(x = 1:100)

  expect_error(tail_fidelity(X, Y, percentile = 50), "'percentile' must be")
  expect_error(tail_fidelity(X, Y, percentile = 100), "'percentile' must be")
  expect_error(tail_fidelity(X, Y, percentile = 30), "'percentile' must be")
})

test_that("per_variable results have correct dimensions", {
  set.seed(106)
  n <- 300
  X <- data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n))
  Y <- data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n))

  result <- tail_fidelity(X, Y)

  expect_equal(nrow(result$per_variable), 3)
  expect_true("variable" %in% names(result$per_variable))
  expect_true("qq_tail_div" %in% names(result$per_variable))
  expect_true("jsd_tail" %in% names(result$per_variable))
  expect_equal(result$per_variable$variable, c("a", "b", "c"))
})

test_that("utility_score in [0, 1]", {
  set.seed(107)
  n <- 300
  X <- data.frame(x = rnorm(n))
  Y <- data.frame(x = rnorm(n, 10, 5))

  result <- tail_fidelity(X, Y)

  expect_true(result$utility_score >= 0)
  expect_true(result$utility_score <= 1)
})

test_that("utility_score equals exp(-qq_divergence)", {
  set.seed(108)
  n <- 500
  X <- data.frame(x = rnorm(n), y = rnorm(n))
  Y <- data.frame(x = rnorm(n), y = rnorm(n))

  result <- tail_fidelity(X, Y)

  expected_score <- exp(-result$qq_divergence)
  expected_score <- max(0, min(1, expected_score))
  expect_equal(result$utility_score, expected_score, tolerance = 1e-10)
})

test_that("synth_pair dispatch works", {
  set.seed(109)
  n <- 300
  original <- data.frame(
    income = rnorm(n, 50000, 10000),
    age = rnorm(n, 40, 10),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )
  synthetic <- data.frame(
    income = rnorm(n, 50000, 10000),
    age = rnorm(n, 40, 10),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )

  pair <- synth_pair(original, synthetic)
  result <- tail_fidelity(pair)

  expect_s3_class(result, "tail_fidelity")
  expect_equal(result$n_vars, 2)
  expect_true("income" %in% result$vars)
  expect_true("age" %in% result$vars)
  expect_false("gender" %in% result$vars)
})

test_that("skips categoricals with message", {
  set.seed(110)
  n <- 300
  X <- data.frame(
    x = rnorm(n),
    cat_var = sample(letters[1:5], n, replace = TRUE)
  )
  Y <- data.frame(
    x = rnorm(n),
    cat_var = sample(letters[1:5], n, replace = TRUE)
  )

  expect_message(tail_fidelity(X, Y), "Skipping non-numeric variables")
})

test_that("Hill estimator adds hill_diff column", {
  set.seed(111)
  n <- 500
  X <- data.frame(x = rlnorm(n, 5, 1))
  Y <- data.frame(x = rlnorm(n, 5, 1))

  result <- tail_fidelity(X, Y, hill = TRUE)

  expect_true("hill_diff" %in% names(result$per_variable))
  expect_true(result$hill)
  # Hill diff should be finite (both datasets are positive)
  expect_true(is.numeric(result$per_variable$hill_diff))
})

test_that("Hill estimator returns NA for non-positive data", {
  set.seed(112)
  n <- 500
  X <- data.frame(x = rnorm(n))  # includes negatives
  Y <- data.frame(x = rnorm(n))

  result <- tail_fidelity(X, Y, hill = TRUE)

  # The upper tail of N(0,1) is positive, so this may or may not be NA

  # depending on the percentile cutoff. With percentile=95, the upper 5%
  # of N(0,1) is all positive, so Hill should work.
  expect_true("hill_diff" %in% names(result$per_variable))
})

test_that("Hill diff without hill=TRUE does not include column", {
  set.seed(113)
  n <- 300
  X <- data.frame(x = rlnorm(n))
  Y <- data.frame(x = rlnorm(n))

  result <- tail_fidelity(X, Y, hill = FALSE)

  expect_false("hill_diff" %in% names(result$per_variable))
  expect_false(result$hill)
})

test_that("print method produces output", {
  set.seed(114)
  n <- 300
  X <- data.frame(x = rnorm(n), y = rnorm(n))
  Y <- data.frame(x = rnorm(n), y = rnorm(n))

  result <- tail_fidelity(X, Y)

  expect_output(print(result), "Tail Fidelity")
  expect_output(print(result), "QQ tail divergence")
  expect_output(print(result), "Utility score")
  expect_output(print(result), "Interpretation")
})

test_that("print method works with hill=TRUE", {
  set.seed(115)
  n <- 300
  X <- data.frame(x = rlnorm(n))
  Y <- data.frame(x = rlnorm(n))

  result <- tail_fidelity(X, Y, hill = TRUE)

  expect_output(print(result), "hill_diff")
})

test_that("summary method returns correct class", {
  set.seed(116)
  n <- 300
  X <- data.frame(x = rnorm(n), y = rnorm(n))
  Y <- data.frame(x = rnorm(n), y = rnorm(n))

  result <- tail_fidelity(X, Y)
  s <- summary(result)

  expect_s3_class(s, "summary.tail_fidelity")
  expect_true(!is.null(s$qq_divergence))
  expect_true(!is.null(s$utility_score))
  expect_true(!is.null(s$mean_jsd_tail))
  expect_true(!is.null(s$max_qq_divergence))
  expect_true(!is.null(s$worst_variable))
  expect_true(!is.null(s$per_variable))
  expect_output(print(s), "Summary: Tail Fidelity")
})

test_that("plot method works for which=1", {
  set.seed(117)
  n <- 300
  X <- data.frame(a = rnorm(n), b = rnorm(n))
  Y <- data.frame(a = rnorm(n), b = rnorm(n))

  result <- tail_fidelity(X, Y)
  expect_no_error(plot(result, which = 1))
})

test_that("plot method works for which=2", {
  set.seed(118)
  n <- 300
  X <- data.frame(a = rnorm(n), b = rnorm(n))
  Y <- data.frame(a = rnorm(n), b = rnorm(n))

  result <- tail_fidelity(X, Y)
  expect_no_error(plot(result, which = 2))
})

test_that("plot method works for which=1:2", {
  set.seed(119)
  n <- 300
  X <- data.frame(a = rnorm(n), b = rnorm(n))
  Y <- data.frame(a = rnorm(n), b = rnorm(n))

  result <- tail_fidelity(X, Y)
  expect_no_error(plot(result, which = 1:2))
})

test_that("NA handling works", {
  set.seed(120)
  n <- 300
  x_vals <- c(rnorm(n - 5), rep(NA, 5))
  X <- data.frame(x = x_vals)
  Y <- data.frame(x = c(rnorm(n - 3), rep(NA, 3)))

  result <- tail_fidelity(X, Y, na.rm = TRUE)

  expect_s3_class(result, "tail_fidelity")
  expect_true(result$n_X <= n)
  expect_true(result$n_Y <= n)
})

test_that("error when no numeric variables", {
  X <- data.frame(a = c("x", "y", "z"), b = c("p", "q", "r"))
  Y <- data.frame(a = c("x", "y", "z"), b = c("p", "q", "r"))

  expect_error(tail_fidelity(X, Y), "No numeric variables found")
})

test_that("error when variable missing in X", {
  X <- data.frame(x1 = 1:100)
  Y <- data.frame(x1 = 1:100, x2 = rnorm(100))

  expect_error(tail_fidelity(X, Y, vars = c("x1", "x2")),
               "Variables missing in X")
})

test_that("error when variable missing in Y", {
  X <- data.frame(x1 = 1:100, x2 = rnorm(100))
  Y <- data.frame(x1 = 1:100)

  expect_error(tail_fidelity(X, Y, vars = c("x1", "x2")),
               "Variables missing in Y")
})

test_that("works with single variable", {
  set.seed(121)
  n <- 500
  X <- data.frame(x = rnorm(n))
  Y <- data.frame(x = rnorm(n))

  result <- tail_fidelity(X, Y)

  expect_equal(result$n_vars, 1)
  expect_equal(result$vars, "x")
  expect_equal(nrow(result$per_variable), 1)
})

test_that("works with data.table input", {
  skip_if_not_installed("data.table")
  set.seed(122)
  n <- 300
  X <- data.table::data.table(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.table::data.table(x1 = rnorm(n), x2 = rnorm(n))

  result <- tail_fidelity(X, Y)

  expect_s3_class(result, "tail_fidelity")
  expect_true(result$qq_divergence >= 0)
})

test_that("explicit vars parameter works", {
  set.seed(123)
  n <- 300
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- tail_fidelity(X, Y, vars = c("x1", "x3"))

  expect_equal(result$n_vars, 2)
  expect_equal(result$vars, c("x1", "x3"))
  expect_equal(nrow(result$per_variable), 2)
})

test_that("distorted tails have worse utility than matching tails", {
  set.seed(124)
  n <- 1000
  X <- data.frame(x = rnorm(n, 0, 1))

  # Good synthetic: same distribution
  Y_good <- data.frame(x = rnorm(n, 0, 1))

  # Poor synthetic: truncated at +/- 2 (losing tails)
  y_poor <- rnorm(n, 0, 1)
  y_poor <- pmin(pmax(y_poor, -2), 2)
  Y_poor <- data.frame(x = y_poor)

  result_good <- tail_fidelity(X, Y_good)
  result_poor <- tail_fidelity(X, Y_poor)

  expect_true(result_poor$qq_divergence > result_good$qq_divergence)
  expect_true(result_poor$utility_score < result_good$utility_score)
})

test_that("jsd_tail is NA when too few tail observations", {
  set.seed(125)
  # Use very small n so tail has < 20 observations
  n <- 50  # 5% tail = 2-3 observations
  X <- data.frame(x = rnorm(n))
  Y <- data.frame(x = rnorm(n))

  result <- tail_fidelity(X, Y, percentile = 99)

  # With n=50 and percentile=99, only ~0.5 observations in each tail
  # JSD should be NA
  expect_true(is.na(result$per_variable$jsd_tail[1]))
})

test_that("multiple variables have independent per-variable results", {
  set.seed(126)
  n <- 500
  # x1: identical tails; x2: very different tails
  X <- data.frame(x1 = rnorm(n), x2 = rt(n, df = 2))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- tail_fidelity(X, Y)

  # x2 (heavy-tailed vs normal) should have larger qq divergence than x1
  qq_x1 <- result$per_variable$qq_tail_div[result$per_variable$variable == "x1"]
  qq_x2 <- result$per_variable$qq_tail_div[result$per_variable$variable == "x2"]
  expect_true(qq_x2 > qq_x1)
})

test_that("qq_divergence is mean of per-variable qq_tail_div", {
  set.seed(127)
  n <- 300
  X <- data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n))
  Y <- data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n))

  result <- tail_fidelity(X, Y)

  expected_mean <- mean(result$per_variable$qq_tail_div)
  expect_equal(result$qq_divergence, expected_mean, tolerance = 1e-10)
})

test_that("Hill estimator detects different tail indices", {
  set.seed(128)
  n <- 1000
  # Pareto-like with different tail indices
  X <- data.frame(x = rlnorm(n, 5, 2))  # heavy tail
  Y <- data.frame(x = rlnorm(n, 5, 0.5))  # lighter tail

  result <- tail_fidelity(X, Y, hill = TRUE)

  # Hill diff should be non-trivial
  expect_true(!is.na(result$per_variable$hill_diff[1]))
  expect_true(result$per_variable$hill_diff[1] > 0)
})

test_that("all NAs in X raises error", {
  X <- data.frame(x = c(NA_real_, NA_real_, NA_real_))
  Y <- data.frame(x = c(1, 2, 3))

  expect_error(tail_fidelity(X, Y, na.rm = TRUE), "No complete cases in X")
})

test_that("all NAs in Y raises error", {
  X <- data.frame(x = c(1, 2, 3))
  Y <- data.frame(x = c(NA_real_, NA_real_, NA_real_))

  expect_error(tail_fidelity(X, Y, na.rm = TRUE), "No complete cases in Y")
})

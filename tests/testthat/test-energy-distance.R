# Tests for energy_distance (Energy Distance for Multivariate Numeric Data)

library(testthat)

test_that("energy_distance returns correct S3 class structure", {
  set.seed(123)
  n <- 80
  X <- data.frame(
    income = rnorm(n, 50000, 10000),
    age = rnorm(n, 40, 10),
    score = rnorm(n, 100, 15)
  )
  Y <- data.frame(
    income = rnorm(n, 50000, 10000),
    age = rnorm(n, 40, 10),
    score = rnorm(n, 100, 15)
  )

  result <- energy_distance(X, Y, seed = 42)

  expect_s3_class(result, "energy_distance")
  expect_true("energy_distance" %in% names(result))
  expect_true("energy_distance_normalized" %in% names(result))
  expect_true("mean_dist_XY" %in% names(result))
  expect_true("mean_dist_XX" %in% names(result))
  expect_true("mean_dist_YY" %in% names(result))
  expect_true("n_X" %in% names(result))
  expect_true("n_Y" %in% names(result))
  expect_true("n_vars" %in% names(result))
  expect_true("vars" %in% names(result))
  expect_true("standardized" %in% names(result))
  expect_true("utility_score" %in% names(result))
  expect_true("sampled" %in% names(result))
})

test_that("energy_distance is non-negative", {
  set.seed(123)
  n <- 60
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- energy_distance(X, Y, seed = 42)

  expect_true(result$energy_distance >= 0)
  expect_true(result$energy_distance_normalized >= 0)
})

test_that("energy_distance is zero (or near zero) for identical data", {
  set.seed(123)
  n <- 50
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- X  # Identical

  result <- energy_distance(X, Y, seed = 42)

  expect_equal(result$energy_distance, 0, tolerance = 1e-10)
})

test_that("energy_distance is small for similar distributions", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    income = rnorm(n, 50000, 10000),
    age = rnorm(n, 40, 10)
  )
  Y <- data.frame(
    income = rnorm(n, 50000, 10000),
    age = rnorm(n, 40, 10)
  )

  result <- energy_distance(X, Y, seed = 42)

  # Energy distance should be small for data from the same distribution
  expect_true(result$energy_distance < 0.5)
})

test_that("energy_distance is larger for different distributions", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    x1 = rnorm(n, 0, 1),
    x2 = rnorm(n, 0, 1)
  )
  Y <- data.frame(
    x1 = rnorm(n, 5, 1),
    x2 = rnorm(n, 5, 1)
  )

  result_diff <- energy_distance(X, Y, seed = 42)

  # Same distribution for comparison
  Y_same <- data.frame(
    x1 = rnorm(n, 0, 1),
    x2 = rnorm(n, 0, 1)
  )
  result_same <- energy_distance(X, Y_same, seed = 42)

  expect_true(result_diff$energy_distance > result_same$energy_distance)
})

test_that("energy_distance utility_score equals exp(-energy_distance)", {
  set.seed(123)
  n <- 60
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- energy_distance(X, Y, seed = 42)

  expect_equal(result$utility_score, exp(-result$energy_distance), tolerance = 1e-10)
})

test_that("energy_distance works with a single variable", {
  set.seed(123)
  n <- 80
  X <- data.frame(x = rnorm(n, 10, 2))
  Y <- data.frame(x = rnorm(n, 10, 2))

  result <- energy_distance(X, Y, seed = 42)

  expect_s3_class(result, "energy_distance")
  expect_equal(result$n_vars, 1)
  expect_equal(result$vars, "x")
})

test_that("energy_distance auto-detects numeric variables", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    income = rnorm(n, 50000, 10000),
    age = rnorm(n, 40, 10),
    gender = sample(c("M", "F"), n, replace = TRUE)  # should be excluded
  )
  Y <- data.frame(
    income = rnorm(n, 50000, 10000),
    age = rnorm(n, 40, 10),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )

  result <- energy_distance(X, Y, seed = 42)

  expect_equal(result$n_vars, 2)
  expect_true("income" %in% result$vars)
  expect_true("age" %in% result$vars)
  expect_false("gender" %in% result$vars)
})

test_that("energy_distance works with explicit vars parameter", {
  set.seed(123)
  n <- 60
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- energy_distance(X, Y, vars = c("x1", "x2"), seed = 42)

  expect_equal(result$n_vars, 2)
  expect_equal(result$vars, c("x1", "x2"))
})

test_that("energy_distance errors when no numeric variables found", {
  X <- data.frame(a = c("x", "y"), b = c("p", "q"))
  Y <- data.frame(a = c("x", "y"), b = c("p", "q"))

  expect_error(energy_distance(X, Y), "No numeric variables found or specified")
})

test_that("energy_distance errors when variables missing in X", {
  X <- data.frame(x1 = 1:10)
  Y <- data.frame(x1 = 1:10, x2 = rnorm(10))

  expect_error(energy_distance(X, Y, vars = c("x1", "x2")), "Variables missing in X")
})

test_that("energy_distance errors when variables missing in Y", {
  X <- data.frame(x1 = 1:10, x2 = rnorm(10))
  Y <- data.frame(x1 = 1:10)

  expect_error(energy_distance(X, Y, vars = c("x1", "x2")), "Variables missing in Y")
})

test_that("energy_distance handles NA values with na.rm = TRUE", {
  set.seed(123)
  n <- 50
  X <- data.frame(x1 = c(rnorm(n - 2), NA, NA), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = c(rnorm(n - 1), NA))

  result <- energy_distance(X, Y, na.rm = TRUE, seed = 42)

  expect_s3_class(result, "energy_distance")
  # Number of used records should be less than original
  expect_true(result$n_X <= n)
  expect_true(result$n_Y <= n)
})

test_that("energy_distance errors when all cases have NA", {
  X <- data.frame(x1 = c(NA_real_, NA_real_, NA_real_))
  Y <- data.frame(x1 = c(1, 2, 3))

  expect_error(energy_distance(X, Y, na.rm = TRUE), "No complete cases in X")
})

test_that("energy_distance standardize parameter works", {
  set.seed(123)
  n <- 80
  # Use very different scales and shifted means so standardization matters
  X <- data.frame(x1 = rnorm(n, 1000, 100), x2 = rnorm(n, 0, 1))
  Y <- data.frame(x1 = rnorm(n, 1050, 100), x2 = rnorm(n, 0.5, 1))

  result_std <- energy_distance(X, Y, standardize = TRUE, seed = 42)
  result_raw <- energy_distance(X, Y, standardize = FALSE, seed = 42)

  expect_true(result_std$standardized)
  expect_false(result_raw$standardized)
  # With different scales, standardization should change the result
  expect_false(isTRUE(all.equal(result_std$energy_distance, result_raw$energy_distance)))
})

test_that("energy_distance sampling works for large datasets", {
  set.seed(123)
  n <- 1500  # Larger than default n_sample = 1000
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- energy_distance(X, Y, n_sample = 100, seed = 42)

  expect_s3_class(result, "energy_distance")
  expect_true(result$sampled)
  expect_equal(result$n_X, 100)
  expect_equal(result$n_Y, 100)
  expect_equal(result$n_X_original, n)
  expect_equal(result$n_Y_original, n)
})

test_that("energy_distance with n_sample = NULL uses all data", {
  set.seed(123)
  n <- 80
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- energy_distance(X, Y, n_sample = NULL, seed = 42)

  expect_false(result$sampled)
  expect_equal(result$n_X, n)
  expect_equal(result$n_Y, n)
})

test_that("energy_distance seed parameter ensures reproducibility", {
  set.seed(1)
  n <- 1500
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result1 <- energy_distance(X, Y, n_sample = 100, seed = 42)
  result2 <- energy_distance(X, Y, n_sample = 100, seed = 42)

  expect_equal(result1$energy_distance, result2$energy_distance)
})

test_that("energy_distance print method works", {
  set.seed(123)
  n <- 60
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- energy_distance(X, Y, seed = 42)

  expect_output(print(result), "Energy Distance")
  expect_output(print(result), "Normalized")
  expect_output(print(result), "Utility")
  expect_output(print(result), "Interpretation")
})

test_that("energy_distance summary method returns correct class", {
  set.seed(123)
  n <- 60
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- energy_distance(X, Y, seed = 42)
  s <- summary(result)

  expect_s3_class(s, "summary.energy_distance")
  expect_true(!is.null(s$energy_distance))
  expect_true(!is.null(s$energy_distance_normalized))
  expect_true(!is.null(s$between_within_ratio))
  expect_true(!is.null(s$utility_score))
  expect_true(!is.null(s$vars))
  expect_output(print(s), "Summary: Energy Distance")
})

test_that("energy_distance plot method works for all which values", {
  set.seed(123)
  n <- 60
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- energy_distance(X, Y, seed = 42)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

test_that("energy_distance works via synth_pair dispatch", {
  set.seed(123)
  n <- 80
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
  result <- energy_distance(pair, seed = 42)

  expect_s3_class(result, "energy_distance")
  # Should only use the numeric variables detected in synth_pair
  expect_equal(result$n_vars, 2)
  expect_true("income" %in% result$vars)
  expect_true("age" %in% result$vars)
})

test_that("energy_distance distance components are consistent", {
  set.seed(123)
  n <- 60
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- energy_distance(X, Y, seed = 42)

  # Energy distance = 2 * E[d(X,Y)] - E[d(X,X')] - E[d(Y,Y')]
  expected_ed <- 2 * result$mean_dist_XY - result$mean_dist_XX - result$mean_dist_YY
  expected_ed <- max(0, expected_ed)

  expect_equal(result$energy_distance, expected_ed, tolerance = 1e-10)
})

test_that("energy_distance works with data.table input", {
  skip_if_not_installed("data.table")
  set.seed(123)
  n <- 60
  X <- data.table::data.table(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.table::data.table(x1 = rnorm(n), x2 = rnorm(n))

  result <- energy_distance(X, Y, seed = 42)

  expect_s3_class(result, "energy_distance")
  expect_true(result$energy_distance >= 0)
})

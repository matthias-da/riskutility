# Tests for copula_fidelity (Empirical Copula Dependence Comparison)

library(testthat)

test_that("copula_fidelity returns correct S3 class structure", {
  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- copula_fidelity(X, Y)

  expect_s3_class(result, "copula_fidelity")
  expect_true("mean_cvm" %in% names(result))
  expect_true("utility_score" %in% names(result))
  expect_true("pairwise" %in% names(result))
  expect_true("n_vars" %in% names(result))
  expect_true("vars" %in% names(result))
  expect_true("n_X" %in% names(result))
  expect_true("n_Y" %in% names(result))
  expect_true("n_grid" %in% names(result))
})

test_that("near zero CvM for identical copulas (independent variables, same distribution)", {
  set.seed(42)
  n <- 500
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- copula_fidelity(X, Y)

  # Independent variables from the same distribution should have very small CvM
  expect_true(result$mean_cvm < 0.005)
  # Utility should be high
  expect_true(result$utility_score > 0.5)
})

test_that("detects distorted dependence (correlated original, independent synthetic)", {
  set.seed(123)
  n <- 400

  # Original with strong dependence
  x1 <- rnorm(n)
  x2 <- 0.9 * x1 + rnorm(n, sd = 0.4)
  x3 <- -0.7 * x1 + rnorm(n, sd = 0.7)
  X <- data.frame(x1 = x1, x2 = x2, x3 = x3)

  # Synthetic with broken dependence (independent)
  Y_bad <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  # Synthetic that preserves dependence
  s1 <- rnorm(n)
  s2 <- 0.9 * s1 + rnorm(n, sd = 0.4)
  s3 <- -0.7 * s1 + rnorm(n, sd = 0.7)
  Y_good <- data.frame(x1 = s1, x2 = s2, x3 = s3)

  result_bad <- copula_fidelity(X, Y_bad)
  result_good <- copula_fidelity(X, Y_good)

  # Distorted dependence should have much larger CvM

  expect_true(result_bad$mean_cvm > result_good$mean_cvm)
  # Utility should be lower for broken dependence
  expect_true(result_bad$utility_score < result_good$utility_score)
})

test_that("pairwise matrix has correct dimensions (3 vars -> 3 pairs)", {
  set.seed(123)
  n <- 200
  X <- data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n))
  Y <- data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n))

  result <- copula_fidelity(X, Y)

  # 3 choose 2 = 3 pairs
  expect_equal(nrow(result$pairwise), 3)
  expect_equal(ncol(result$pairwise), 3)
  expect_true(all(c("var1", "var2", "cvm_distance") %in% names(result$pairwise)))

  # Check all pairs are present
  pair_set <- paste(result$pairwise$var1, result$pairwise$var2, sep = "-")
  expect_true("a-b" %in% pair_set)
  expect_true("a-c" %in% pair_set)
  expect_true("b-c" %in% pair_set)
})

test_that("pairwise dimensions with 4 vars -> 6 pairs", {
  set.seed(456)
  n <- 200
  X <- data.frame(v1 = rnorm(n), v2 = rnorm(n), v3 = rnorm(n), v4 = rnorm(n))
  Y <- data.frame(v1 = rnorm(n), v2 = rnorm(n), v3 = rnorm(n), v4 = rnorm(n))

  result <- copula_fidelity(X, Y)

  # 4 choose 2 = 6 pairs
  expect_equal(nrow(result$pairwise), 6)
})

test_that("skips categorical variables with message", {
  set.seed(123)
  n <- 200
  X <- data.frame(
    income = rnorm(n, 50000, 10000),
    age = rnorm(n, 40, 10),
    score = rnorm(n, 100, 15),
    gender = sample(c("M", "F"), n, replace = TRUE),
    status = sample(c("A", "B", "C"), n, replace = TRUE)
  )
  Y <- data.frame(
    income = rnorm(n, 50000, 10000),
    age = rnorm(n, 40, 10),
    score = rnorm(n, 100, 15),
    gender = sample(c("M", "F"), n, replace = TRUE),
    status = sample(c("A", "B", "C"), n, replace = TRUE)
  )

  expect_message(
    result <- copula_fidelity(X, Y),
    "Skipping non-numeric variables"
  )

  # Should only use the 3 numeric variables
  expect_equal(result$n_vars, 3)
  expect_true(all(c("income", "age", "score") %in% result$vars))
  expect_false("gender" %in% result$vars)
  expect_false("status" %in% result$vars)
})

test_that("utility_score is in [0, 1]", {
  set.seed(123)
  n <- 200

  # Similar data
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  result1 <- copula_fidelity(X, Y)

  expect_true(result1$utility_score >= 0)
  expect_true(result1$utility_score <= 1)

  # Very different data
  set.seed(456)
  x1 <- rnorm(n)
  x2 <- x1 + rnorm(n, sd = 0.1)
  X2 <- data.frame(x1 = x1, x2 = x2)
  Y2 <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  result2 <- copula_fidelity(X2, Y2)

  expect_true(result2$utility_score >= 0)
  expect_true(result2$utility_score <= 1)
})

test_that("utility_score formula is correct: 1 / (1 + mean_cvm * 100)", {
  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- copula_fidelity(X, Y)

  expected_utility <- 1 / (1 + result$mean_cvm * 100)
  expect_equal(result$utility_score, expected_utility, tolerance = 1e-12)
})

test_that("synth_pair dispatch works", {
  set.seed(123)
  n <- 200
  original <- data.frame(
    income = rnorm(n, 50000, 10000),
    age = rnorm(n, 40, 10),
    score = rnorm(n, 100, 15),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )
  synthetic <- data.frame(
    income = rnorm(n, 50000, 10000),
    age = rnorm(n, 40, 10),
    score = rnorm(n, 100, 15),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )

  pair <- synth_pair(original, synthetic)
  result <- copula_fidelity(pair)

  expect_s3_class(result, "copula_fidelity")
  # Should only use the numeric variables detected by synth_pair
  expect_equal(result$n_vars, 3)
  expect_true("income" %in% result$vars)
  expect_true("age" %in% result$vars)
  expect_true("score" %in% result$vars)
})

test_that("vars parameter restricts variables", {
  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n))

  result <- copula_fidelity(X, Y, vars = c("x1", "x3"))

  expect_equal(result$n_vars, 2)
  expect_equal(result$vars, c("x1", "x3"))
  # 2 choose 2 = 1 pair
  expect_equal(nrow(result$pairwise), 1)
  expect_equal(result$pairwise$var1[1], "x1")
  expect_equal(result$pairwise$var2[1], "x3")
})

test_that("error with < 2 numeric variables", {
  set.seed(123)
  n <- 50
  X <- data.frame(x1 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n))

  expect_error(copula_fidelity(X, Y), "At least 2 numeric variables")
})

test_that("error with only categorical variables", {
  X <- data.frame(a = c("x", "y", "z"), b = c("p", "q", "r"))
  Y <- data.frame(a = c("x", "y", "z"), b = c("p", "q", "r"))

  expect_error(copula_fidelity(X, Y), "At least 2 numeric variables")
})

test_that("error with 1 numeric + categoricals", {
  set.seed(123)
  n <- 50
  X <- data.frame(x1 = rnorm(n), cat = sample(letters[1:3], n, replace = TRUE))
  Y <- data.frame(x1 = rnorm(n), cat = sample(letters[1:3], n, replace = TRUE))

  expect_error(copula_fidelity(X, Y), "At least 2 numeric variables")
})

test_that("error when variables missing in X", {
  X <- data.frame(x1 = 1:10, x2 = rnorm(10))
  Y <- data.frame(x1 = 1:10, x2 = rnorm(10), x3 = rnorm(10))

  expect_error(copula_fidelity(X, Y, vars = c("x1", "x3")), "Variables missing in X")
})

test_that("error when variables missing in Y", {
  X <- data.frame(x1 = 1:10, x2 = rnorm(10), x3 = rnorm(10))
  Y <- data.frame(x1 = 1:10, x2 = rnorm(10))

  expect_error(copula_fidelity(X, Y, vars = c("x1", "x3")), "Variables missing in Y")
})

test_that("print method works", {
  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- copula_fidelity(X, Y)

  expect_output(print(result), "Copula Fidelity")
  expect_output(print(result), "Mean CvM distance")
  expect_output(print(result), "Utility score")
  expect_output(print(result), "Interpretation")
})

test_that("summary method returns correct class", {
  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- copula_fidelity(X, Y)
  s <- summary(result)

  expect_s3_class(s, "summary.copula_fidelity")
  expect_true(!is.null(s$mean_cvm))
  expect_true(!is.null(s$utility_score))
  expect_true(!is.null(s$max_cvm))
  expect_true(!is.null(s$min_cvm))
  expect_true(!is.null(s$sd_cvm))
  expect_true(!is.null(s$n_pairs))
  expect_true(!is.null(s$vars))
  expect_equal(s$n_pairs, 3)
})

test_that("print.summary.copula_fidelity works", {
  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- copula_fidelity(X, Y)
  s <- summary(result)

  expect_output(print(s), "Summary: Copula Fidelity")
  expect_output(print(s), "CvM Distance Summary")
  expect_output(print(s), "Utility score")
})

test_that("plot method runs without error", {
  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- copula_fidelity(X, Y)

  expect_no_error(plot(result, which = 1))
})

test_that("NA handling with na.rm = TRUE", {
  set.seed(123)
  n <- 200
  x1 <- c(rnorm(n - 5), rep(NA, 5))
  x2 <- c(rnorm(n - 3), rep(NA, 3))
  X <- data.frame(x1 = x1, x2 = x2)
  Y <- data.frame(x1 = rnorm(n), x2 = c(rnorm(n - 2), NA, NA))

  result <- copula_fidelity(X, Y, na.rm = TRUE)

  expect_s3_class(result, "copula_fidelity")
  # n_X should reflect removed NA rows
  expect_true(result$n_X < n)
  expect_true(result$n_Y < n)
})

test_that("error when all cases have NA", {
  X <- data.frame(x1 = c(NA_real_, NA_real_), x2 = c(1, 2))
  Y <- data.frame(x1 = c(1, 2), x2 = c(3, 4))

  expect_error(copula_fidelity(X, Y, na.rm = TRUE), "No complete cases in X")
})

test_that("n_grid parameter is respected", {
  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result_low <- copula_fidelity(X, Y, n_grid = 20)
  result_high <- copula_fidelity(X, Y, n_grid = 80)

  expect_equal(result_low$n_grid, 20)
  expect_equal(result_high$n_grid, 80)
  # Results should differ (different grid resolution)
  expect_false(isTRUE(all.equal(result_low$mean_cvm, result_high$mean_cvm)))
})

test_that("mean_cvm is the mean of pairwise distances", {
  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- copula_fidelity(X, Y)

  expect_equal(result$mean_cvm, mean(result$pairwise$cvm_distance), tolerance = 1e-12)
})

test_that("CvM distances are non-negative", {
  set.seed(789)
  n <- 300
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- copula_fidelity(X, Y)

  expect_true(all(result$pairwise$cvm_distance >= 0))
  expect_true(result$mean_cvm >= 0)
})

test_that("identical data gives CvM near zero", {
  set.seed(123)
  n <- 300
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- X  # Exact copy

  result <- copula_fidelity(X, Y)

  # CvM should be exactly zero for identical data (same ranks -> same empirical copula)
  expect_equal(result$mean_cvm, 0, tolerance = 1e-10)
  expect_equal(result$utility_score, 1, tolerance = 1e-10)
})

test_that("data.table input works", {
  skip_if_not_installed("data.table")
  set.seed(123)
  n <- 200
  X <- data.table::data.table(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.table::data.table(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- copula_fidelity(X, Y)

  expect_s3_class(result, "copula_fidelity")
  expect_true(result$mean_cvm >= 0)
})

test_that("different sample sizes between X and Y work", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(300), x2 = rnorm(300))
  Y <- data.frame(x1 = rnorm(200), x2 = rnorm(200))

  result <- copula_fidelity(X, Y)

  expect_s3_class(result, "copula_fidelity")
  expect_equal(result$n_X, 300)
  expect_equal(result$n_Y, 200)
})

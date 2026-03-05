# Tests for nndr (Nearest Neighbor Distance Ratio)

test_that("nndr returns correct S3 class structure", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    age = rnorm(n, 40, 10),
    income = rnorm(n, 50000, 15000)
  )
  Y <- data.frame(
    age = rnorm(80, 40, 10),
    income = rnorm(80, 50000, 15000)
  )

  result <- nndr(X, Y, method = "euclidean", seed = 42)

  expect_s3_class(result, "nndr")
  expect_true("nndr_train" %in% names(result))
  expect_true("nndr_holdout" %in% names(result))
  expect_true("nndr_ratio" %in% names(result))
  expect_true("mean_nndr_train" %in% names(result))
  expect_true("mean_nndr_holdout" %in% names(result))
  expect_true("n_suspicious" %in% names(result))
  expect_true("privacy_pass" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("vars" %in% names(result))
})

test_that("nndr values are non-negative", {
  set.seed(123)
  n <- 80
  X <- data.frame(
    age = rnorm(n, 40, 10),
    income = rnorm(n, 50000, 15000)
  )
  Y <- data.frame(
    age = rnorm(60, 40, 10),
    income = rnorm(60, 50000, 15000)
  )

  result <- nndr(X, Y, method = "euclidean", seed = 42)

  # Remove NAs before checking

  train_clean <- result$nndr_train[!is.na(result$nndr_train)]
  holdout_clean <- result$nndr_holdout[!is.na(result$nndr_holdout)]

  expect_true(all(train_clean >= 0))
  expect_true(all(holdout_clean >= 0))
  expect_true(result$mean_nndr_train >= 0)
  expect_true(result$mean_nndr_holdout >= 0)
})

test_that("nndr works with gower method", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    age = rnorm(n, 40, 10),
    gender = sample(c("M", "F"), n, replace = TRUE),
    income = rnorm(n, 50000, 15000)
  )
  Y <- data.frame(
    age = rnorm(40, 40, 10),
    gender = sample(c("M", "F"), 40, replace = TRUE),
    income = rnorm(40, 50000, 15000)
  )

  result <- nndr(X, Y, method = "gower", seed = 42)

  expect_s3_class(result, "nndr")
  expect_equal(result$method, "gower")
  expect_true(result$mean_nndr_train >= 0)
})

test_that("nndr works with euclidean method", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  Y <- data.frame(
    x1 = rnorm(40),
    x2 = rnorm(40)
  )

  result <- nndr(X, Y, method = "euclidean", seed = 42)

  expect_s3_class(result, "nndr")
  expect_equal(result$method, "euclidean")
})

test_that("nndr works with explicit holdout", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  Y <- data.frame(
    x1 = rnorm(40),
    x2 = rnorm(40)
  )
  holdout <- data.frame(
    x1 = rnorm(30),
    x2 = rnorm(30)
  )

  result <- nndr(X, Y, holdout = holdout, method = "euclidean")

  expect_s3_class(result, "nndr")
  # With explicit holdout, all of X is used as training
  expect_equal(result$n_train, nrow(X))
  expect_equal(result$n_holdout, nrow(holdout))
})

test_that("nndr works with holdout_fraction", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  Y <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50)
  )

  result <- nndr(X, Y, holdout_fraction = 0.3, method = "euclidean", seed = 42)

  expect_s3_class(result, "nndr")
  # holdout should be 30% of X
  expect_equal(result$n_holdout, floor(n * 0.3))
  expect_equal(result$n_train, n - floor(n * 0.3))
})

test_that("nndr print method works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- nndr(X, Y, method = "euclidean", seed = 42)

  expect_output(print(result), "Nearest Neighbor Distance Ratio")
  expect_output(print(result), "Mean NNDR")
  expect_output(print(result), "Privacy Assessment")
})

test_that("nndr summary method returns correct class", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- nndr(X, Y, method = "euclidean", seed = 42)
  s <- summary(result)

  expect_s3_class(s, "summary.nndr")
  expect_true(!is.null(s$nndr_ratio))
  expect_true(!is.null(s$privacy_pass))
  expect_true(!is.null(s$quantiles_train))
  expect_true(!is.null(s$quantiles_holdout))
  expect_output(print(s), "Summary: Nearest Neighbor Distance Ratio")
})

test_that("nndr input validation works", {
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  Y <- data.frame(x1 = rnorm(10), x2 = rnorm(10))

  expect_error(nndr(1:10, Y), "X must be a data frame")
  expect_error(nndr(X, 1:10), "Y must be a data frame")
  expect_error(nndr(X, Y, holdout = "not a df"), "holdout must be a data frame")
})

test_that("nndr errors when no common variables", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(b = 1:10)

  expect_error(nndr(X, Y), "No common variables")
})

test_that("nndr euclidean errors with non-numeric data", {
  X <- data.frame(x = c("a", "b", "c", "d", "e"))
  Y <- data.frame(x = c("a", "b", "c"))

  expect_error(nndr(X, Y, method = "euclidean"), "numeric")
})

test_that("nndr NNDR ratio is computed correctly", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- nndr(X, Y, method = "euclidean", seed = 42)

  # NNDR ratio = p5_train / p5_holdout
  expect_equal(result$nndr_ratio, result$p5_train / result$p5_holdout)
})

test_that("nndr privacy_pass is logical", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- nndr(X, Y, method = "euclidean", seed = 42)

  expect_type(result$privacy_pass, "logical")
})

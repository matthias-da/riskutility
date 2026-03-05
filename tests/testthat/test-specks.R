# Tests for specks (SPECKS - Propensity Score via KS Test)

test_that("specks returns correct S3 class structure", {
  set.seed(123)
  n <- 200
  X <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    income = rnorm(n, 50000, 10000)
  )
  Y <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    income = rnorm(n, 50000, 10000)
  )

  result <- specks(X, Y, method = "cart")

  expect_s3_class(result, "specks")
  expect_true("specks" %in% names(result))
  expect_true("pMSE" %in% names(result))
  expect_true("pMSE_ratio" %in% names(result))
  expect_true("ks_pvalue" %in% names(result))
  expect_true("propensity_original" %in% names(result))
  expect_true("propensity_synthetic" %in% names(result))
  expect_true("utility_score" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("n_original" %in% names(result))
  expect_true("n_synthetic" %in% names(result))
})

test_that("specks value is in range [0, 1]", {
  set.seed(123)
  n <- 200
  X <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    income = rnorm(n, 50000, 10000)
  )
  Y <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    income = rnorm(n, 50000, 10000)
  )

  result <- specks(X, Y, method = "cart")

  expect_true(result$specks >= 0 && result$specks <= 1)
  expect_true(result$utility_score >= 0 && result$utility_score <= 1)
})

test_that("specks utility_score equals 1 - specks", {
  set.seed(123)
  n <- 200
  X <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )
  Y <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )

  result <- specks(X, Y, method = "cart")

  expect_equal(result$utility_score, 1 - result$specks)
})

test_that("specks is low for similar data", {
  set.seed(123)
  n <- 300
  X <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )
  # Y drawn from same distribution
  Y <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )

  result <- specks(X, Y, method = "cart")

  # SPECKS should be relatively low for indistinguishable data
  expect_true(result$specks < 0.5)
})

test_that("specks is higher for different data", {
  set.seed(123)
  n <- 300
  X <- data.frame(
    age = sample(18:30, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE, prob = c(0.5, 0.5))
  )
  # Very different distribution
  Y <- data.frame(
    age = sample(60:80, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE, prob = c(0.9, 0.1))
  )

  result <- specks(X, Y, method = "cart")

  # SPECKS should be higher for distinguishable data
  expect_true(result$specks > 0.1)
})

test_that("specks print method works", {
  set.seed(123)
  n <- 200
  X <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )
  Y <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )

  result <- specks(X, Y, method = "cart")

  expect_output(print(result), "SPECKS")
  expect_output(print(result), "KS statistic")
  expect_output(print(result), "pMSE")
})

test_that("specks summary method returns correct class", {
  set.seed(123)
  n <- 200
  X <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )
  Y <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE)
  )

  result <- specks(X, Y, method = "cart")
  s <- summary(result)

  expect_s3_class(s, "summary.specks")
  expect_true(!is.null(s$specks))
  expect_true(!is.null(s$pMSE))
  expect_true(!is.null(s$prop_original_summary))
  expect_true(!is.null(s$prop_synthetic_summary))
  expect_output(print(s), "Summary: SPECKS")
})

test_that("specks works with cart method", {
  skip_if_not_installed("rpart")
  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- specks(X, Y, method = "cart")
  expect_s3_class(result, "specks")
  expect_equal(result$method, "cart")
})

test_that("specks works with logit method", {
  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- specks(X, Y, method = "logit")
  expect_s3_class(result, "specks")
  expect_equal(result$method, "logit")
})

test_that("specks input validation works", {
  X <- data.frame(x1 = 1:10)
  Y <- data.frame(x1 = 1:10)

  expect_error(specks(1:10, Y), "X must be a data frame")
  expect_error(specks(X, 1:10), "Y must be a data frame")
})

test_that("specks errors when no common variables", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(b = 1:10)

  expect_error(specks(X, Y), "No common variables")
})

test_that("specks propensity scores have correct length", {
  set.seed(123)
  n_x <- 100
  n_y <- 150
  X <- data.frame(x1 = rnorm(n_x))
  Y <- data.frame(x1 = rnorm(n_y))

  result <- specks(X, Y, method = "logit")

  expect_equal(length(result$propensity_original), n_x)
  expect_equal(length(result$propensity_synthetic), n_y)
  expect_equal(result$n_original, n_x)
  expect_equal(result$n_synthetic, n_y)
})

test_that("specks plot method works", {
  set.seed(123)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- specks(X, Y, method = "cart")

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 3))
  expect_no_error(plot(result, which = 1:3))
})

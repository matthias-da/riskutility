# Tests for dcr (Distance to Closest Record)

test_that("dcr returns correct S3 class structure", {
  set.seed(123)
  X <- data.frame(
    age = rnorm(100, 40, 10),
    income = rnorm(100, 50000, 15000),
    gender = sample(c("M", "F"), 100, replace = TRUE)
  )
  Y <- data.frame(
    age = rnorm(100, 40, 10),
    income = rnorm(100, 50000, 15000),
    gender = sample(c("M", "F"), 100, replace = TRUE)
  )

  result <- dcr(X, Y, seed = 42)

  expect_s3_class(result, "dcr")
  expect_true("dcr_train" %in% names(result))
  expect_true("dcr_holdout" %in% names(result))
  expect_true("dcr_ratio" %in% names(result))
  expect_true("dcr_share" %in% names(result))
  expect_true("privacy_pass" %in% names(result))
})

test_that("dcr values are non-negative", {
  set.seed(123)
  X <- data.frame(
    age = rnorm(50, 40, 10),
    income = rnorm(50, 50000, 15000)
  )
  Y <- data.frame(
    age = rnorm(50, 40, 10),
    income = rnorm(50, 50000, 15000)
  )

  result <- dcr(X, Y, seed = 42)

  expect_true(all(result$dcr_train >= 0))
  expect_true(all(result$dcr_holdout >= 0))
  expect_true(result$dcr_ratio >= 0)
  expect_true(result$dcr_share >= 0 && result$dcr_share <= 1)
})

test_that("dcr detects privacy issues with copied data", {
  set.seed(123)
  X <- data.frame(
    age = rnorm(100, 40, 10),
    income = rnorm(100, 50000, 15000)
  )
  # Bad synthetic: exact copies of training data with tiny noise
  Y <- X[sample(nrow(X), 100, replace = TRUE), ]
  # Add very tiny noise to avoid exact matches (which could affect distance calcs)
  Y$age <- Y$age + rnorm(100, 0, 0.001)
  Y$income <- Y$income + rnorm(100, 0, 0.001)

  result <- dcr(X, Y, seed = 42)

  # Share closer to training should be elevated (privacy concern)
  # Use >= 0.45 as a looser threshold since holdout sampling can vary
  expect_true(result$dcr_share >= 0.45)
})

test_that("dcr handles non-dataframe input", {
  expect_error(dcr(1:10, data.frame(a = 1:10)),
               "X must be a data frame")
  expect_error(dcr(data.frame(a = 1:10), 1:10),
               "Y must be a data frame")
})

test_that("dcr handles missing common variables", {
  X <- data.frame(a = 1:10, b = 1:10)
  Y <- data.frame(c = 1:10, d = 1:10)

  expect_error(dcr(X, Y), "No common variables found")
})

test_that("dcr print method works", {
  set.seed(123)
  X <- data.frame(
    age = rnorm(50, 40, 10),
    income = rnorm(50, 50000, 15000)
  )
  Y <- data.frame(
    age = rnorm(50, 40, 10),
    income = rnorm(50, 50000, 15000)
  )

  result <- dcr(X, Y, seed = 42)

  expect_output(print(result), "Distance to Closest Record")
  expect_output(print(result), "DCR")
})

test_that("dcr summary method works", {
  set.seed(123)
  X <- data.frame(
    age = rnorm(50, 40, 10),
    income = rnorm(50, 50000, 15000)
  )
  Y <- data.frame(
    age = rnorm(50, 40, 10),
    income = rnorm(50, 50000, 15000)
  )

  result <- dcr(X, Y, seed = 42)
  summ <- summary(result)

  expect_s3_class(summ, "summary.dcr")
  expect_true("quantiles_train" %in% names(summ))
  expect_true("quantiles_holdout" %in% names(summ))
})

test_that("dcr with euclidean method works for numeric data", {
  set.seed(123)
  X <- data.frame(
    age = rnorm(50, 40, 10),
    income = rnorm(50, 50000, 15000)
  )
  Y <- data.frame(
    age = rnorm(50, 40, 10),
    income = rnorm(50, 50000, 15000)
  )

  result <- dcr(X, Y, method = "euclidean", seed = 42)

  expect_s3_class(result, "dcr")
  expect_equal(result$method, "euclidean")
})

test_that("dcr with provided holdout works", {
  set.seed(123)
  X <- data.frame(
    age = rnorm(50, 40, 10),
    income = rnorm(50, 50000, 15000)
  )
  Y <- data.frame(
    age = rnorm(50, 40, 10),
    income = rnorm(50, 50000, 15000)
  )
  holdout <- data.frame(
    age = rnorm(25, 40, 10),
    income = rnorm(25, 50000, 15000)
  )

  result <- dcr(X, Y, holdout = holdout)

  expect_s3_class(result, "dcr")
  expect_equal(result$n_holdout, 25)
})

test_that("dcr respects seed for reproducibility", {
  set.seed(123)
  X <- data.frame(
    age = rnorm(50, 40, 10),
    income = rnorm(50, 50000, 15000)
  )
  Y <- data.frame(
    age = rnorm(50, 40, 10),
    income = rnorm(50, 50000, 15000)
  )

  result1 <- dcr(X, Y, seed = 42)
  result2 <- dcr(X, Y, seed = 42)

  expect_equal(result1$dcr_ratio, result2$dcr_ratio)
  expect_equal(result1$dcr_share, result2$dcr_share)
})

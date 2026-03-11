# Tests for .distance_risk_prepare()

test_that(".distance_risk_prepare splits holdout correctly", {
  set.seed(1)
  X <- data.frame(a = 1:100, b = rnorm(100))
  Y <- data.frame(a = 101:200, b = rnorm(100))
  prep <- riskutility:::.distance_risk_prepare(X, Y, holdout_fraction = 0.3,
                                                seed = 42)

  expect_equal(nrow(prep$train) + nrow(prep$holdout), 100)
  expect_equal(nrow(prep$synthetic), 100)
  expect_true(prep$was_split)
  expect_true(all(prep$vars %in% c("a", "b")))
  # Reproducible with same seed
  prep2 <- riskutility:::.distance_risk_prepare(X, Y, holdout_fraction = 0.3,
                                                 seed = 42)
  expect_identical(prep$train, prep2$train)
})

test_that(".distance_risk_prepare uses explicit holdout", {
  X <- data.frame(a = 1:50, b = rnorm(50))
  Y <- data.frame(a = 51:100, b = rnorm(50))
  H <- data.frame(a = 101:120, b = rnorm(20))
  prep <- riskutility:::.distance_risk_prepare(X, Y, holdout = H)

  expect_equal(nrow(prep$train), 50)  # X unchanged
  expect_equal(nrow(prep$holdout), 20)
  expect_false(prep$was_split)
})

test_that(".distance_risk_prepare intersects vars", {
  X <- data.frame(a = 1:10, b = 1:10, c = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10, d = 1:10)
  prep <- riskutility:::.distance_risk_prepare(X, Y)

  expect_equal(sort(prep$vars), c("a", "b"))
})

test_that(".distance_risk_prepare applies vars filter", {
  X <- data.frame(a = 1:10, b = 1:10, c = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10, c = 1:10)
  prep <- riskutility:::.distance_risk_prepare(X, Y, vars = "a")

  expect_equal(prep$vars, "a")
})

test_that(".distance_risk_prepare removes NAs when na.rm = TRUE", {
  X <- data.frame(a = c(1:9, NA), b = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10)
  prep <- riskutility:::.distance_risk_prepare(X, Y, na.rm = TRUE,
                                                holdout_fraction = 0.5,
                                                seed = 1)

  total <- nrow(prep$train) + nrow(prep$holdout)
  expect_true(total <= 9)  # NA row removed before split
})

test_that(".distance_risk_prepare keeps NAs when na.rm = FALSE", {
  X <- data.frame(a = c(1:9, NA), b = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10)
  prep <- riskutility:::.distance_risk_prepare(X, Y, na.rm = FALSE,
                                                holdout_fraction = 0.5,
                                                seed = 1)

  total <- nrow(prep$train) + nrow(prep$holdout)
  expect_equal(total, 10)
})

test_that(".distance_risk_prepare validates holdout_fraction", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = 1:10)
  expect_error(riskutility:::.distance_risk_prepare(X, Y,
                holdout_fraction = 0), "holdout_fraction")
  expect_error(riskutility:::.distance_risk_prepare(X, Y,
                holdout_fraction = 1), "holdout_fraction")
})

test_that(".distance_risk_prepare enforces min_holdout", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = 1:10)
  # With min_holdout = 2 and very small fraction, should still get >= 2
  prep <- riskutility:::.distance_risk_prepare(X, Y,
            holdout_fraction = 0.05, min_holdout = 2, seed = 1)
  expect_true(nrow(prep$holdout) >= 2)
})

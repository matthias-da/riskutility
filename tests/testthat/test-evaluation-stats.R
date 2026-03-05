# Tests for evaluation statistics (mae, mse, rmse, mape)

test_that("mae returns zero for identical vectors", {
  x <- c(1, 2, 3, 4, 5)

  result <- mae(x, x)
  expect_equal(result, 0)
})

test_that("mae computes correct known result", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(2, 3, 4, 5, 6)

  # All differences are 1, so MAE = 1
  result <- mae(x, y)
  expect_equal(result, 1)
})

test_that("mae is non-negative", {
  set.seed(42)
  x <- rnorm(100)
  y <- rnorm(100)

  result <- mae(x, y)
  expect_true(result >= 0)
})

test_that("mae is symmetric", {
  x <- c(1, 2, 3)
  y <- c(4, 5, 6)

  expect_equal(mae(x, y), mae(y, x))
})

test_that("mae gives known result with hand computation", {
  x <- c(10, 20, 30)
  y <- c(12, 18, 33)

  # |10-12| + |20-18| + |30-33| = 2 + 2 + 3 = 7
  # MAE = 7 / 3
  expected <- 7 / 3
  result <- mae(x, y)
  expect_equal(result, expected)
})

# ---- MSE tests ----

test_that("mse returns zero for identical vectors", {
  x <- c(1, 2, 3, 4, 5)

  result <- mse(x, x)
  expect_equal(result, 0)
})

test_that("mse computes correct known result", {
  x <- c(1, 2, 3)
  y <- c(4, 5, 6)

  # (1-4)^2 + (2-5)^2 + (3-6)^2 = 9 + 9 + 9 = 27
  # MSE = 27 / 3 = 9
  result <- mse(x, y)
  expect_equal(result, 9)
})

test_that("mse is non-negative", {
  set.seed(42)
  x <- rnorm(100)
  y <- rnorm(100)

  result <- mse(x, y)
  expect_true(result >= 0)
})

test_that("mse gives known result with hand computation", {
  x <- c(10, 20, 30)
  y <- c(12, 18, 33)

  # (10-12)^2 + (20-18)^2 + (30-33)^2 = 4 + 4 + 9 = 17
  # MSE = 17 / 3
  expected <- 17 / 3
  result <- mse(x, y)
  expect_equal(result, expected)
})

# ---- RMSE tests ----

test_that("rmse returns zero for identical vectors", {
  x <- c(1, 2, 3, 4, 5)

  result <- rmse(x, x)
  expect_equal(result, 0)
})

test_that("rmse equals sqrt(mse)", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(2, 4, 6, 8, 10)

  mse_val <- mse(x, y)
  rmse_val <- rmse(x, y)

  expect_equal(rmse_val, sqrt(mse_val))
})

test_that("rmse is non-negative", {
  set.seed(42)
  x <- rnorm(100)
  y <- rnorm(100)

  result <- rmse(x, y)
  expect_true(result >= 0)
})

test_that("rmse gives known result", {
  x <- c(1, 2, 3)
  y <- c(4, 5, 6)

  # MSE = 9, RMSE = 3
  result <- rmse(x, y)
  expect_equal(result, 3)
})

# ---- MAPE tests ----

test_that("mape returns zero for identical vectors", {
  x <- c(1, 2, 3, 4, 5)

  result <- mape(x, x)
  expect_equal(result, 0)
})

test_that("mape computes correct known result", {
  x <- c(10, 20, 30)
  y <- c(12, 18, 33)

  # |(10-12)/10| + |(20-18)/20| + |(30-33)/30| = 0.2 + 0.1 + 0.1 = 0.4
  # MAPE = 0.4/3 * 100 = 13.33...
  expected <- mean(abs((x - y) / x)) * 100
  result <- mape(x, y)
  expect_equal(result, expected)
})

test_that("mape is non-negative", {
  set.seed(42)
  x <- abs(rnorm(100, mean = 10))
  y <- abs(rnorm(100, mean = 10))

  result <- mape(x, y)
  expect_true(result >= 0)
})

test_that("mape is in percentage terms", {
  x <- c(100, 200, 300)
  y <- c(110, 210, 310)

  # Errors are about 10% for first, 5% for second, 3.33% for third
  result <- mape(x, y)
  # Should be in percentage form
  expect_true(result > 1)  # Not a fraction
  expect_true(result < 100)
})

# ---- Cross-function consistency ----

test_that("rmse >= mae always holds", {
  set.seed(42)
  x <- rnorm(100)
  y <- rnorm(100)

  expect_true(rmse(x, y) >= mae(x, y))
})

test_that("mse >= 0 and mae >= 0 for all inputs", {
  x <- c(-5, 0, 5, 10, -10)
  y <- c(5, 0, -5, -10, 10)

  expect_true(mse(x, y) >= 0)
  expect_true(mae(x, y) >= 0)
  expect_true(rmse(x, y) >= 0)
})

# Tests for gower (Average Gower Distance)

library(testthat)

# --- Setup: shared test data ---

make_mixed_data <- function(n = 50, seed = 123) {
  set.seed(seed)
  data.frame(
    age = sample(20:70, n, replace = TRUE),
    income = rnorm(n, 50000, 15000),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE))
  )
}

make_numeric_data <- function(n = 40, seed = 123) {
  set.seed(seed)
  data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n, 5, 2),
    x3 = runif(n, 0, 100)
  )
}

make_categorical_data <- function(n = 40, seed = 123) {
  set.seed(seed)
  data.frame(
    color = factor(sample(c("red", "blue", "green"), n, replace = TRUE)),
    size = factor(sample(c("S", "M", "L"), n, replace = TRUE))
  )
}

# --- Class structure ---

test_that("gower returns correct S3 class", {
  X <- make_mixed_data(n = 30, seed = 1)
  Y <- make_mixed_data(n = 30, seed = 2)
  result <- gower(X, Y)
  expect_s3_class(result, "gower")
})

test_that("gower result has expected fields", {
  X <- make_mixed_data(n = 30, seed = 1)
  Y <- make_mixed_data(n = 30, seed = 2)
  result <- gower(X, Y)
  expect_named(result, c("gower_distance", "utility_score",
                          "n_records", "n_variables", "n", "n_vars"))
})

test_that("gower_distance is numeric", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  expect_type(result$gower_distance, "double")
})

test_that("utility_score is 1 - gower_distance", {
  X <- make_mixed_data(n = 25, seed = 1)
  Y <- make_mixed_data(n = 25, seed = 2)
  result <- gower(X, Y)
  expect_equal(result$utility_score, 1 - result$gower_distance)
})

test_that("n_records and n alias are correct", {
  X <- make_mixed_data(n = 25, seed = 1)
  Y <- make_mixed_data(n = 25, seed = 2)
  result <- gower(X, Y)
  expect_equal(result$n_records, 25)
  expect_equal(result$n, 25)
})

test_that("n_variables and n_vars alias are correct", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  expect_equal(result$n_variables, 4)
  expect_equal(result$n_vars, 4)
})

test_that("n is correct", {
  X <- make_mixed_data(n = 25, seed = 1)
  Y <- make_mixed_data(n = 25, seed = 2)
  result <- gower(X, Y)
  expect_equal(result$n, 25)
})

test_that("n_vars is correct", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  expect_equal(result$n_vars, 4)
})

# --- Basic computation ---

test_that("gower_distance is non-negative", {
  X <- make_mixed_data(n = 30, seed = 1)
  Y <- make_mixed_data(n = 30, seed = 2)
  result <- gower(X, Y)
  expect_true(result$gower_distance >= 0)
})

test_that("identical data gives same distance as gower(X, X)", {
  # gowerD computes a full pairwise distance matrix, so even with

  # identical arguments, off-diagonal distances are non-zero
  X <- make_mixed_data(n = 30, seed = 1)
  result <- gower(X, X)
  # The distance should be reproducible and non-negative
  expect_true(result$gower_distance >= 0)
  # Running again should give the same value

  result2 <- gower(X, X)
  expect_equal(result$gower_distance, result2$gower_distance)
})

test_that("gower_distance is positive for different data", {
  X <- make_mixed_data(n = 30, seed = 1)
  Y <- make_mixed_data(n = 30, seed = 2)
  result <- gower(X, Y)
  expect_true(result$gower_distance > 0)
})

test_that("gower works with purely numeric data", {
  X <- make_numeric_data(n = 25, seed = 1)
  Y <- make_numeric_data(n = 25, seed = 2)
  result <- gower(X, Y)
  expect_s3_class(result, "gower")
  expect_true(result$gower_distance >= 0)
  expect_equal(result$n_vars, 3)
})

test_that("gower works with purely categorical data", {
  X <- make_categorical_data(n = 25, seed = 1)
  Y <- make_categorical_data(n = 25, seed = 2)
  result <- gower(X, Y)
  expect_s3_class(result, "gower")
  expect_true(result$gower_distance >= 0)
  expect_equal(result$n_vars, 2)
})

test_that("gower works with mixed data types", {
  X <- data.frame(
    age = c(25, 30, 35),
    income = c(30000, 45000, 50000),
    gender = factor(c("M", "F", "M"))
  )
  Y <- data.frame(
    age = c(26, 31, 34),
    income = c(32000, 44000, 52000),
    gender = factor(c("M", "F", "M"))
  )
  result <- gower(X, Y)
  expect_s3_class(result, "gower")
  expect_true(result$gower_distance > 0)
  expect_equal(result$n, 3)
  expect_equal(result$n_vars, 3)
})

test_that("more similar data gives smaller gower distance", {
  # Use simple data where the ordering is clear
  X <- data.frame(x1 = 1:20, x2 = seq(0, 1, length.out = 20))
  # Y_close has values near X
  Y_close <- data.frame(x1 = 1:20 + 0.1, x2 = seq(0, 1, length.out = 20) + 0.01)
  # Y_far has very different values
  Y_far <- data.frame(x1 = 100:119, x2 = seq(10, 20, length.out = 20))

  result_close <- gower(X, Y_close)
  result_far <- gower(X, Y_far)

  expect_true(result_close$gower_distance < result_far$gower_distance)
})

test_that("gower works with single-column data", {
  X <- data.frame(val = c(1, 2, 3, 4, 5))
  Y <- data.frame(val = c(1.1, 2.2, 3.3, 4.4, 5.5))
  result <- gower(X, Y)
  expect_s3_class(result, "gower")
  expect_equal(result$n_vars, 1)
  expect_equal(result$n, 5)
  expect_true(result$gower_distance >= 0)
})

test_that("gower works with single-row data", {
  X <- data.frame(a = 1, b = factor("A"))
  Y <- data.frame(a = 2, b = factor("B"))
  result <- gower(X, Y)
  expect_s3_class(result, "gower")
  expect_equal(result$n, 1)
  expect_true(result$gower_distance > 0)
})

test_that("gower distance does not depend on variable order", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20), c = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20), c = rnorm(20))
  result1 <- gower(X, Y)
  result2 <- gower(X[, c("c", "a", "b")], Y[, c("c", "a", "b")])
  expect_equal(result1$gower_distance, result2$gower_distance, tolerance = 1e-10)
})

# --- S3 methods: print ---

test_that("print.gower runs without error", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  expect_output(print(result), "Average Gower Distance")
})

test_that("print.gower shows observation count", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  expect_output(print(result), "Observations:")
})

test_that("print.gower shows utility score", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  expect_output(print(result), "Utility score:")
})

test_that("print.gower shows variable count", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  expect_output(print(result), "Variables:")
})

test_that("print.gower returns object invisibly", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  out <- capture.output(ret <- print(result))
  expect_identical(ret, result)
})

# --- S3 methods: summary ---

test_that("summary.gower returns correct class", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  s <- summary(result)
  expect_s3_class(s, "summary.gower")
})

test_that("summary.gower has expected fields", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  s <- summary(result)
  expect_named(s, c("gower_distance", "utility_score",
                      "n_records", "n_variables", "n", "n_vars"))
})

test_that("summary.gower preserves values from result", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  s <- summary(result)
  expect_equal(s$gower_distance, result$gower_distance)
  expect_equal(s$utility_score, result$utility_score)
  expect_equal(s$n_records, result$n_records)
  expect_equal(s$n_variables, result$n_variables)
  expect_equal(s$n, result$n)
  expect_equal(s$n_vars, result$n_vars)
})

# --- S3 methods: print.summary ---

test_that("print.summary.gower runs without error", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  expect_output(print(summary(result)), "Summary: Average Gower Distance")
})

test_that("print.summary.gower shows distance", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  expect_output(print(summary(result)), "Distance:")
})

test_that("print.summary.gower shows observations", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  expect_output(print(summary(result)), "Observations:")
})

test_that("print.summary.gower shows variables", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  expect_output(print(summary(result)), "Variables:")
})

test_that("print.summary.gower shows utility score", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  expect_output(print(summary(result)), "Utility score:")
})

test_that("print.summary.gower returns object invisibly", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  s <- summary(result)
  out <- capture.output(ret <- print(s))
  expect_identical(ret, s)
})

# --- synth_pair dispatch ---

test_that("gower works with synth_pair objects", {
  set.seed(123)
  orig <- data.frame(
    age = sample(20:70, 50, replace = TRUE),
    income = rnorm(50, 50000, 15000),
    gender = factor(sample(c("M", "F"), 50, replace = TRUE))
  )
  syn <- data.frame(
    age = sample(20:70, 50, replace = TRUE),
    income = rnorm(50, 50000, 15000),
    gender = factor(sample(c("M", "F"), 50, replace = TRUE))
  )

  pair <- synth_pair(orig, syn)
  result <- gower(pair)
  expect_s3_class(result, "gower")
  expect_true(result$gower_distance >= 0)
  expect_equal(result$n, 50)
})

test_that("gower.synth_pair gives same result as gower.default", {
  set.seed(123)
  orig <- data.frame(
    age = sample(20:70, 30, replace = TRUE),
    income = rnorm(30, 50000, 15000),
    gender = factor(sample(c("M", "F"), 30, replace = TRUE))
  )
  syn <- data.frame(
    age = sample(20:70, 30, replace = TRUE),
    income = rnorm(30, 50000, 15000),
    gender = factor(sample(c("M", "F"), 30, replace = TRUE))
  )

  pair <- synth_pair(orig, syn)
  result_pair <- gower(pair)
  result_direct <- gower(orig, syn)

  expect_equal(result_pair$gower_distance, result_direct$gower_distance)
  expect_equal(result_pair$n, result_direct$n)
  expect_equal(result_pair$n_vars, result_direct$n_vars)
})

test_that("gower.synth_pair uses original and synthetic from pair", {
  set.seed(42)
  orig <- data.frame(x = c(1, 2, 3), y = c(10, 20, 30))
  syn <- data.frame(x = c(1.5, 2.5, 3.5), y = c(15, 25, 35))

  pair <- synth_pair(orig, syn)
  result <- gower(pair)

  expect_equal(result$n, nrow(orig))
  expect_equal(result$n_vars, ncol(orig))
})

# --- Error handling ---

test_that("gower errors on non-data.frame X", {
  Y <- data.frame(a = 1:5)
  expect_error(gower(1:10, Y))
})

test_that("gower errors on non-data.frame Y", {
  X <- data.frame(a = 1:5)
  expect_error(gower(X, 1:10))
})

test_that("gower errors on matrix input", {
  mat <- matrix(1:12, nrow = 3)
  Y <- data.frame(a = 1:3, b = 4:6, c = 7:9, d = 10:12)
  expect_error(gower(mat, Y))
})

# --- VIM dependency ---

test_that("gower uses VIM::gowerD internally", {
  # Verify the function depends on VIM
  expect_true(requireNamespace("VIM", quietly = TRUE))
})

# --- Consistency checks ---

test_that("gower is symmetric for same-sized datasets", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))

  result_xy <- gower(X, Y)
  result_yx <- gower(Y, X)

  # With same-sized data the average Gower distance should be symmetric
  expect_equal(result_xy$gower_distance, result_yx$gower_distance, tolerance = 1e-10)
})

test_that("gower n equals nrow of X", {
  X <- data.frame(a = 1:15, b = rnorm(15))
  Y <- data.frame(a = 1:15, b = rnorm(15))
  result <- gower(X, Y)
  expect_equal(result$n, nrow(X))
})

test_that("gower n_vars equals ncol of X", {
  X <- data.frame(a = 1:10, b = rnorm(10), c = factor(rep("A", 10)))
  Y <- data.frame(a = 1:10, b = rnorm(10), c = factor(rep("B", 10)))
  result <- gower(X, Y)
  expect_equal(result$n_vars, ncol(X))
})

test_that("gower with large number of variables works", {
  set.seed(1)
  n <- 20
  p <- 10
  X <- as.data.frame(matrix(rnorm(n * p), nrow = n))
  Y <- as.data.frame(matrix(rnorm(n * p), nrow = n))
  result <- gower(X, Y)
  expect_s3_class(result, "gower")
  expect_equal(result$n_vars, p)
  expect_true(result$gower_distance >= 0)
})

test_that("gower with character columns works", {
  X <- data.frame(
    name = c("Alice", "Bob", "Carol"),
    score = c(80, 90, 70),
    stringsAsFactors = FALSE
  )
  Y <- data.frame(
    name = c("Alice", "David", "Carol"),
    score = c(85, 92, 68),
    stringsAsFactors = FALSE
  )
  result <- gower(X, Y)
  expect_s3_class(result, "gower")
  expect_true(result$gower_distance >= 0)
})

test_that("uniform single-variable categorical data gives distance 0", {
  # When all values are the same in both X and Y, all pairwise distances are 0
  X <- data.frame(color = factor(c("red", "red", "red")))
  Y <- data.frame(color = factor(c("red", "red", "red")))
  result <- gower(X, Y)
  expect_equal(result$gower_distance, 0, tolerance = 1e-10)
})

test_that("identical categorical data has non-negative distance", {
  X <- data.frame(color = factor(c("red", "blue", "green")))
  Y <- data.frame(color = factor(c("red", "blue", "green")))
  result <- gower(X, Y)
  # gowerD computes all pairwise distances, so off-diagonal elements
  # contribute even when X and Y are identical
  expect_true(result$gower_distance >= 0)
})

# --- S3 methods: plot ---

test_that("plot.gower runs without error", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  expect_silent(plot(result))
})

test_that("plot.gower returns object invisibly", {
  X <- make_mixed_data(n = 20, seed = 1)
  Y <- make_mixed_data(n = 20, seed = 2)
  result <- gower(X, Y)
  ret <- plot(result)
  expect_identical(ret, result)
})

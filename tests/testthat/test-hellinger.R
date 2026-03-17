# Tests for hellinger (Hellinger Distance for Categorical Distributions)

library(testthat)

test_that("hellinger returns correct S3 class structure", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    gender = sample(c("Male", "Female"), n, replace = TRUE),
    region = sample(c("North", "South", "East"), n, replace = TRUE),
    education = sample(c("High", "Medium", "Low"), n, replace = TRUE)
  )
  Y <- data.frame(
    gender = sample(c("Male", "Female"), n, replace = TRUE),
    region = sample(c("North", "South", "East"), n, replace = TRUE),
    education = sample(c("High", "Medium", "Low"), n, replace = TRUE)
  )

  result <- hellinger(X, Y)

  expect_s3_class(result, "hellinger")
  expect_true("per_variable" %in% names(result))
  expect_true("hellinger_mean" %in% names(result))
  expect_true("hellinger_max" %in% names(result))
  expect_true("hellinger_min" %in% names(result))
  expect_true("n_vars" %in% names(result))
  expect_true("vars" %in% names(result))
  expect_true("utility_score" %in% names(result))
  expect_true("n_X" %in% names(result))
  expect_true("n_Y" %in% names(result))
})

test_that("hellinger values are bounded in [0, 1]", {
  set.seed(42)
  n <- 100
  X <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S"), n, replace = TRUE)
  )
  Y <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S"), n, replace = TRUE)
  )

  result <- hellinger(X, Y)

  expect_true(result$hellinger_mean >= 0 && result$hellinger_mean <= 1)
  expect_true(result$hellinger_max >= 0 && result$hellinger_max <= 1)
  expect_true(result$hellinger_min >= 0 && result$hellinger_min <= 1)
  expect_true(all(result$per_variable$hellinger >= 0))
  expect_true(all(result$per_variable$hellinger <= 1))
})

test_that("hellinger is zero for identical distributions", {
  n <- 100
  X <- data.frame(
    gender = rep(c("M", "F"), each = 50),
    region = rep(c("N", "S", "E", "W"), each = 25)
  )
  # Y is identical to X
  Y <- X

  result <- hellinger(X, Y)

  expect_equal(result$hellinger_mean, 0)
  expect_equal(result$hellinger_max, 0)
  expect_equal(result$hellinger_min, 0)
  expect_equal(result$utility_score, 1)
})

test_that("hellinger is 1 for completely disjoint distributions", {
  X <- data.frame(
    color = rep("red", 50)
  )
  Y <- data.frame(
    color = rep("blue", 50)
  )

  result <- hellinger(X, Y)

  expect_equal(result$hellinger_mean, 1)
  expect_equal(result$utility_score, 0)
})

test_that("hellinger utility_score equals 1 - hellinger_mean", {
  set.seed(123)
  n <- 80
  X <- data.frame(
    a = sample(c("x", "y", "z"), n, replace = TRUE),
    b = sample(c("p", "q"), n, replace = TRUE)
  )
  Y <- data.frame(
    a = sample(c("x", "y", "z"), n, replace = TRUE),
    b = sample(c("p", "q"), n, replace = TRUE)
  )

  result <- hellinger(X, Y)

  expect_equal(result$utility_score, 1 - result$hellinger_mean)
})

test_that("hellinger is low for similar distributions", {
  set.seed(123)
  n <- 500
  X <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
  )
  Y <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
  )

  result <- hellinger(X, Y)

  # With large samples drawn from the same distribution, Hellinger should be small
  expect_true(result$hellinger_mean < 0.1)
})

test_that("hellinger is higher for different distributions", {
  set.seed(123)
  n <- 200
  X <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE, prob = c(0.5, 0.5)),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
  )
  Y <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE, prob = c(0.95, 0.05)),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE, prob = c(0.8, 0.1, 0.05, 0.05))
  )

  result <- hellinger(X, Y)

  # Hellinger should be notably higher when distributions differ

  expect_true(result$hellinger_mean > 0.05)
})

test_that("hellinger works with a single variable", {
  set.seed(123)
  n <- 80
  X <- data.frame(color = sample(c("red", "blue", "green"), n, replace = TRUE))
  Y <- data.frame(color = sample(c("red", "blue", "green"), n, replace = TRUE))

  result <- hellinger(X, Y)

  expect_s3_class(result, "hellinger")
  expect_equal(result$n_vars, 1)
  expect_equal(result$hellinger_mean, result$hellinger_max)
  expect_equal(result$hellinger_mean, result$hellinger_min)
})

test_that("hellinger auto-detects categorical variables", {
  set.seed(123)
  n <- 80
  X <- data.frame(
    name = sample(c("A", "B"), n, replace = TRUE),
    score = rnorm(n),  # numeric, should be excluded
    grade = factor(sample(c("pass", "fail"), n, replace = TRUE))
  )
  Y <- data.frame(
    name = sample(c("A", "B"), n, replace = TRUE),
    score = rnorm(n),
    grade = factor(sample(c("pass", "fail"), n, replace = TRUE))
  )

  result <- hellinger(X, Y)

  expect_equal(result$n_vars, 2)
  expect_true("name" %in% result$vars)
  expect_true("grade" %in% result$vars)
  expect_false("score" %in% result$vars)
})

test_that("hellinger works with explicit vars parameter", {
  set.seed(123)
  n <- 80
  X <- data.frame(
    a = sample(c("x", "y"), n, replace = TRUE),
    b = sample(c("p", "q"), n, replace = TRUE),
    c = sample(c("1", "2"), n, replace = TRUE)
  )
  Y <- data.frame(
    a = sample(c("x", "y"), n, replace = TRUE),
    b = sample(c("p", "q"), n, replace = TRUE),
    c = sample(c("1", "2"), n, replace = TRUE)
  )

  result <- hellinger(X, Y, vars = c("a", "b"))

  expect_equal(result$n_vars, 2)
  expect_equal(result$vars, c("a", "b"))
})

test_that("hellinger errors when no categorical variables found", {
  X <- data.frame(x = 1:10, y = rnorm(10))
  Y <- data.frame(x = 1:10, y = rnorm(10))

  expect_error(hellinger(X, Y), "No categorical variables found or specified")
})

test_that("hellinger errors when variables missing in X", {
  X <- data.frame(a = c("x", "y"))
  Y <- data.frame(a = c("x", "y"), b = c("p", "q"))

  expect_error(hellinger(X, Y, vars = c("a", "b")), "Variables missing in X")
})

test_that("hellinger errors when variables missing in Y", {
  X <- data.frame(a = c("x", "y"), b = c("p", "q"))
  Y <- data.frame(a = c("x", "y"))

  expect_error(hellinger(X, Y, vars = c("a", "b")), "Variables missing in Y")
})

test_that("hellinger per_variable data.frame is correctly structured", {
  set.seed(123)
  n <- 80
  X <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E"), n, replace = TRUE)
  )
  Y <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E"), n, replace = TRUE)
  )

  result <- hellinger(X, Y)
  pv <- result$per_variable

  expect_true(is.data.frame(pv))
  expect_equal(nrow(pv), 2)
  expect_true("variable" %in% names(pv))
  expect_true("hellinger" %in% names(pv))
  expect_true("n_categories_X" %in% names(pv))
  expect_true("n_categories_Y" %in% names(pv))
  expect_true("n_categories_union" %in% names(pv))
})

test_that("hellinger handles NA values with na.rm = TRUE", {
  X <- data.frame(
    gender = c("M", "F", NA, "M", "F"),
    region = c("N", "S", "N", NA, "S")
  )
  Y <- data.frame(
    gender = c("M", "F", "M", "F", "M"),
    region = c("N", "S", "N", "S", "N")
  )

  result <- hellinger(X, Y, na.rm = TRUE)

  expect_s3_class(result, "hellinger")
  # Should compute without error; NAs are removed per variable
  expect_true(!is.na(result$hellinger_mean))
})

test_that("hellinger handles different number of categories between X and Y", {
  X <- data.frame(color = c("red", "blue", "green", "red", "blue"))
  Y <- data.frame(color = c("red", "blue", "red", "blue", "yellow"))

  result <- hellinger(X, Y)

  expect_s3_class(result, "hellinger")
  # Union of categories should be 4
  expect_equal(result$per_variable$n_categories_union[1], 4)
  expect_true(result$hellinger_mean > 0)
})

test_that("hellinger print method works", {
  set.seed(123)
  n <- 50
  X <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S"), n, replace = TRUE)
  )
  Y <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S"), n, replace = TRUE)
  )

  result <- hellinger(X, Y)

  expect_output(print(result), "Hellinger Distance")
  expect_output(print(result), "Mean Hellinger distance")
  expect_output(print(result), "Utility score")
  expect_output(print(result), "Interpretation")
})

test_that("hellinger summary method returns correct class", {
  set.seed(123)
  n <- 80
  X <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E"), n, replace = TRUE)
  )
  Y <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E"), n, replace = TRUE)
  )

  result <- hellinger(X, Y)
  s <- summary(result)

  expect_s3_class(s, "summary.hellinger")
  expect_true(!is.null(s$hellinger_mean))
  expect_true(!is.null(s$hellinger_sd))
  expect_true(!is.null(s$worst_variable))
  expect_true(!is.null(s$best_variable))
  expect_true(!is.null(s$per_variable))
  expect_output(print(s), "Summary: Hellinger Distance")
})

test_that("hellinger plot method works for all which values", {
  set.seed(123)
  n <- 80
  X <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E"), n, replace = TRUE),
    education = sample(c("High", "Low"), n, replace = TRUE)
  )
  Y <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E"), n, replace = TRUE),
    education = sample(c("High", "Low"), n, replace = TRUE)
  )

  result <- hellinger(X, Y)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

test_that("hellinger works via synth_pair dispatch", {
  set.seed(123)
  n <- 80
  original <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E"), n, replace = TRUE),
    income = rnorm(n, 50000, 10000)
  )
  synthetic <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E"), n, replace = TRUE),
    income = rnorm(n, 50000, 10000)
  )

  pair <- synth_pair(original, synthetic)
  result <- hellinger(pair)

  expect_s3_class(result, "hellinger")
  # Should only use the categorical variables detected in synth_pair
  expect_equal(result$n_vars, 2)
  expect_true("gender" %in% result$vars)
  expect_true("region" %in% result$vars)
})

test_that("hellinger is symmetric: H(X,Y) == H(Y,X)", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    color = sample(c("red", "blue", "green"), n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  )
  Y <- data.frame(
    color = sample(c("red", "blue", "green"), n, replace = TRUE, prob = c(0.3, 0.5, 0.2))
  )

  result_xy <- hellinger(X, Y)
  result_yx <- hellinger(Y, X)

  expect_equal(result_xy$hellinger_mean, result_yx$hellinger_mean, tolerance = 1e-10)
})

test_that("hellinger records correct dataset sizes", {
  X <- data.frame(a = rep("x", 30))
  Y <- data.frame(a = rep("x", 50))

  result <- hellinger(X, Y)

  expect_equal(result$n_X, 30)
  expect_equal(result$n_Y, 50)
})

test_that("hellinger works with factor variables", {
  set.seed(123)
  n <- 80
  X <- data.frame(
    color = factor(sample(c("red", "blue"), n, replace = TRUE)),
    size = factor(sample(c("S", "M", "L"), n, replace = TRUE), ordered = TRUE)
  )
  Y <- data.frame(
    color = factor(sample(c("red", "blue"), n, replace = TRUE)),
    size = factor(sample(c("S", "M", "L"), n, replace = TRUE), ordered = TRUE)
  )

  result <- hellinger(X, Y)

  expect_s3_class(result, "hellinger")
  expect_equal(result$n_vars, 2)
})

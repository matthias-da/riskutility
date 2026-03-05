# Tests for ci_proximity (Confidence Interval Proximity)

test_that("ci_proximity returns correct S3 class structure", {
  set.seed(123)
  X <- data.frame(
    income = rnorm(200, 50000, 10000),
    age = rnorm(200, 40, 10),
    score = rnorm(200, 100, 15)
  )
  Y <- data.frame(
    income = rnorm(200, 50000, 10000),
    age = rnorm(200, 40, 10),
    score = rnorm(200, 100, 15)
  )

  result <- ci_proximity(X, Y)

  expect_s3_class(result, "ci_proximity")
  expect_true("per_variable" %in% names(result))
  expect_true("proximity_mean" %in% names(result))
  expect_true("overlap_mean" %in% names(result))
  expect_true("relative_error_mean" %in% names(result))
  expect_true("n_vars" %in% names(result))
  expect_true("vars" %in% names(result))
  expect_true("conf.level" %in% names(result))
})

test_that("ci_proximity gives high proximity for identical data", {
  set.seed(123)
  X <- data.frame(
    income = rnorm(500, 50000, 10000),
    age = rnorm(500, 40, 10)
  )

  result <- ci_proximity(X, X)

  # Identical data should have perfect proximity
  expect_equal(result$proximity_mean, 1)
  expect_equal(result$overlap_mean, 1)
  expect_equal(result$relative_error_mean, 0)
})

test_that("ci_proximity gives high proximity for similar data", {
  set.seed(123)
  X <- data.frame(
    income = rnorm(500, 50000, 10000),
    age = rnorm(500, 40, 10)
  )
  Y <- data.frame(
    income = rnorm(500, 50000, 10000),
    age = rnorm(500, 40, 10)
  )

  result <- ci_proximity(X, Y)

  # Similar data should have high proximity
  expect_true(result$proximity_mean > 0.8)
  expect_true(result$overlap_mean > 0.5)
})

test_that("ci_proximity gives lower proximity for different data", {
  set.seed(123)
  X <- data.frame(
    income = rnorm(500, 50000, 10000),
    age = rnorm(500, 40, 10)
  )
  Y <- data.frame(
    income = rnorm(500, 80000, 10000),
    age = rnorm(500, 60, 10)
  )

  result_diff <- ci_proximity(X, Y)

  # Now similar data
  Y_similar <- data.frame(
    income = rnorm(500, 50000, 10000),
    age = rnorm(500, 40, 10)
  )

  result_similar <- ci_proximity(X, Y_similar)

  # Different data should have lower proximity than similar data
  expect_true(result_diff$proximity_mean < result_similar$proximity_mean)
})

test_that("ci_proximity per_variable has correct structure", {
  set.seed(123)
  X <- data.frame(
    income = rnorm(200, 50000, 10000),
    age = rnorm(200, 40, 10)
  )
  Y <- data.frame(
    income = rnorm(200, 50000, 10000),
    age = rnorm(200, 40, 10)
  )

  result <- ci_proximity(X, Y)

  expect_true(is.data.frame(result$per_variable))
  expect_equal(nrow(result$per_variable), 2)
  expect_true("variable" %in% names(result$per_variable))
  expect_true("mean_X" %in% names(result$per_variable))
  expect_true("mean_Y" %in% names(result$per_variable))
  expect_true("overlap" %in% names(result$per_variable))
  expect_true("relative_error" %in% names(result$per_variable))
  expect_true("proximity" %in% names(result$per_variable))
})

test_that("ci_proximity print method works", {
  set.seed(123)
  X <- data.frame(income = rnorm(200, 50000, 10000))
  Y <- data.frame(income = rnorm(200, 50000, 10000))

  result <- ci_proximity(X, Y)

  expect_output(print(result), "Confidence Interval Proximity")
  expect_output(print(result), "proximity score")
})

test_that("ci_proximity summary method returns correct class", {
  set.seed(123)
  X <- data.frame(
    income = rnorm(200, 50000, 10000),
    age = rnorm(200, 40, 10)
  )
  Y <- data.frame(
    income = rnorm(200, 50000, 10000),
    age = rnorm(200, 40, 10)
  )

  result <- ci_proximity(X, Y)
  s <- summary(result)

  expect_s3_class(s, "summary.ci_proximity")
  expect_true(!is.null(s$proximity_mean))
  expect_true(!is.null(s$proximity_sd))
  expect_true(!is.null(s$worst_variable))
  expect_true(!is.null(s$best_variable))
  expect_output(print(s), "Summary: Confidence Interval Proximity")
})

test_that("ci_proximity uses specified vars", {
  set.seed(123)
  X <- data.frame(
    income = rnorm(200, 50000, 10000),
    age = rnorm(200, 40, 10),
    score = rnorm(200, 100, 15)
  )
  Y <- data.frame(
    income = rnorm(200, 50000, 10000),
    age = rnorm(200, 40, 10),
    score = rnorm(200, 100, 15)
  )

  result <- ci_proximity(X, Y, vars = c("income", "age"))

  expect_equal(result$n_vars, 2)
  expect_equal(result$vars, c("income", "age"))
  expect_equal(nrow(result$per_variable), 2)
})

test_that("ci_proximity uses specified confidence level", {
  set.seed(123)
  X <- data.frame(income = rnorm(200, 50000, 10000))
  Y <- data.frame(income = rnorm(200, 50000, 10000))

  result_95 <- ci_proximity(X, Y, conf.level = 0.95)
  result_99 <- ci_proximity(X, Y, conf.level = 0.99)

  expect_equal(result_95$conf.level, 0.95)
  expect_equal(result_99$conf.level, 0.99)
})

test_that("ci_proximity auto-detects numeric variables", {
  set.seed(123)
  X <- data.frame(
    income = rnorm(200, 50000, 10000),
    gender = sample(c("M", "F"), 200, replace = TRUE),
    age = rnorm(200, 40, 10)
  )
  Y <- data.frame(
    income = rnorm(200, 50000, 10000),
    gender = sample(c("M", "F"), 200, replace = TRUE),
    age = rnorm(200, 40, 10)
  )

  result <- ci_proximity(X, Y)

  # Should only use numeric variables (income, age)
  expect_equal(result$n_vars, 2)
  expect_true("income" %in% result$vars)
  expect_true("age" %in% result$vars)
  expect_false("gender" %in% result$vars)
})

test_that("ci_proximity errors with no numeric variables", {
  X <- data.frame(a = c("x", "y", "z"))
  Y <- data.frame(a = c("x", "y", "z"))

  expect_error(ci_proximity(X, Y), "No numeric variables")
})

test_that("ci_proximity overlap is in [0, 1]", {
  set.seed(123)
  X <- data.frame(
    income = rnorm(200, 50000, 10000),
    age = rnorm(200, 40, 10)
  )
  Y <- data.frame(
    income = rnorm(200, 50000, 10000),
    age = rnorm(200, 40, 10)
  )

  result <- ci_proximity(X, Y)

  expect_true(all(result$per_variable$overlap >= 0))
  expect_true(all(result$per_variable$overlap <= 1))
})

test_that("ci_proximity proximity is in [0, 1]", {
  set.seed(123)
  X <- data.frame(income = rnorm(200, 50000, 10000))
  Y <- data.frame(income = rnorm(200, 50000, 10000))

  result <- ci_proximity(X, Y)

  expect_true(all(result$per_variable$proximity >= 0))
  expect_true(all(result$per_variable$proximity <= 1))
})

test_that("ci_proximity plot method works", {
  set.seed(123)
  X <- data.frame(
    income = rnorm(200, 50000, 10000),
    age = rnorm(200, 40, 10)
  )
  Y <- data.frame(
    income = rnorm(200, 50000, 10000),
    age = rnorm(200, 40, 10)
  )

  result <- ci_proximity(X, Y)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 3))
  expect_no_error(plot(result, which = 1:3))
})

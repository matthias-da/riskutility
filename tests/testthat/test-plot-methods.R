# Smoke tests for all plot methods

# Shared test data
setup_test_data <- function() {
  set.seed(42)
  n <- 100
  X <- data.frame(
    age = sample(20:60, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    income = rnorm(n, 50000, 10000)
  )
  Y <- data.frame(
    age = sample(20:60, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    income = rnorm(n, 48000, 12000)
  )
  list(X = X, Y = Y)
}

test_that("plot.dcap works for all which values", {
  d <- setup_test_data()
  X <- d$X
  Y <- d$Y
  X$target <- sample(c("low", "medium", "high"), nrow(X), replace = TRUE)
  Y$target <- sample(c("low", "medium", "high"), nrow(Y), replace = TRUE)

  result <- dcap(X, Y, key_vars = c("age", "gender"), target_var = "target")

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

test_that("plot.tcap works for all which values", {
  d <- setup_test_data()
  X <- d$X
  Y <- d$Y
  X$target <- sample(c("low", "medium", "high"), nrow(X), replace = TRUE)
  Y$target <- sample(c("low", "medium", "high"), nrow(Y), replace = TRUE)

  result <- tcap(X, Y, key_vars = c("age", "gender"), target_var = "target")

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 3))
  expect_no_error(plot(result, which = 1:3))
})

test_that("plot.weap works for all which values", {
  set.seed(42)
  Y <- data.frame(
    age = sample(c("young", "middle", "old"), 100, replace = TRUE),
    gender = sample(c("M", "F"), 100, replace = TRUE),
    income = sample(c("low", "medium", "high"), 100, replace = TRUE)
  )

  result <- weap(Y, key_vars = c("age", "gender"), target_var = "income")

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 3))
  expect_no_error(plot(result, which = 1:3))
})

test_that("plot.disco works for all which values", {
  d <- setup_test_data()
  X <- d$X
  Y <- d$Y
  X$target <- sample(c("low", "medium", "high"), nrow(X), replace = TRUE)
  Y$target <- sample(c("low", "medium", "high"), nrow(Y), replace = TRUE)

  result <- disco(X, Y, key_vars = c("age", "gender"), target_var = "target")

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 3))
  expect_no_error(plot(result, which = 1:3))
})

test_that("plot.dcr works for all which values", {
  set.seed(42)
  n <- 80
  X <- data.frame(
    x1 = rnorm(n, 40, 10),
    x2 = rnorm(n, 50000, 15000)
  )
  Y <- data.frame(
    x1 = rnorm(60, 40, 10),
    x2 = rnorm(60, 50000, 15000)
  )

  result <- dcr(X, Y, method = "euclidean", seed = 42, null_test = FALSE)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 3))
  expect_no_error(plot(result, which = 1:3))
})

test_that("plot.nndr works for all which values", {
  set.seed(42)
  n <- 80
  X <- data.frame(
    x1 = rnorm(n, 40, 10),
    x2 = rnorm(n, 50000, 15000)
  )
  Y <- data.frame(
    x1 = rnorm(60, 40, 10),
    x2 = rnorm(60, 50000, 15000)
  )

  result <- nndr(X, Y, method = "euclidean", seed = 42)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 3))
  expect_no_error(plot(result, which = 1:3))
})

test_that("plot.ims works for all which values", {
  d <- setup_test_data()
  X <- d$X
  Y <- d$Y

  result <- ims(X, Y)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

test_that("plot.rapid works for which = 1:3", {
  skip_if_not_installed("rpart")
  set.seed(42)
  n <- 100
  X <- data.frame(
    age = rnorm(n, 40, 10),
    gender = sample(c("M", "F"), n, replace = TRUE),
    income = rnorm(n, 50000, 10000)
  )
  Y <- data.frame(
    age = rnorm(n, 40, 10),
    gender = sample(c("M", "F"), n, replace = TRUE),
    income = rnorm(n, 48000, 12000)
  )

  result <- rapid(X, Y,
                  key_vars = c("age", "gender"),
                  target_var = "income",
                  model_type = "lm",
                  return_all_records = TRUE)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 3))
})

test_that("plot.hellinger works for all which values", {
  set.seed(42)
  X <- data.frame(
    gender = sample(c("Male", "Female"), 100, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), 100, replace = TRUE)
  )
  Y <- data.frame(
    gender = sample(c("Male", "Female"), 100, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), 100, replace = TRUE)
  )

  result <- hellinger(X, Y)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

test_that("plot.energy_distance works for all which values", {
  set.seed(42)
  X <- data.frame(
    income = rnorm(100, 50000, 10000),
    age = rnorm(100, 40, 10)
  )
  Y <- data.frame(
    income = rnorm(100, 50000, 10000),
    age = rnorm(100, 40, 10)
  )

  result <- energy_distance(X, Y, seed = 42)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

test_that("plot.propscore works for all which values", {
  d <- setup_test_data()
  X <- d$X
  Y <- d$Y

  result <- propscore(X, Y, form = ~ age + gender + income)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

test_that("plot.kanonymity works for all which values", {
  set.seed(42)
  data <- data.frame(
    age = sample(c("young", "middle", "old"), 100, replace = TRUE),
    gender = sample(c("M", "F"), 100, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), 100, replace = TRUE)
  )

  result <- kanonymity(data, key_vars = c("age", "gender", "region"), k = 3)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

test_that("plot.ldiversity works for all which values", {
  set.seed(42)
  data <- data.frame(
    age = sample(c("young", "middle", "old"), 200, replace = TRUE),
    gender = sample(c("M", "F"), 200, replace = TRUE),
    disease = sample(c("healthy", "cold", "flu", "covid"), 200, replace = TRUE)
  )

  result <- ldiversity(data,
                       key_vars = c("age", "gender"),
                       sensitive_var = "disease",
                       l = 2)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

test_that("plot.suda works for all which values", {
  set.seed(42)
  data <- data.frame(
    age = sample(c("young", "middle", "old"), 100, replace = TRUE),
    gender = sample(c("M", "F"), 100, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), 100, replace = TRUE),
    education = sample(c("low", "medium", "high"), 100, replace = TRUE)
  )

  result <- suda(data, key_vars = c("age", "gender", "region", "education"))

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

test_that("plot.individual_risk works for all which values", {
  set.seed(42)
  data <- data.frame(
    age = sample(c("young", "middle", "old"), 200, replace = TRUE),
    gender = sample(c("M", "F"), 200, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), 200, replace = TRUE)
  )

  result <- individual_risk(data,
                            key_vars = c("age", "gender", "region"),
                            method = "sample",
                            threshold = 0.1)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

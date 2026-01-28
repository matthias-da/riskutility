# Tests for RAPID (Risk of Attribute Prediction-Induced Disclosure)

test_that("rapid() works with numeric target and lm model", {
  set.seed(123)
  n <- 100

  original <- data.frame(
    age = sample(20:60, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    income = round(rnorm(n, 50000, 10000))
  )

  synthetic <- data.frame(
    age = original$age + sample(-3:3, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    income = round(original$income * runif(n, 0.9, 1.1))
  )

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = c("age", "gender"),
    target_var = "income",
    model_type = "lm",
    num_epsilon = 10,
    verbose = FALSE
  )

  expect_s3_class(result, "rapid")
  expect_true(result$rapid >= 0 && result$rapid <= 1)
  expect_equal(result$model_type, "lm")
  expect_true(!result$is_categorical)
  expect_equal(result$n_total, n)
  expect_true(result$n_at_risk <= result$n_total)
  expect_true(!is.null(result$model_metrics$mae))
  expect_true(!is.null(result$model_metrics$rmse))
})

test_that("rapid() works with categorical target and rf model", {
  skip_if_not_installed("ranger")

  set.seed(123)
  n <- 100

  original <- data.frame(
    age = sample(20:60, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    health = factor(sample(c("Good", "Fair", "Poor"), n,
                            replace = TRUE, prob = c(0.6, 0.3, 0.1)))
  )

  synthetic <- data.frame(
    age = original$age + sample(-3:3, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    health = factor(sample(c("Good", "Fair", "Poor"), n,
                            replace = TRUE, prob = c(0.55, 0.35, 0.1)))
  )

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = c("age", "gender"),
    target_var = "health",
    model_type = "rf",
    cat_tau = 0.3,
    cat_eval_method = "RCS_marginal",
    verbose = FALSE
  )

  expect_s3_class(result, "rapid")
  expect_true(result$rapid >= 0 && result$rapid <= 1)
  expect_equal(result$model_type, "rf")
  expect_true(result$is_categorical)
  expect_true(!is.null(result$model_metrics$accuracy))
})

test_that("rapid() works with RCS_conditional method", {
  skip_if_not_installed("ranger")

  set.seed(123)
  n <- 100

  original <- data.frame(
    age = sample(20:60, n, replace = TRUE),
    status = factor(sample(c("A", "B"), n, replace = TRUE))
  )

  synthetic <- data.frame(
    age = original$age + sample(-2:2, n, replace = TRUE),
    status = factor(sample(c("A", "B"), n, replace = TRUE))
  )

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "age",
    target_var = "status",
    model_type = "rf",
    cat_tau = 1.25,
    cat_eval_method = "RCS_conditional",
    verbose = FALSE
  )

  expect_s3_class(result, "rapid")
  expect_equal(result$method, "RCS_conditional")
})

test_that("rapid() works with NCE method", {
  skip_if_not_installed("ranger")

  set.seed(123)
  n <- 100

  original <- data.frame(
    x = rnorm(n),
    y = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )

  synthetic <- data.frame(
    x = rnorm(n),
    y = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "x",
    target_var = "y",
    model_type = "rf",
    cat_tau = 0.5,
    cat_eval_method = "NCE",
    verbose = FALSE
  )

  expect_s3_class(result, "rapid")
  expect_equal(result$method, "NCE")
})

test_that("rapid() works with cart model", {
  skip_if_not_installed("rpart")

  set.seed(123)
  n <- 100

  original <- data.frame(
    age = sample(20:60, n, replace = TRUE),
    income = round(rnorm(n, 50000, 10000))
  )

  synthetic <- data.frame(
    age = original$age + sample(-3:3, n, replace = TRUE),
    income = round(original$income * runif(n, 0.9, 1.1))
  )

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "age",
    target_var = "income",
    model_type = "cart",
    num_epsilon = 15,
    verbose = FALSE
  )

  expect_s3_class(result, "rapid")
  expect_equal(result$model_type, "cart")
})

test_that("rapid() returns all records when requested", {
  set.seed(123)
  n <- 50

  original <- data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )

  synthetic <- data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "x",
    target_var = "y",
    model_type = "lm",
    num_epsilon = 20,
    return_all_records = TRUE,
    verbose = FALSE
  )

  expect_equal(nrow(result$records), n)
  expect_true("at_risk" %in% names(result$records))
})

test_that("rapid() print method works", {
  set.seed(123)
  n <- 50

  original <- data.frame(x = rnorm(n), y = rnorm(n))
  synthetic <- data.frame(x = rnorm(n), y = rnorm(n))

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "x",
    target_var = "y",
    model_type = "lm",
    verbose = FALSE
  )

  expect_output(print(result), "RAPID Disclosure Risk Assessment")
  expect_output(print(result), "Model: lm")
})

test_that("rapid() summary method works", {
  set.seed(123)
  n <- 50

  original <- data.frame(x = rnorm(n), y = rnorm(n))
  synthetic <- data.frame(x = rnorm(n), y = rnorm(n))

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "x",
    target_var = "y",
    model_type = "lm",
    verbose = FALSE
  )

  summ <- summary(result)
  expect_s3_class(summ, "summary.rapid")
  expect_true(!is.null(summ$risk_level))
  expect_output(print(summ), "Summary: RAPID")
})

test_that("rapid() validates inputs correctly", {
  set.seed(123)
  n <- 50

  original <- data.frame(x = rnorm(n), y = rnorm(n))
  synthetic <- data.frame(x = rnorm(n), y = rnorm(n))

  # Missing variable

  expect_error(
    rapid(original, synthetic, key_vars = "nonexistent", target_var = "y"),
    "Variables missing"
  )

  # lm with categorical target
  original$cat <- factor(sample(c("A", "B"), n, replace = TRUE))
  synthetic$cat <- factor(sample(c("A", "B"), n, replace = TRUE))

  expect_error(
    rapid(original, synthetic, key_vars = "x", target_var = "cat", model_type = "lm"),
    "not suitable for categorical"
  )

  # logit with more than 2 levels
  original$cat3 <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  synthetic$cat3 <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  expect_error(
    rapid(original, synthetic, key_vars = "x", target_var = "cat3", model_type = "logit"),
    "only supports binary"
  )
})

test_that("rapid() works with synth_pair", {
  set.seed(123)
  n <- 50

  original <- data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )
  synthetic <- data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )

  pair <- synth_pair(original, synthetic,
                     key_vars = "x",
                     target_var = "y")

  result <- rapid(pair, model_type = "lm", verbose = FALSE)
  expect_s3_class(result, "rapid")
})

test_that("rapid() handles absolute epsilon type", {
  set.seed(123)
  n <- 50

  original <- data.frame(
    x = rnorm(n),
    y = rnorm(n, 100, 10)
  )
  synthetic <- data.frame(
    x = rnorm(n),
    y = rnorm(n, 100, 10)
  )

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "x",
    target_var = "y",
    model_type = "lm",
    num_epsilon = 5,
    num_epsilon_type = "absolute",
    verbose = FALSE
  )

  expect_s3_class(result, "rapid")
  expect_true(grepl("absolute", result$method))
})

test_that("rapid() handles NA values in sensitive attribute", {
  set.seed(123)
  n <- 50

  original <- data.frame(
    x = rnorm(n),
    y = c(rnorm(n - 5), rep(NA, 5))
  )
  synthetic <- data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )

  # With drop strategy
  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "x",
    target_var = "y",
    model_type = "lm",
    na_strategy = "drop",
    verbose = FALSE
  )

  expect_s3_class(result, "rapid")
  expect_equal(result$n_total, n - 5)

  # With constant strategy
  result2 <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "x",
    target_var = "y",
    model_type = "lm",
    na_strategy = "constant",
    na_constant_value = 0,
    verbose = FALSE
  )

  expect_s3_class(result2, "rapid")
  expect_equal(result2$n_total, n)
})

# Tests for tcap (Targeted Correct Attribution Probability)

test_that("tcap returns correct S3 class structure", {
  set.seed(42)
  X <- data.frame(
    age = sample(20:60, 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE),
    disease = sample(c("none", "A", "B"), 50, replace = TRUE)
  )
  Y <- X
  Y$disease <- sample(Y$disease)

  result <- tcap(X, Y, key_vars = c("age", "gender"), target_var = "disease")

  expect_s3_class(result, "tcap")
  expect_true("tcap_scores" %in% names(result))
  expect_true("tcap_mean" %in% names(result))
  expect_true("tcap_max" %in% names(result))
  expect_true("tcap_median" %in% names(result))
  expect_true("tcap_certain" %in% names(result))
  expect_true("n_certain" %in% names(result))
  expect_true("baseline" %in% names(result))
})

test_that("tcap values are in valid range [0, 1]", {
  set.seed(42)
  X <- data.frame(
    age = sample(20:40, 30, replace = TRUE),
    gender = sample(c("M", "F"), 30, replace = TRUE),
    disease = sample(c("A", "B"), 30, replace = TRUE)
  )
  Y <- X
  Y$disease <- sample(Y$disease)

  result <- tcap(X, Y, key_vars = c("age", "gender"), target_var = "disease")

  expect_true(result$tcap_mean >= 0 && result$tcap_mean <= 1)
  expect_true(result$baseline >= 0 && result$baseline <= 1)
  expect_true(all(result$tcap_scores[!is.na(result$tcap_scores)] >= 0))
  expect_true(all(result$tcap_scores[!is.na(result$tcap_scores)] <= 1))
})

test_that("tcap handles missing variables correctly", {
  X <- data.frame(a = 1:5, b = letters[1:5])
  Y <- data.frame(a = 1:5, c = letters[1:5])

  expect_error(tcap(X, Y, key_vars = "a", target_var = "b"),
               "Variables missing in Y")
})

test_that("tcap handles non-dataframe input", {
  expect_error(tcap(1:10, data.frame(a = 1:10), key_vars = "a", target_var = "a"),
               "X must be a data frame")
  expect_error(tcap(data.frame(a = 1:10), 1:10, key_vars = "a", target_var = "a"),
               "Y must be a data frame")
})

test_that("tcap print method works", {
  set.seed(42)
  X <- data.frame(
    age = sample(20:60, 20, replace = TRUE),
    gender = sample(c("M", "F"), 20, replace = TRUE),
    disease = sample(c("A", "B"), 20, replace = TRUE)
  )
  Y <- X
  Y$disease <- sample(Y$disease)

  result <- tcap(X, Y, key_vars = c("age", "gender"), target_var = "disease")

  expect_output(print(result), "Correct Attribution Probability")
  expect_output(print(result), "Mean CAP")
})

test_that("tcap summary returns risk categories", {
  set.seed(42)
  X <- data.frame(
    age = sample(20:60, 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE),
    disease = sample(c("A", "B"), 50, replace = TRUE)
  )
  Y <- X
  Y$disease <- sample(Y$disease)

  result <- tcap(X, Y, key_vars = c("age", "gender"), target_var = "disease")
  summ <- summary(result)

  expect_s3_class(summ, "summary.tcap")
  expect_true("n_high_risk" %in% names(summ))
  expect_true("n_medium_risk" %in% names(summ))
  expect_true("n_low_risk" %in% names(summ))
})

test_that("tcap with gower method works", {
  set.seed(42)
  X <- data.frame(
    age = sample(20:60, 30, replace = TRUE),
    gender = sample(c("M", "F"), 30, replace = TRUE),
    disease = sample(c("A", "B", "C"), 30, replace = TRUE)
  )
  Y <- X
  Y$disease <- sample(Y$disease)

  result <- tcap(X, Y,
                 key_vars = c("age", "gender"),
                 target_var = "disease",
                 method = "gower",
                 gower_threshold = 0.15)

  expect_s3_class(result, "tcap")
  expect_equal(result$method, "gower")
})

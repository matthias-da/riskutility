# Tests for dcap (Correct Attribution Probability)

test_that("dcap returns correct S3 class structure", {
  set.seed(123)
  X <- data.frame(
    age = sample(20:60, 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE),
    income = sample(c("low", "medium", "high"), 50, replace = TRUE)
  )
  Y <- X
  Y$income <- sample(Y$income)

  result <- dcap(X, Y, key_vars = c("age", "gender"), target_var = "income")

  expect_s3_class(result, "dcap")
  expect_named(result, c("cap_scores", "n_matches", "dcap", "dcap_median",
                         "n_matched", "n_unmatched", "n_total", "baseline",
                         "key_vars", "target_var", "method", "gower_threshold"))
})

test_that("dcap values are in valid range [0, 1]", {
  set.seed(123)
  X <- data.frame(
    age = sample(20:40, 30, replace = TRUE),
    gender = sample(c("M", "F"), 30, replace = TRUE),
    income = sample(c("low", "high"), 30, replace = TRUE)
  )
  Y <- X
  Y$income <- sample(Y$income)

  result <- dcap(X, Y, key_vars = c("age", "gender"), target_var = "income")

  expect_true(result$dcap >= 0 && result$dcap <= 1)
  expect_true(result$baseline >= 0 && result$baseline <= 1)
  expect_true(all(result$cap_scores[!is.na(result$cap_scores)] >= 0))
  expect_true(all(result$cap_scores[!is.na(result$cap_scores)] <= 1))
})

test_that("dcap handles missing variables correctly", {
  X <- data.frame(a = 1:5, b = letters[1:5])
  Y <- data.frame(a = 1:5, c = letters[1:5])

  expect_error(dcap(X, Y, key_vars = "a", target_var = "b"),
               "Variables missing in Y")
  expect_error(dcap(Y, X, key_vars = "a", target_var = "c"),
               "Variables missing in Y")
})

test_that("dcap handles non-dataframe input", {
  expect_error(dcap(1:10, data.frame(a = 1:10), key_vars = "a", target_var = "a"),
               "X must be a data frame")
  expect_error(dcap(data.frame(a = 1:10), 1:10, key_vars = "a", target_var = "a"),
               "Y must be a data frame")
})

test_that("dcap returns high CAP for identical data", {
  set.seed(123)
  X <- data.frame(
    age = sample(20:30, 20, replace = TRUE),
    gender = sample(c("M", "F"), 20, replace = TRUE),
    income = sample(c("low", "high"), 20, replace = TRUE)
  )
  Y <- X  # Identical data

  result <- dcap(X, Y, key_vars = c("age", "gender"), target_var = "income")

  # When data is identical, matched records should have CAP = 1
  matched_caps <- result$cap_scores[!is.na(result$cap_scores)]
  expect_true(mean(matched_caps) >= 0.5)
})

test_that("dcap with gower method works", {
  set.seed(123)
  X <- data.frame(
    age = sample(20:60, 30, replace = TRUE),
    gender = sample(c("M", "F"), 30, replace = TRUE),
    income = sample(c("low", "medium", "high"), 30, replace = TRUE)
  )
  Y <- X
  Y$income <- sample(Y$income)

  result <- dcap(X, Y,
                 key_vars = c("age", "gender"),
                 target_var = "income",
                 method = "gower",
                 gower_threshold = 0.2)

  expect_s3_class(result, "dcap")
  expect_equal(result$method, "gower")
  expect_equal(result$gower_threshold, 0.2)
})

test_that("dcap print method works", {
  set.seed(123)
  X <- data.frame(
    age = sample(20:60, 20, replace = TRUE),
    gender = sample(c("M", "F"), 20, replace = TRUE),
    income = sample(c("low", "high"), 20, replace = TRUE)
  )
  Y <- X
  Y$income <- sample(Y$income)

  result <- dcap(X, Y, key_vars = c("age", "gender"), target_var = "income")

  expect_output(print(result), "Correct Attribution Probability")
  expect_output(print(result), "DCAP")
})

test_that("dcap summary method works", {
  set.seed(123)
  X <- data.frame(
    age = sample(20:60, 20, replace = TRUE),
    gender = sample(c("M", "F"), 20, replace = TRUE),
    income = sample(c("low", "high"), 20, replace = TRUE)
  )
  Y <- X
  Y$income <- sample(Y$income)

  result <- dcap(X, Y, key_vars = c("age", "gender"), target_var = "income")
  summ <- summary(result)

  expect_s3_class(summ, "summary.dcap")
  expect_true("risk_ratio" %in% names(summ))
  expect_true("cap_quantiles" %in% names(summ))
})

test_that("dcap handles NA values correctly", {
  X <- data.frame(
    age = c(25, 30, NA, 40, 45),
    gender = c("M", "F", "M", NA, "F"),
    income = c("low", "high", "low", "high", "low")
  )
  Y <- data.frame(
    age = c(25, 30, 35, 40, 45),
    gender = c("M", "F", "M", "F", "F"),
    income = c("low", "high", "medium", "high", "low")
  )

  # With na.rm = TRUE (default), should work
  result <- dcap(X, Y, key_vars = c("age", "gender"), target_var = "income")
  expect_s3_class(result, "dcap")
  expect_equal(result$n_total, 3)  # 2 rows with NA removed from X
})

test_that("dcap handles continuous target with binning", {
  set.seed(123)
  X <- data.frame(
    age = sample(20:60, 30, replace = TRUE),
    gender = sample(c("M", "F"), 30, replace = TRUE),
    salary = rnorm(30, 50000, 10000)
  )
  Y <- X
  Y$salary <- Y$salary + rnorm(30, 0, 5000)

  result <- dcap(X, Y,
                 key_vars = c("age", "gender"),
                 target_var = "salary",
                 cont_bins = 5)

  expect_s3_class(result, "dcap")
})

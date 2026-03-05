# Tests for weap (Within Equivalence Class Attribution Probability)

test_that("weap returns correct S3 class structure", {
  set.seed(123)
  Y <- data.frame(
    age = sample(c("young", "middle", "old"), 200, replace = TRUE),
    gender = sample(c("M", "F"), 200, replace = TRUE),
    region = sample(c("N", "S"), 200, replace = TRUE),
    income = sample(c("low", "medium", "high"), 200, replace = TRUE)
  )

  result <- weap(Y, key_vars = c("age", "gender", "region"), target_var = "income")

  expect_s3_class(result, "weap")
  expect_named(result, c("weap_scores", "weap_mean", "weap_median",
                          "equivalence_classes", "n_disclosive",
                          "pct_disclosive", "n_total", "n_unique_keys",
                          "key_vars", "target_var"))
})

test_that("weap values are in valid range [0, 1]", {
  set.seed(42)
  Y <- data.frame(
    age = sample(c("young", "old"), 100, replace = TRUE),
    gender = sample(c("M", "F"), 100, replace = TRUE),
    income = sample(c("low", "high"), 100, replace = TRUE)
  )

  result <- weap(Y, key_vars = c("age", "gender"), target_var = "income")

  expect_true(result$weap_mean >= 0 && result$weap_mean <= 1)
  expect_true(result$weap_median >= 0 && result$weap_median <= 1)
  expect_true(all(result$weap_scores >= 0))
  expect_true(all(result$weap_scores <= 1))
  expect_true(result$pct_disclosive >= 0 && result$pct_disclosive <= 100)
})

test_that("weap is 1 for all records when key uniquely determines target", {
  # Every record has a unique key combination and same target value
  Y <- data.frame(
    age = c("young", "young", "old", "old"),
    gender = c("M", "F", "M", "F"),
    income = c("low", "low", "low", "low")
  )

  result <- weap(Y, key_vars = c("age", "gender"), target_var = "income")

  # All records have WEAP = 1 because every equivalence class has only one target value

  expect_true(all(result$weap_scores == 1))
  expect_equal(result$n_disclosive, nrow(Y))
  expect_equal(result$pct_disclosive, 100)
})

test_that("weap is high when records have unique key combinations", {
  # Create data where each record has unique key combination
  Y <- data.frame(
    age = as.character(1:50),
    gender = rep(c("M", "F"), 25),
    income = sample(c("low", "medium", "high"), 50, replace = TRUE)
  )

  result <- weap(Y, key_vars = c("age", "gender"), target_var = "income")

  # With unique keys, every WEAP should be 1
  expect_true(all(result$weap_scores == 1))
  expect_equal(result$weap_mean, 1)
  expect_equal(result$n_unique_keys, 50)
})

test_that("weap print method works", {
  set.seed(123)
  Y <- data.frame(
    age = sample(c("young", "old"), 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE),
    income = sample(c("low", "high"), 50, replace = TRUE)
  )

  result <- weap(Y, key_vars = c("age", "gender"), target_var = "income")

  expect_output(print(result), "Within Equivalence Class Attribution Probability")
  expect_output(print(result), "Mean WEAP")
  expect_output(print(result), "Median WEAP")
})

test_that("weap summary method returns correct class", {
  set.seed(123)
  Y <- data.frame(
    age = sample(c("young", "old"), 100, replace = TRUE),
    gender = sample(c("M", "F"), 100, replace = TRUE),
    income = sample(c("low", "medium", "high"), 100, replace = TRUE)
  )

  result <- weap(Y, key_vars = c("age", "gender"), target_var = "income")
  s <- summary(result)

  expect_s3_class(s, "summary.weap")
  expect_true(!is.null(s$weap_mean))
  expect_true(!is.null(s$weap_sd))
  expect_true(!is.null(s$weap_quantiles))
  expect_true(!is.null(s$n_singleton_keys))
  expect_true(!is.null(s$n_disclosive_classes))
  expect_output(print(s), "Summary: Within Equivalence Class Attribution Probability")
})

test_that("weap plot method works for all which values", {
  set.seed(123)
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

test_that("weap handles NA with na.rm = TRUE", {
  Y <- data.frame(
    age = c("young", "old", NA, "young", "old"),
    gender = c("M", "F", "M", NA, "F"),
    income = c("low", "high", "low", "medium", "high")
  )

  result <- weap(Y, key_vars = c("age", "gender"), target_var = "income", na.rm = TRUE)

  # Should only process complete cases (3 records: rows 1, 2, 5)
  expect_equal(result$n_total, 3)
  expect_s3_class(result, "weap")
})

test_that("weap stops when all cases have NA", {
  Y <- data.frame(
    age = c(NA, NA),
    gender = c("M", "F"),
    income = c("low", "high")
  )

  expect_error(
    weap(Y, key_vars = c("age", "gender"), target_var = "income", na.rm = TRUE),
    "No complete cases"
  )
})

test_that("weap errors on non-data.frame input", {
  expect_error(
    weap(matrix(1:10, ncol = 2), key_vars = "V1", target_var = "V2"),
    "Y must be a data frame"
  )
})

test_that("weap errors when variables are missing", {
  Y <- data.frame(a = 1:5, b = letters[1:5])

  expect_error(
    weap(Y, key_vars = c("a", "nonexistent"), target_var = "b"),
    "Variables missing in Y"
  )

  expect_error(
    weap(Y, key_vars = "a", target_var = "nonexistent"),
    "Variables missing in Y"
  )
})

test_that("weap equivalence classes are computed correctly", {
  # 3 distinct key combinations
  Y <- data.frame(
    age = c("young", "young", "old", "old", "young", "old"),
    gender = c("M", "M", "F", "F", "M", "F"),
    income = c("low", "high", "low", "high", "low", "low")
  )

  result <- weap(Y, key_vars = c("age", "gender"), target_var = "income")

  expect_equal(result$n_unique_keys, 2)
  expect_equal(result$n_total, 6)

  # Check equivalence class sizes
  ec <- result$equivalence_classes
  expect_true(all(ec$size > 0))
  expect_equal(sum(ec$size), result$n_total)
})

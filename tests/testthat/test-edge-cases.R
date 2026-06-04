# Tests for edge cases: all-NA, single-record, high-dimensional, etc.

library(testthat)

# ============================================================================
# All-NA Data Edge Cases
# ============================================================================

test_that("dcap handles NA values correctly", {
  set.seed(123)

  X <- data.frame(
    age = c(25, 30, 35, NA, 40),
    gender = c("M", "F", "M", "F", NA),
    income = c("low", "medium", "high", "low", "medium")
  )

  Y <- data.frame(
    age = c(25, 30, 35, 40, 45),
    gender = c("M", "F", "M", "F", "M"),
    income = c("low", "medium", "high", "low", "medium")
  )

  # Should work with na.rm = TRUE (default)
  result <- dcap(X, Y,
                 key_vars = c("age", "gender"),
                 target_var = "income",
                 na.rm = TRUE)

  expect_s3_class(result, "dcap")
  expect_true(result$n_total < 5)  # Some records removed due to NA
})

test_that("functions error on all-NA data", {
  # Create NA data with correct types to match Y
  X_all_na <- data.frame(
    age = as.integer(rep(NA, 10)),
    gender = as.character(rep(NA, 10)),
    income = as.character(rep(NA, 10)),
    stringsAsFactors = FALSE
  )

  Y <- data.frame(
    age = 1:10,
    gender = rep(c("M", "F"), 5),
    income = rep(c("low", "high"), 5),
    stringsAsFactors = FALSE
  )

  expect_error(
    dcap(X_all_na, Y,
         key_vars = c("age", "gender"),
         target_var = "income"),
    "complete cases|NA|no data|empty"
  )
})

test_that("hellinger handles NA values", {
  set.seed(123)

  X <- data.frame(
    gender = c("M", "F", "M", NA, "F"),
    region = c("N", "S", NA, "E", "W")
  )

  Y <- data.frame(
    gender = c("M", "F", "M", "F", "M"),
    region = c("N", "S", "E", "E", "W")
  )

  result <- hellinger(X, Y, na.rm = TRUE)

  expect_s3_class(result, "hellinger")
  expect_true(!is.na(result$hellinger_mean))
})

# ============================================================================
# Single-Record Edge Cases
# ============================================================================

test_that("dcap handles single-record data", {
  X <- data.frame(
    age = 25,
    gender = "M",
    income = "low"
  )

  Y <- data.frame(
    age = 25,
    gender = "M",
    income = "low"
  )

  result <- dcap(X, Y,
                 key_vars = c("age", "gender"),
                 target_var = "income")

  expect_s3_class(result, "dcap")
  expect_equal(result$n_total, 1)
})

test_that("kanonymity handles single-record data", {
  X <- data.frame(
    age = 25,
    gender = "M",
    region = "N"
  )

  result <- kanonymity(X, key_vars = c("age", "gender", "region"))

  expect_s3_class(result, "kanonymity")
  expect_equal(result$k_level, 1)  # Single record is unique
})

test_that("energy_distance handles small samples", {
  set.seed(123)

  X <- data.frame(income = rnorm(5), age = rnorm(5))
  Y <- data.frame(income = rnorm(5), age = rnorm(5))

  result <- energy_distance(X, Y, n_sample = NULL, seed = 42)

  expect_s3_class(result, "energy_distance")
  expect_equal(result$n_X, 5)
})

# ============================================================================
# High-Dimensional Edge Cases
# ============================================================================

test_that("dcap handles many key variables", {
  set.seed(123)
  n <- 100

  # Create data with many key variables
  X <- data.frame(
    var1 = sample(1:5, n, replace = TRUE),
    var2 = sample(1:5, n, replace = TRUE),
    var3 = sample(1:5, n, replace = TRUE),
    var4 = sample(1:5, n, replace = TRUE),
    var5 = sample(1:5, n, replace = TRUE),
    target = sample(c("A", "B", "C"), n, replace = TRUE)
  )

  Y <- data.frame(
    var1 = sample(1:5, n, replace = TRUE),
    var2 = sample(1:5, n, replace = TRUE),
    var3 = sample(1:5, n, replace = TRUE),
    var4 = sample(1:5, n, replace = TRUE),
    var5 = sample(1:5, n, replace = TRUE),
    target = sample(c("A", "B", "C"), n, replace = TRUE)
  )

  result <- dcap(X, Y,
                 key_vars = c("var1", "var2", "var3", "var4", "var5"),
                 target_var = "target")

  expect_s3_class(result, "dcap")
})

test_that("suda handles many variables", {
  set.seed(123)
  n <- 50

  # Create data with many variables (limit max_k to keep it tractable)
  data <- data.frame(
    v1 = sample(LETTERS[1:3], n, replace = TRUE),
    v2 = sample(LETTERS[1:3], n, replace = TRUE),
    v3 = sample(LETTERS[1:3], n, replace = TRUE),
    v4 = sample(LETTERS[1:3], n, replace = TRUE),
    v5 = sample(LETTERS[1:3], n, replace = TRUE)
  )

  result <- suda(data,
                 key_vars = c("v1", "v2", "v3", "v4", "v5"),
                 max_k = 3)  # Limit complexity

  expect_s3_class(result, "suda")
})

# ============================================================================
# Empty or Minimal Data
# ============================================================================

test_that("functions handle empty after NA removal gracefully", {
  X <- data.frame(
    age = c(NA, NA),
    gender = c(NA, NA),
    income = c(NA, NA)
  )

  Y <- data.frame(
    age = c(25, 30),
    gender = c("M", "F"),
    income = c("low", "high")
  )

  expect_error(
    dcap(X, Y,
         key_vars = c("age", "gender"),
         target_var = "income",
         na.rm = TRUE)
  )
})

test_that("compare_ks_test handles identical distributions", {
  set.seed(123)
  data <- rnorm(100)

  X <- data.frame(value = data)
  Y <- data.frame(value = data)  # Identical

  result <- compare_ks_test(X, Y, num_var = "value")

  expect_equal(result$ks_statistic, 0)
  expect_equal(result$p_value, 1)
})

# ============================================================================
# Boundary Values
# ============================================================================

test_that("dcap returns 1 for perfect match", {
  set.seed(123)

  X <- data.frame(
    age = c(25, 30, 35),
    gender = c("M", "F", "M"),
    income = c("low", "medium", "high")
  )

  Y <- X  # Identical data

  result <- dcap(X, Y,
                 key_vars = c("age", "gender"),
                 target_var = "income")

  expect_equal(result$cap, 1)
})

test_that("hellinger returns 0 for identical categorical distributions", {
  X <- data.frame(
    gender = rep(c("M", "F"), each = 50),
    region = rep(c("N", "S", "E", "W"), 25)
  )

  Y <- X  # Identical

  result <- hellinger(X, Y)

  expect_equal(result$hellinger_mean, 0)
})

test_that("energy_distance returns 0 for identical distributions", {
  set.seed(123)
  data <- data.frame(income = rnorm(100), age = rnorm(100))

  result <- energy_distance(data, data, n_sample = NULL, seed = 42)

  expect_equal(result$energy_distance, 0, tolerance = 1e-10)
})

# ============================================================================
# Special Character and Factor Handling
# ============================================================================

test_that("dcap handles factor variables", {
  set.seed(123)

  X <- data.frame(
    age = factor(c("young", "middle", "old")),
    gender = factor(c("M", "F", "M")),
    income = factor(c("low", "medium", "high"))
  )

  Y <- data.frame(
    age = factor(c("young", "middle", "old")),
    gender = factor(c("M", "F", "F")),
    income = factor(c("low", "medium", "low"))
  )

  result <- dcap(X, Y,
                 key_vars = c("age", "gender"),
                 target_var = "income")

  expect_s3_class(result, "dcap")
})

test_that("hellinger handles new categories in synthetic", {
  X <- data.frame(
    color = c("red", "blue", "green")
  )

  Y <- data.frame(
    color = c("red", "blue", "yellow")  # yellow is new
  )

  result <- hellinger(X, Y)

  expect_s3_class(result, "hellinger")
  # Should handle the asymmetric categories
  expect_true(result$hellinger_mean > 0)
})

# ============================================================================
# Large Category Counts
# ============================================================================
test_that("hellinger handles many categories", {
  set.seed(123)

  # Many unique values
  X <- data.frame(
    id = as.character(1:100)
  )

  Y <- data.frame(
    id = as.character(51:150)  # 50% overlap
  )

  result <- hellinger(X, Y)

  expect_s3_class(result, "hellinger")
})

# ============================================================================
# Unbalanced Sample Sizes
# ============================================================================

test_that("dcap handles unbalanced sample sizes", {
  set.seed(123)

  X <- data.frame(
    age = sample(20:60, 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE),
    income = sample(c("low", "medium", "high"), 50, replace = TRUE)
  )

  Y <- data.frame(
    age = sample(20:60, 500, replace = TRUE),  # 10x larger
    gender = sample(c("M", "F"), 500, replace = TRUE),
    income = sample(c("low", "medium", "high"), 500, replace = TRUE)
  )

  result <- dcap(X, Y,
                 key_vars = c("age", "gender"),
                 target_var = "income")

  expect_s3_class(result, "dcap")
  expect_equal(result$n_total, 50)  # Should match X size
})

test_that("energy_distance handles unbalanced samples", {
  set.seed(123)

  X <- data.frame(income = rnorm(50), age = rnorm(50))
  Y <- data.frame(income = rnorm(500), age = rnorm(500))

  result <- energy_distance(X, Y, n_sample = NULL, seed = 42)

  expect_s3_class(result, "energy_distance")
})

# ============================================================================
# Mixed Type Data
# ============================================================================

test_that("synth_pair correctly separates numeric and categorical", {
  set.seed(123)

  data <- data.frame(
    num1 = rnorm(50),
    num2 = runif(50),
    int1 = sample(1:10, 50, replace = TRUE),
    cat1 = sample(c("A", "B"), 50, replace = TRUE),
    cat2 = factor(sample(c("X", "Y", "Z"), 50, replace = TRUE)),
    stringsAsFactors = FALSE
  )

  pair <- synth_pair(data, data)

  expect_true(all(c("num1", "num2", "int1") %in% pair$num_vars))
  expect_true(all(c("cat1", "cat2") %in% pair$cat_vars))
})

# ============================================================================
# Reproducibility
# ============================================================================

test_that("energy_distance is reproducible with seed", {
  set.seed(123)
  X <- data.frame(income = rnorm(100), age = rnorm(100))
  Y <- data.frame(income = rnorm(100), age = rnorm(100))

  result1 <- energy_distance(X, Y, seed = 42)
  result2 <- energy_distance(X, Y, seed = 42)

  expect_equal(result1$energy_distance, result2$energy_distance)
})

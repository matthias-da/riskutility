# Tests for utility functions: propscore, hellinger, energy_distance

library(testthat)

# ============================================================================
# Test Data Setup
# ============================================================================

create_utility_test_data <- function(n = 200, seed = 123) {
  set.seed(seed)

  X <- data.frame(
    income = rnorm(n, mean = 50000, sd = 10000),
    age = rnorm(n, mean = 40, sd = 10),
    score = rnorm(n, mean = 100, sd = 15),
    gender = sample(c("Male", "Female"), n, replace = TRUE),
    region = sample(c("North", "South", "East", "West"), n, replace = TRUE),
    education = sample(c("High", "Medium", "Low"), n, replace = TRUE),
    weight = runif(n, 0.5, 1.5),
    stringsAsFactors = FALSE
  )

  # Good synthetic (similar distribution)
  Y_good <- data.frame(
    income = rnorm(n, mean = 50000, sd = 10000),
    age = rnorm(n, mean = 40, sd = 10),
    score = rnorm(n, mean = 100, sd = 15),
    gender = sample(c("Male", "Female"), n, replace = TRUE),
    region = sample(c("North", "South", "East", "West"), n, replace = TRUE),
    education = sample(c("High", "Medium", "Low"), n, replace = TRUE),
    weight = runif(n, 0.5, 1.5),
    stringsAsFactors = FALSE
  )

  # Poor synthetic (shifted distributions)
  Y_poor <- data.frame(
    income = rnorm(n, mean = 70000, sd = 20000),
    age = rnorm(n, mean = 50, sd = 15),
    score = rnorm(n, mean = 80, sd = 20),
    gender = sample(c("Male", "Female"), n, replace = TRUE, prob = c(0.9, 0.1)),
    region = sample(c("North", "South", "East", "West"), n, replace = TRUE,
                    prob = c(0.7, 0.1, 0.1, 0.1)),
    education = sample(c("High", "Medium", "Low"), n, replace = TRUE,
                       prob = c(0.1, 0.1, 0.8)),
    weight = runif(n, 0.5, 1.5),
    stringsAsFactors = FALSE
  )

  list(X = X, Y_good = Y_good, Y_poor = Y_poor)
}

# ============================================================================
# hellinger
# ============================================================================

test_that("hellinger returns correct S3 class", {
  data <- create_utility_test_data()

  result <- hellinger(data$X, data$Y_good)

  expect_s3_class(result, "hellinger")
})

test_that("hellinger returns bounded values", {
  data <- create_utility_test_data()

  result <- hellinger(data$X, data$Y_good)

  expect_true(result$hellinger_mean >= 0)
  expect_true(result$hellinger_mean <= 1)
  expect_true(result$hellinger_max >= 0)
  expect_true(result$hellinger_max <= 1)
})

test_that("hellinger detects poor vs good synthetic", {
  data <- create_utility_test_data()

  result_good <- hellinger(data$X, data$Y_good)
  result_poor <- hellinger(data$X, data$Y_poor)

  # Poor synthetic should have higher Hellinger distance
  expect_true(result_poor$hellinger_mean > result_good$hellinger_mean)
})

test_that("hellinger auto-detects categorical variables", {
  data <- create_utility_test_data()

  result <- hellinger(data$X, data$Y_good)

  # Should have detected gender, region, education
  expect_true(length(result$vars) >= 3)
})

test_that("hellinger works with specified variables", {
  data <- create_utility_test_data()

  result <- hellinger(data$X, data$Y_good, vars = c("gender", "region"))

  expect_equal(result$n_vars, 2)
  expect_equal(sort(result$vars), sort(c("gender", "region")))
})

test_that("hellinger handles weights", {
  data <- create_utility_test_data()

  result <- hellinger(data$X, data$Y_good,
                      weight_X = "weight",
                      weight_Y = "weight")

  expect_s3_class(result, "hellinger")
  expect_true(!is.na(result$hellinger_mean))
})

test_that("hellinger print method works", {
  data <- create_utility_test_data()
  result <- hellinger(data$X, data$Y_good)

  expect_output(print(result), "Hellinger Distance")
})

test_that("hellinger summary method works", {
  data <- create_utility_test_data()
  result <- hellinger(data$X, data$Y_good)

  summ <- summary(result)

  expect_s3_class(summ, "summary.hellinger")
  expect_output(print(summ), "Summary")
})

test_that("hellinger plot method works", {
  data <- create_utility_test_data()
  result <- hellinger(data$X, data$Y_good)

  expect_no_error(plot(result, which = 1))
})

test_that("hellinger validates input", {
  data <- create_utility_test_data()

  # No categorical variables
  X_num <- data$X[, c("income", "age", "score")]
  Y_num <- data$Y_good[, c("income", "age", "score")]

  expect_error(hellinger(X_num, Y_num), "categorical")

  # Missing variable
  expect_error(hellinger(data$X, data$Y_good, vars = c("nonexistent")))
})

# ============================================================================
# energy_distance
# ============================================================================

test_that("energy_distance returns correct S3 class", {
  data <- create_utility_test_data()

  result <- energy_distance(data$X, data$Y_good, seed = 42)

  expect_s3_class(result, "energy_distance")
})

test_that("energy_distance is non-negative", {
  data <- create_utility_test_data()

  result <- energy_distance(data$X, data$Y_good, seed = 42)

  expect_true(result$energy_distance >= 0)
  expect_true(result$energy_distance_normalized >= 0)
})

test_that("energy_distance detects poor vs good synthetic", {
  data <- create_utility_test_data()

  result_good <- energy_distance(data$X, data$Y_good, seed = 42)
  result_poor <- energy_distance(data$X, data$Y_poor, seed = 42)

  # Poor synthetic should have higher energy distance
  expect_true(result_poor$energy_distance > result_good$energy_distance)
})

test_that("energy_distance auto-detects numeric variables", {
  data <- create_utility_test_data()

  result <- energy_distance(data$X, data$Y_good, seed = 42)

  # Should have detected income, age, score (and maybe weight)
  expect_true(result$n_vars >= 3)
})

test_that("energy_distance works with specified variables", {
  data <- create_utility_test_data()

  result <- energy_distance(data$X, data$Y_good,
                            vars = c("income", "age"),
                            seed = 42)

  expect_equal(result$n_vars, 2)
})

test_that("energy_distance handles standardization", {
  data <- create_utility_test_data()

  result_std <- energy_distance(data$X, data$Y_good,
                                standardize = TRUE, seed = 42)
  result_raw <- energy_distance(data$X, data$Y_good,
                                standardize = FALSE, seed = 42)

  expect_true(result_std$standardized)
  expect_false(result_raw$standardized)
  # Both should be valid non-negative values
  expect_true(result_std$energy_distance >= 0)
  expect_true(result_raw$energy_distance >= 0)
  # With variables on very different scales (income ~50000 vs age ~40),

  # standardization typically produces different results
  # But we mainly verify the flag is set correctly
})

test_that("energy_distance handles sampling for large datasets", {
  data <- create_utility_test_data(n = 2000)

  result <- energy_distance(data$X, data$Y_good,
                            n_sample = 500, seed = 42)

  expect_equal(result$n_X, 500)
  expect_equal(result$n_Y, 500)
  expect_true(result$sampled)
})

test_that("energy_distance print method works", {
  data <- create_utility_test_data()
  result <- energy_distance(data$X, data$Y_good, seed = 42)

  expect_output(print(result), "Energy Distance")
})

test_that("energy_distance summary method works", {
  data <- create_utility_test_data()
  result <- energy_distance(data$X, data$Y_good, seed = 42)

  summ <- summary(result)

  expect_s3_class(summ, "summary.energy_distance")
  expect_output(print(summ), "Summary")
})

test_that("energy_distance validates input", {
  data <- create_utility_test_data()

  # No numeric variables
  X_cat <- data$X[, c("gender", "region", "education")]
  Y_cat <- data$Y_good[, c("gender", "region", "education")]

  expect_error(energy_distance(X_cat, Y_cat), "numeric")

  # Missing variable
  expect_error(energy_distance(data$X, data$Y_good, vars = c("nonexistent")))
})

# ============================================================================
# propscore
# ============================================================================

test_that("propscore returns correct S3 class", {
  data <- create_utility_test_data()

  result <- propscore(data$X, data$Y_good, form = ~ income + age + score)

  expect_s3_class(result, "propscore")
})

test_that("propscore returns valid propensity scores", {
  data <- create_utility_test_data()

  result <- propscore(data$X, data$Y_good, form = ~ income + age + score)

  # Mean propensity scores should be between 0 and 1
  expect_true(result$mean_ps_x >= 0 && result$mean_ps_x <= 1)
  expect_true(result$mean_ps_y >= 0 && result$mean_ps_y <= 1)
})

test_that("propscore detects distinguishable distributions", {
  data <- create_utility_test_data()

  result_good <- propscore(data$X, data$Y_good, form = ~ income + age + score)
  result_poor <- propscore(data$X, data$Y_poor, form = ~ income + age + score)

  # For good synthetic, mean propensity scores should be closer to 0.5
  # For poor synthetic, they should be more separated
  diff_good <- abs(result_good$mean_ps_x - result_good$mean_ps_y)
  diff_poor <- abs(result_poor$mean_ps_x - result_poor$mean_ps_y)

  expect_true(diff_poor > diff_good)
})

test_that("propscore works with explicit formula", {
  data <- create_utility_test_data()

  # Subset to numeric only
  X_num <- data$X[, c("income", "age", "score")]
  Y_num <- data$Y_good[, c("income", "age", "score")]

  # Use explicit formula (default formula with "." requires data context)
  result <- propscore(X_num, Y_num, form = ~ income + age + score)

  expect_s3_class(result, "propscore")
})

test_that("propscore print method works", {
  data <- create_utility_test_data()
  result <- propscore(data$X, data$Y_good, form = ~ income + age)

  expect_output(print(result), "propensity")
})

test_that("propscore validates input", {
  data <- create_utility_test_data()

  # Missing formula variable
  expect_error(propscore(data$X, data$Y_good, form = ~ nonexistent))

  # Non-data.frame
  expect_error(propscore("not df", data$Y_good))
})

test_that("propscore handles different NA strategies", {
  data <- create_utility_test_data()

  # Add some NAs
  X_na <- data$X
  X_na$income[1:5] <- NA

  # Should work with remove
  result <- propscore(X_na, data$Y_good,
                      form = ~ income + age,
                      na = "remove")
  expect_s3_class(result, "propscore")

  # Should work with impute
  result2 <- propscore(X_na, data$Y_good,
                       form = ~ income + age,
                       na = "impute")
  expect_s3_class(result2, "propscore")
})

# ============================================================================
# Utility Score Consistency
# ============================================================================

test_that("utility functions agree on good vs poor synthetic", {
  data <- create_utility_test_data(n = 300, seed = 456)

  # All methods should show that Y_good is better than Y_poor

  # Hellinger: lower is better
  h_good <- hellinger(data$X, data$Y_good)$hellinger_mean
  h_poor <- hellinger(data$X, data$Y_poor)$hellinger_mean
  expect_true(h_good < h_poor)


  # Energy distance: lower is better
  e_good <- energy_distance(data$X, data$Y_good, seed = 42)$energy_distance
  e_poor <- energy_distance(data$X, data$Y_poor, seed = 42)$energy_distance
  expect_true(e_good < e_poor)
})

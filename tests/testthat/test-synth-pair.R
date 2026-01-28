# Tests for synth_pair class and S3 dispatch methods

library(testthat)

# ============================================================================
# Test Data Setup
# ============================================================================

create_synth_pair_data <- function(n = 100, seed = 42) {
  set.seed(seed)

  original <- data.frame(
    age = sample(20:70, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE),
    income = sample(c("low", "medium", "high"), n, replace = TRUE),
    score = rnorm(n, 100, 15),
    weight = runif(n, 0.5, 2),
    stringsAsFactors = FALSE
  )

  synthetic <- data.frame(
    age = sample(20:70, n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE),
    income = sample(c("low", "medium", "high"), n, replace = TRUE),
    score = rnorm(n, 100, 15),
    weight = runif(n, 0.5, 2),
    stringsAsFactors = FALSE
  )

  list(original = original, synthetic = synthetic)
}

# ============================================================================
# synth_pair Creation
# ============================================================================

test_that("synth_pair creates valid object", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    key_vars = c("age", "gender", "region"),
    target_var = "income"
  )

  expect_s3_class(pair, "synth_pair")
  expect_equal(pair$n_original, 100)
  expect_equal(pair$n_synthetic, 100)
  expect_equal(pair$key_vars, c("age", "gender", "region"))
  expect_equal(pair$target_var, "income")
})

test_that("synth_pair auto-detects variable types", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic
  )

  # Should detect categorical and numeric
  expect_true("gender" %in% pair$cat_vars || "region" %in% pair$cat_vars)
  expect_true("age" %in% pair$num_vars || "score" %in% pair$num_vars)
})

test_that("synth_pair validates inputs", {
  data <- create_synth_pair_data()

  # Non-data.frame
  expect_error(synth_pair("not df", data$synthetic))
  expect_error(synth_pair(data$original, "not df"))

  # Missing key_vars
  expect_error(
    synth_pair(data$original, data$synthetic,
               key_vars = c("nonexistent"))
  )

  # Missing target_var
  expect_error(
    synth_pair(data$original, data$synthetic,
               target_var = "nonexistent")
  )
})

test_that("synth_pair print method works", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    key_vars = c("age", "gender"),
    target_var = "income"
  )

  expect_output(print(pair), "Synthetic Data Comparison Pair")
})

test_that("synth_pair summary method works", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    key_vars = c("age", "gender"),
    target_var = "income"
  )

  summ <- summary(pair)

  expect_s3_class(summ, "summary.synth_pair")
  expect_equal(summ$n_original, 100)
  expect_output(print(summ), "Summary")
})

# ============================================================================
# S3 Dispatch: dcap.synth_pair
# ============================================================================

test_that("dcap works with synth_pair", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    key_vars = c("age", "gender", "region"),
    target_var = "income"
  )

  result <- dcap(pair)

  expect_s3_class(result, "dcap")
  expect_true(!is.na(result$dcap))
  expect_true(result$dcap >= 0 && result$dcap <= 1)
})

test_that("dcap.synth_pair validates key_vars and target_var", {
  data <- create_synth_pair_data()

  # No key_vars
  pair_no_key <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    target_var = "income"
  )
  expect_error(dcap(pair_no_key), "key_vars")

  # No target_var
  pair_no_target <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    key_vars = c("age", "gender")
  )
  expect_error(dcap(pair_no_target), "target_var")
})

# ============================================================================
# S3 Dispatch: tcap.synth_pair
# ============================================================================

test_that("tcap works with synth_pair", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    key_vars = c("age", "gender", "region"),
    target_var = "income"
  )

  result <- tcap(pair)

  expect_s3_class(result, "tcap")
  expect_true(!is.na(result$tcap_mean))
})

# ============================================================================
# S3 Dispatch: hellinger.synth_pair
# ============================================================================

test_that("hellinger works with synth_pair", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic
  )

  result <- hellinger(pair)

  expect_s3_class(result, "hellinger")
  expect_true(!is.na(result$hellinger_mean))
  expect_true(result$hellinger_mean >= 0 && result$hellinger_mean <= 1)
})

# ============================================================================
# S3 Dispatch: energy_distance.synth_pair
# ============================================================================

test_that("energy_distance works with synth_pair", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic
  )

  result <- energy_distance(pair, seed = 42)

  expect_s3_class(result, "energy_distance")
  expect_true(!is.na(result$energy_distance))
  expect_true(result$energy_distance >= 0)
})

# ============================================================================
# S3 Dispatch: propscore.synth_pair
# ============================================================================

test_that("propscore works with synth_pair", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    vars = c("age", "score")  # Only numeric for simplicity
  )

  result <- propscore(pair, form = ~ age + score)

  expect_s3_class(result, "propscore")
  expect_true(!is.na(result$ps_score))
})

# ============================================================================
# S3 Dispatch: ims.synth_pair
# ============================================================================

test_that("ims works with synth_pair", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    key_vars = c("age", "gender", "region", "income")
  )

  result <- ims(pair)

  expect_s3_class(result, "ims")
  expect_true(!is.na(result$ims))
})

# ============================================================================
# S3 Dispatch: kanonymity.synth_pair
# ============================================================================

test_that("kanonymity works with synth_pair", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    key_vars = c("age", "gender", "region")
  )

  result <- kanonymity(pair)

  expect_s3_class(result, "kanonymity")
  expect_true(!is.na(result$k_level))
  expect_true(result$k_level >= 1)
})

# ============================================================================
# S3 Dispatch: ldiversity.synth_pair
# ============================================================================

test_that("ldiversity works with synth_pair", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    key_vars = c("age", "gender", "region"),
    target_var = "income"
  )

  result <- ldiversity(pair)

  expect_s3_class(result, "ldiversity")
  expect_true(!is.na(result$distinct_l))
})

# ============================================================================
# S3 Dispatch: individual_risk.synth_pair
# ============================================================================

test_that("individual_risk works with synth_pair", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    key_vars = c("age", "gender", "region")
  )

  result <- individual_risk(pair)

  expect_s3_class(result, "individual_risk")
  expect_true(!is.na(result$mean_risk))
  expect_true(result$mean_risk >= 0 && result$mean_risk <= 1)
})

# ============================================================================
# S3 Dispatch: chisq_utility.synth_pair
# ============================================================================

test_that("chisq_utility works with synth_pair", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic
  )

  result <- chisq_utility(pair, vars = c("gender", "region"))

  expect_s3_class(result, "chisq_utility")
  expect_true(!is.na(result$chi2))
  expect_true(!is.na(result$VW))
})

# ============================================================================
# synth_pair with weights
# ============================================================================

test_that("synth_pair handles weight variables correctly", {
  data <- create_synth_pair_data()

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    weight_original = "weight",
    weight_synthetic = "weight"
  )

  expect_equal(pair$weight_original, "weight")
  expect_equal(pair$weight_synthetic, "weight")
  # Weight should be excluded from vars
  expect_false("weight" %in% pair$vars)
})

test_that("synth_pair validates weight variables", {
  data <- create_synth_pair_data()

  expect_error(
    synth_pair(
      original = data$original,
      synthetic = data$synthetic,
      weight_original = "nonexistent"
    )
  )
})

# ============================================================================
# synth_pair with holdout
# ============================================================================

test_that("synth_pair accepts holdout data", {
  data <- create_synth_pair_data()
  holdout <- data$original[1:20, ]

  pair <- synth_pair(
    original = data$original,
    synthetic = data$synthetic,
    holdout = holdout
  )

  expect_equal(nrow(pair$holdout), 20)
})

# Tests for individual_risk (Individual Re-identification Risk)

library(testthat)

# --- Setup: shared test data ---

make_test_data <- function(n = 100, seed = 42) {
  set.seed(seed)
  data.frame(
    age = sample(c("young", "middle", "old"), n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE),
    income = sample(c("low", "medium", "high"), n, replace = TRUE)
  )
}

make_unique_data <- function() {
  # Every record has a unique QI combination
  data.frame(
    id = as.character(1:10),
    val = letters[1:10]
  )
}

make_uniform_data <- function() {
  # All records share the same QI combination
  data.frame(
    key = rep("A", 20),
    val = 1:20
  )
}

# --- Class structure ---

test_that("individual_risk returns correct S3 class", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  expect_s3_class(result, "individual_risk")
})

test_that("individual_risk result has expected fields", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"),
                            threshold = 0.2)

  expect_named(result, c("risk", "risk_sample", "risk_model", "method",
                          "mean_risk", "max_risk", "n_high_risk",
                          "pct_high_risk", "threshold", "n_records",
                          "n_unique", "pct_unique", "n_ec",
                          "risk_distribution", "ec_info",
                          "key_vars", "record_indices"))
  expect_equal(result$threshold, 0.2)
  expect_equal(result$key_vars, c("age", "gender"))
  expect_equal(result$n_records, 100)
})

# --- Sample-based risk ---

test_that("sample risk equals 1/f_k for each record", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"),
                            method = "sample")
  expect_equal(result$method, "sample")

  # Compute expected risk manually
  key_sig <- apply(data[, c("age", "gender"), drop = FALSE], 1, paste, collapse = "|")
  ec_table <- table(key_sig)
  expected_risk <- 1 / as.numeric(ec_table[key_sig])
  expect_equal(result$risk, expected_risk)
})

test_that("sample risk is 1 for all-unique data", {
  data <- make_unique_data()
  result <- individual_risk(data, key_vars = "id", method = "sample")
  expect_true(all(result$risk == 1))
  expect_equal(result$max_risk, 1)
  expect_equal(result$mean_risk, 1)
})

test_that("sample risk is 1/n for uniform data", {
  data <- make_uniform_data()
  result <- individual_risk(data, key_vars = "key", method = "sample")
  expect_true(all(result$risk == 1/20))
  expect_equal(result$mean_risk, 1/20)
})

# --- Model-based risk ---

test_that("model method returns risk_model", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"),
                            method = "model")
  expect_false(is.null(result$risk_model))
})

test_that("combined method returns both sample and model risk", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"),
                            method = "both")
  expect_equal(result$method, "combined")
  expect_false(is.null(result$risk_sample))
  expect_false(is.null(result$risk_model))
})

# --- Risk bounds ---

test_that("risk values are bounded in [0, 1]", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  expect_true(all(result$risk >= 0))
  expect_true(all(result$risk <= 1))
})

test_that("risk vector has same length as n_records", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  expect_equal(length(result$risk), result$n_records)
})

# --- Threshold and high-risk counts ---

test_that("n_high_risk counts records above threshold", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"),
                            threshold = 0.1)
  expect_equal(result$n_high_risk, sum(result$risk > 0.1))
})

test_that("pct_high_risk is consistent with n_high_risk", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"),
                            threshold = 0.2)
  expected_pct <- 100 * result$n_high_risk / result$n_records
  expect_equal(result$pct_high_risk, expected_pct)
})

test_that("all records are high-risk for unique data with low threshold", {
  data <- make_unique_data()
  result <- individual_risk(data, key_vars = "id", threshold = 0.5)
  expect_equal(result$n_high_risk, 10)
})

test_that("no records are high-risk for uniform data", {
  data <- make_uniform_data()
  result <- individual_risk(data, key_vars = "key", threshold = 0.1)
  expect_equal(result$n_high_risk, 0)
})

# --- Risk distribution ---

test_that("risk_distribution contains expected quantiles", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  rd <- result$risk_distribution
  expect_true(all(c("min", "q25", "median", "mean", "q75", "q90",
                     "q95", "max") %in% names(rd)))
  expect_true(rd$min <= rd$median)
  expect_true(rd$median <= rd$max)
})

# --- Unique records ---

test_that("n_unique is correct for all-unique data", {
  data <- make_unique_data()
  result <- individual_risk(data, key_vars = "id")
  expect_equal(result$n_unique, 10)
  expect_equal(result$pct_unique, 100)
})

test_that("n_unique is zero for uniform data", {
  data <- make_uniform_data()
  result <- individual_risk(data, key_vars = "key")
  expect_equal(result$n_unique, 0)
  expect_equal(result$pct_unique, 0)
})

# --- EC info ---

test_that("ec_info is a data.frame with correct columns", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  expect_true(is.data.frame(result$ec_info))
  expect_true(all(c("key", "size", "n_records", "sample_risk") %in%
                    names(result$ec_info)))
})

test_that("ec_info sizes sum to n_records", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  expect_equal(sum(result$ec_info$size), result$n_records)
})

# --- Weight handling ---

test_that("individual_risk accepts weight as column name", {
  data <- make_test_data(n = 50)
  data$wt <- runif(50, 0.5, 2)
  result <- individual_risk(data, key_vars = c("age", "gender"),
                            weight = "wt")
  expect_s3_class(result, "individual_risk")
})

test_that("individual_risk accepts weight as numeric vector", {
  data <- make_test_data(n = 50)
  wt <- runif(50, 0.5, 2)
  result <- individual_risk(data, key_vars = c("age", "gender"),
                            weight = wt)
  expect_s3_class(result, "individual_risk")
})

test_that("individual_risk errors on bad weight column name", {
  data <- make_test_data(n = 50)
  expect_error(individual_risk(data, key_vars = c("age", "gender"),
                               weight = "nonexistent"),
               "not found")
})

test_that("individual_risk errors on wrong-length weight vector", {
  data <- make_test_data(n = 50)
  expect_error(individual_risk(data, key_vars = c("age", "gender"),
                               weight = 1:10),
               "same length")
})

# --- NA handling ---

test_that("individual_risk removes NAs by default", {
  data <- make_test_data(n = 50)
  data$age[1:5] <- NA
  result <- individual_risk(data, key_vars = c("age", "gender"), na.rm = TRUE)
  expect_equal(result$n_records, 45)
})

test_that("individual_risk errors when all cases have NAs", {
  data <- data.frame(a = c(NA, NA, NA), b = c(1, 2, 3))
  expect_error(individual_risk(data, key_vars = "a"),
               "No complete cases")
})

# --- Edge cases ---

test_that("individual_risk works with single key variable", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = "age")
  expect_s3_class(result, "individual_risk")
  expect_equal(result$n_records, 100)
})

test_that("individual_risk works with many key variables", {
  data <- make_test_data()
  result <- individual_risk(data,
                            key_vars = c("age", "gender", "region", "income"))
  expect_s3_class(result, "individual_risk")
})

test_that("more key variables leads to higher or equal mean risk", {
  data <- make_test_data(seed = 42)
  r1 <- individual_risk(data, key_vars = "age", method = "sample")
  r2 <- individual_risk(data, key_vars = c("age", "gender"), method = "sample")
  r3 <- individual_risk(data, key_vars = c("age", "gender", "region"),
                         method = "sample")
  expect_true(r1$mean_risk <= r2$mean_risk + 1e-10)
  expect_true(r2$mean_risk <= r3$mean_risk + 1e-10)
})

# --- Error handling ---

test_that("individual_risk errors on non-data.frame", {
  expect_error(individual_risk(1:10, key_vars = "a"),
               "X must be a data frame")
})

test_that("individual_risk errors on missing key variables", {
  data <- make_test_data()
  expect_error(individual_risk(data, key_vars = c("age", "nonexistent")),
               "Key variables missing")
})

# --- synth_pair method ---

test_that("individual_risk works with synth_pair (default = synthetic)", {
  set.seed(42)
  orig <- data.frame(
    age = sample(c("young", "old"), 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE),
    income = sample(c("low", "high"), 50, replace = TRUE)
  )
  syn <- orig
  syn$income <- sample(syn$income)

  pair <- synth_pair(orig, syn,
                     key_vars = c("age", "gender"),
                     target_var = "income")

  result <- individual_risk(pair)
  expect_s3_class(result, "individual_risk")
  expect_equal(result$n_records, 50)
  expect_equal(result$key_vars, c("age", "gender"))
})

test_that("individual_risk.synth_pair data='original' uses original data", {
  set.seed(42)
  orig <- data.frame(
    age = sample(c("young", "old"), 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE)
  )
  # Synthetic: all same keys -> low risk; original: varied keys -> higher risk
  syn <- data.frame(
    age = rep("young", 50),
    gender = rep("M", 50)
  )

  pair <- synth_pair(orig, syn, key_vars = c("age", "gender"))

  result_syn <- individual_risk(pair, data = "synthetic")
  result_orig <- individual_risk(pair, data = "original")

  # Synthetic is uniform -> risk = 1/50 for all
  expect_equal(result_syn$mean_risk, 1/50)
  # Original has higher risk (more diverse ECs)
  expect_true(result_orig$mean_risk > result_syn$mean_risk)
})

test_that("individual_risk.synth_pair errors without key_vars", {
  pair <- synth_pair(data.frame(a = 1:5), data.frame(a = 1:5))
  expect_error(individual_risk(pair), "key_vars")
})

# --- high_risk_records helper ---

test_that("high_risk_records returns correct records", {
  data <- make_unique_data()
  result <- individual_risk(data, key_vars = "id", threshold = 0.5)
  hr <- high_risk_records(result)
  expect_true(is.data.frame(hr))
  expect_true(all(c("index", "risk") %in% names(hr)))
  expect_equal(nrow(hr), 10)
  expect_true(all(hr$risk > 0.5))
})

test_that("high_risk_records with custom threshold", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  hr <- high_risk_records(result, threshold = 0.5)
  expect_true(is.data.frame(hr))
  expect_true(all(c("index", "risk") %in% names(hr)))
  if (nrow(hr) > 0) {
    expect_true(all(hr$risk > 0.5))
  }
})

test_that("high_risk_records errors on non-individual_risk object", {
  expect_error(high_risk_records(list()), "individual_risk")
})

# --- S3 methods: print, summary, plot ---

test_that("print.individual_risk runs without error", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  expect_output(print(result), "Individual Re-identification Risk")
  expect_output(print(result), "Key variables")
  expect_output(print(result), "Mean risk")
})

test_that("summary.individual_risk returns correct class", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  s <- summary(result)
  expect_s3_class(s, "summary.individual_risk")
  expect_true(!is.null(s$risk_bands))
  expect_true(!is.null(s$risk_distribution))
  expect_true(!is.null(s$smallest_ec))
})

test_that("print.summary.individual_risk runs without error", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  expect_output(print(summary(result)), "Summary: Individual Re-identification Risk")
  expect_output(print(summary(result)), "Risk Bands")
})

test_that("summary risk_bands sum to n_records", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  s <- summary(result)
  expect_equal(sum(s$risk_bands), result$n_records)
})

test_that("plot.individual_risk which=1 runs without error", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  expect_silent(plot(result, which = 1))
})

test_that("plot.individual_risk which=2 runs without error", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  expect_silent(plot(result, which = 2))
})

test_that("plot.individual_risk which=1:2 runs without error", {
  data <- make_test_data()
  result <- individual_risk(data, key_vars = c("age", "gender"))
  expect_silent(plot(result, which = 1:2))
})

# Tests for kanonymity (k-Anonymity Assessment)

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

test_that("kanonymity returns correct S3 class", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"))
  expect_s3_class(result, "kanonymity")
})

test_that("kanonymity result has expected fields", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"), k = 3)

  expect_named(result, c("k_level", "k_threshold", "satisfies_k",
                          "n_violating", "pct_violating", "n_records",
                          "n_ec", "ec_size_distribution",
                          "violating_records", "equivalence_classes",
                          "risk_summary", "key_vars"))
  expect_equal(result$k_threshold, 3)
  expect_equal(result$key_vars, c("age", "gender"))
  expect_equal(result$n_records, 100)
})

# --- Basic computation ---

test_that("k_level is minimum equivalence class size", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"))
  ec_sizes <- result$equivalence_classes$size
  expect_equal(result$k_level, min(ec_sizes))
})

test_that("satisfies_k is TRUE when k_level >= k", {
  data <- make_uniform_data()
  result <- kanonymity(data, key_vars = "key", k = 5)
  expect_true(result$satisfies_k)
  expect_equal(result$k_level, 20)
})

test_that("satisfies_k is FALSE when k_level < k", {
  data <- make_unique_data()
  result <- kanonymity(data, key_vars = "id", k = 5)
  expect_false(result$satisfies_k)
  expect_equal(result$k_level, 1)
})

test_that("k_level is 1 for all-unique data", {
  data <- make_unique_data()
  result <- kanonymity(data, key_vars = "id", k = 2)
  expect_equal(result$k_level, 1)
  expect_false(result$satisfies_k)
  expect_equal(result$n_violating, 10)
  expect_equal(result$pct_violating, 100)
})

test_that("k_level equals n for uniform data", {
  data <- make_uniform_data()
  result <- kanonymity(data, key_vars = "key", k = 2)
  expect_equal(result$k_level, 20)
  expect_true(result$satisfies_k)
  expect_equal(result$n_violating, 0)
  expect_equal(result$pct_violating, 0)
})

# --- Different k values ---

test_that("increasing k increases violations", {
  data <- make_test_data(n = 100, seed = 42)
  r2 <- kanonymity(data, key_vars = c("age", "gender", "region"), k = 2)
  r5 <- kanonymity(data, key_vars = c("age", "gender", "region"), k = 5)
  r10 <- kanonymity(data, key_vars = c("age", "gender", "region"), k = 10)
  expect_true(r2$n_violating <= r5$n_violating)
  expect_true(r5$n_violating <= r10$n_violating)
})

test_that("more key variables generally leads to lower k_level", {
  data <- make_test_data(seed = 42)
  r1 <- kanonymity(data, key_vars = "age")
  r2 <- kanonymity(data, key_vars = c("age", "gender"))
  r3 <- kanonymity(data, key_vars = c("age", "gender", "region"))
  expect_true(r1$k_level >= r2$k_level)
  expect_true(r2$k_level >= r3$k_level)
})

# --- Equivalence classes ---

test_that("equivalence classes data frame has correct structure", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"), k = 3)
  ec <- result$equivalence_classes
  expect_true(is.data.frame(ec))
  expect_true(all(c("key", "size", "violates_k") %in% names(ec)))
  expect_true(all(ec$violates_k == (ec$size < 3)))
})

test_that("equivalence class sizes sum to n_records", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"))
  expect_equal(sum(result$equivalence_classes$size), result$n_records)
})

test_that("n_ec matches number of equivalence classes", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"))
  expect_equal(result$n_ec, nrow(result$equivalence_classes))
})

# --- Risk summary ---

test_that("risk_summary contains expected components", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"))
  rs <- result$risk_summary
  expect_true(all(c("n_unique", "pct_unique", "n_small_ec",
                     "mean_ec_size", "median_ec_size") %in% names(rs)))
})

test_that("risk_summary n_unique counts unique records correctly", {
  data <- make_unique_data()
  result <- kanonymity(data, key_vars = "id")
  expect_equal(result$risk_summary$n_unique, 10)
})

# --- NA handling ---

test_that("kanonymity removes NAs by default", {
  data <- make_test_data(n = 50)
  data$age[1:5] <- NA
  result <- kanonymity(data, key_vars = c("age", "gender"), na.rm = TRUE)
  expect_equal(result$n_records, 45)
})

test_that("kanonymity errors when all cases have NAs", {
  data <- data.frame(a = c(NA, NA, NA), b = c(1, 2, 3))
  expect_error(kanonymity(data, key_vars = "a"),
               "No complete cases")
})

# --- Edge cases ---

test_that("kanonymity works with single key variable", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = "age")
  expect_s3_class(result, "kanonymity")
  expect_equal(result$n_records, 100)
})

test_that("kanonymity works with many key variables", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender", "region", "income"))
  expect_s3_class(result, "kanonymity")
})

test_that("kanonymity works with single-record dataset", {
  data <- data.frame(a = "x", b = "y")
  result <- kanonymity(data, key_vars = c("a", "b"), k = 2)
  expect_equal(result$n_records, 1)
  expect_equal(result$k_level, 1)
  expect_false(result$satisfies_k)
})

# --- Error handling ---

test_that("kanonymity errors on non-data.frame", {
  expect_error(kanonymity(1:10, key_vars = "a"),
               "X must be a data frame")
})

test_that("kanonymity errors on missing key variables", {
  data <- make_test_data()
  expect_error(kanonymity(data, key_vars = c("age", "nonexistent")),
               "Key variables missing")
})

# --- synth_pair method ---

test_that("kanonymity works with synth_pair objects (default = synthetic)", {
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

  result <- kanonymity(pair)
  expect_s3_class(result, "kanonymity")
  expect_equal(result$n_records, 50)
  expect_equal(result$key_vars, c("age", "gender"))
})

test_that("kanonymity.synth_pair data='original' uses original data", {
  set.seed(42)
  # Build distinct original and synthetic so results differ
  orig <- data.frame(
    age = sample(c("young", "old"), 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE)
  )
  syn <- data.frame(
    age = rep("young", 50),
    gender = rep("M", 50)
  )

  pair <- synth_pair(orig, syn, key_vars = c("age", "gender"))

  result_syn <- kanonymity(pair, data = "synthetic")
  result_orig <- kanonymity(pair, data = "original")

  # Synthetic is uniform, so k_level = 50
  expect_equal(result_syn$k_level, 50)
  # Original should have lower k_level
  expect_true(result_orig$k_level < 50)
})

test_that("kanonymity.synth_pair errors without key_vars", {
  pair <- synth_pair(data.frame(a = 1:5), data.frame(a = 1:5))
  expect_error(kanonymity(pair), "key_vars")
})

# --- S3 methods: print, summary, plot ---

test_that("print.kanonymity runs without error", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"))
  expect_output(print(result), "k-Anonymity Assessment")
  expect_output(print(result), "Key variables")
})

test_that("print.kanonymity shows violations when present", {
  data <- make_unique_data()
  result <- kanonymity(data, key_vars = "id", k = 5)
  expect_output(print(result), "Violations")
})

test_that("print.kanonymity shows no violations when none", {
  data <- make_uniform_data()
  result <- kanonymity(data, key_vars = "key", k = 2)
  expect_output(print(result), "No violations")
})

test_that("summary.kanonymity returns correct class", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"))
  s <- summary(result)
  expect_s3_class(s, "summary.kanonymity")
  expect_true(!is.null(s$ec_size_distribution))
  expect_true(!is.null(s$smallest_ec))
})

test_that("print.summary.kanonymity runs without error", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"))
  expect_output(print(summary(result)), "Summary: k-Anonymity Assessment")
})

test_that("plot.kanonymity which=1 runs without error", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"))
  expect_no_error(suppressWarnings(plot(result, which = 1)))
})

test_that("plot.kanonymity which=2 runs without error", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"))
  expect_silent(plot(result, which = 2))
})

test_that("plot.kanonymity which=1:2 runs without error", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"))
  expect_no_error(suppressWarnings(plot(result, which = 1:2)))
})

# --- Consistency checks ---

test_that("pct_violating is consistent with n_violating", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender", "region"), k = 3)
  expected_pct <- 100 * result$n_violating / result$n_records
  expect_equal(result$pct_violating, expected_pct)
})

test_that("violating_records count matches n_violating", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender", "region"), k = 5)
  expect_equal(length(result$violating_records), result$n_violating)
})

test_that("ec_size_distribution sums to n_ec", {
  data <- make_test_data()
  result <- kanonymity(data, key_vars = c("age", "gender"))
  expect_equal(sum(result$ec_size_distribution), result$n_ec)
})

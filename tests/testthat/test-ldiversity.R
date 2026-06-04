# Tests for ldiversity (l-Diversity Assessment)

library(testthat)

# --- Setup: shared test data ---

make_test_data <- function(n = 100, seed = 42) {
  set.seed(seed)
  data.frame(
    age = sample(c("young", "middle", "old"), n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S"), n, replace = TRUE),
    disease = sample(c("healthy", "cold", "flu", "covid"), n, replace = TRUE,
                     prob = c(0.4, 0.25, 0.2, 0.15))
  )
}

make_homogeneous_data <- function() {
  # One EC has a single sensitive value -> distinct_l = 1
  data.frame(
    key = c(rep("A", 10), rep("B", 10)),
    sens = c(rep("x", 10), sample(c("x", "y", "z"), 10, replace = TRUE))
  )
}

make_diverse_data <- function() {
  # Every EC has at least 3 distinct sensitive values
  set.seed(42)
  data.frame(
    key = rep(c("A", "B", "C"), each = 30),
    sens = rep(c("x", "y", "z"), times = 30)
  )
}

# --- Class structure ---

test_that("ldiversity returns correct S3 class", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_s3_class(result, "ldiversityRisk")
})

test_that("ldiversity result has expected fields", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease", l = 3)

  expect_named(result, c("distinct_l", "entropy_l", "recursive_cl",
                          "l_threshold", "c_param",
                          "satisfies_distinct_l", "satisfies_entropy_l",
                          "n_violating_distinct", "n_violating_entropy",
                          "pct_violating_distinct", "pct_violating_entropy",
                          "n_records", "n_ec", "per_ec",
                          "key_vars", "sensitive_var"))
  expect_equal(result$l_threshold, 3)
  expect_equal(result$key_vars, c("age", "gender"))
  expect_equal(result$sensitive_var, "disease")
  expect_equal(result$n_records, 100)
})

# --- Basic computation ---

test_that("distinct_l is minimum distinct values across ECs", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_equal(result$distinct_l, min(result$per_ec$n_distinct))
})

test_that("entropy_l is minimum entropy equivalent across ECs", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_equal(result$entropy_l, min(result$per_ec$entropy_l_equiv))
})

test_that("distinct l = 1 when EC has single sensitive value", {
  data <- data.frame(
    key = c(rep("A", 10), rep("B", 10)),
    sens = c(rep("x", 10), rep("y", 10))
  )
  result <- ldiversity(data, key_vars = "key", sensitive_var = "sens", l = 2)
  expect_equal(result$distinct_l, 1)
  expect_false(result$satisfies_distinct_l)
})

test_that("distinct l >= 3 when every EC has 3+ distinct values", {
  data <- make_diverse_data()
  result <- ldiversity(data, key_vars = "key", sensitive_var = "sens", l = 3)
  expect_true(result$distinct_l >= 3)
  expect_true(result$satisfies_distinct_l)
})

# --- Different l values ---

test_that("higher l threshold leads to more or equal violations", {
  data <- make_test_data(n = 200, seed = 42)
  r2 <- ldiversity(data, key_vars = "age", sensitive_var = "disease", l = 2)
  r4 <- ldiversity(data, key_vars = "age", sensitive_var = "disease", l = 4)
  # Higher l means more ECs can violate
  expect_true(r2$n_violating_distinct <= r4$n_violating_distinct)
})

test_that("entropy l-diversity is stricter than distinct l-diversity", {
  data <- make_test_data(n = 200, seed = 42)
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease", l = 2)
  # Entropy l-diversity is generally stricter for the same l
  expect_true(result$n_violating_entropy >= result$n_violating_distinct)
})

# --- Recursive (c,l)-diversity ---

test_that("recursive_cl is boolean", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease", l = 2, c = 2)
  expect_true(is.logical(result$recursive_cl))
})

test_that("recursive_cl is FALSE when one value dominates all ECs", {
  data <- data.frame(
    key = c(rep("A", 10), rep("B", 10)),
    sens = c(rep("x", 9), "y", rep("x", 9), "y")
  )
  result <- ldiversity(data, key_vars = "key", sensitive_var = "sens",
                       l = 2, c = 2)
  # r1 = 9, rest_sum = 1, so 9 < 2*1 is FALSE
  expect_false(result$recursive_cl)
})

# --- Per-EC data frame ---

test_that("per_ec data frame has correct columns", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_true(is.data.frame(result$per_ec))
  expect_true(all(c("key", "size", "n_distinct", "entropy",
                     "entropy_l_equiv", "recursive_cl",
                     "violates_distinct", "violates_entropy") %in%
                    names(result$per_ec)))
})

test_that("per_ec sizes sum to n_records", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_equal(sum(result$per_ec$size), result$n_records)
})

test_that("n_ec matches nrow of per_ec", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_equal(result$n_ec, nrow(result$per_ec))
})

test_that("entropy is non-negative for all ECs", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_true(all(result$per_ec$entropy >= 0))
})

test_that("entropy is zero when EC has only one distinct value", {
  data <- data.frame(
    key = c(rep("A", 10), rep("B", 10)),
    sens = c(rep("x", 10), rep("y", 5), rep("z", 5))
  )
  result <- ldiversity(data, key_vars = "key", sensitive_var = "sens")
  ec_a <- result$per_ec[result$per_ec$key == "A", ]
  expect_equal(ec_a$entropy, 0)
  expect_equal(ec_a$n_distinct, 1)
})

# --- NA handling ---

test_that("ldiversity removes NAs by default", {
  data <- make_test_data(n = 50)
  data$disease[1:5] <- NA
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease", na.rm = TRUE)
  expect_equal(result$n_records, 45)
})

test_that("ldiversity errors when all cases have NAs", {
  data <- data.frame(a = c(NA, NA, NA), b = c(NA, NA, NA))
  expect_error(ldiversity(data, key_vars = "a", sensitive_var = "b"),
               "No complete cases")
})

# --- Edge cases ---

test_that("ldiversity works with single key variable", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = "age", sensitive_var = "disease")
  expect_s3_class(result, "ldiversityRisk")
  expect_equal(result$n_records, 100)
})

test_that("ldiversity works with many key variables", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender", "region"),
                       sensitive_var = "disease")
  expect_s3_class(result, "ldiversityRisk")
})

test_that("ldiversity handles single EC", {
  data <- data.frame(
    key = rep("A", 20),
    sens = sample(c("x", "y", "z"), 20, replace = TRUE)
  )
  result <- ldiversity(data, key_vars = "key", sensitive_var = "sens")
  expect_equal(result$n_ec, 1)
})

# --- Error handling ---

test_that("ldiversity errors on non-data.frame", {
  expect_error(ldiversity(1:10, key_vars = "a", sensitive_var = "b"),
               "X must be a data frame")
})

test_that("ldiversity errors on missing key variables", {
  data <- make_test_data()
  expect_error(ldiversity(data, key_vars = c("age", "nonexistent"),
                          sensitive_var = "disease"),
               "Variables missing")
})

test_that("ldiversity errors on missing sensitive variable", {
  data <- make_test_data()
  expect_error(ldiversity(data, key_vars = c("age", "gender"),
                          sensitive_var = "nonexistent"),
               "Variables missing")
})

# --- synth_pair method ---

test_that("ldiversity works with synth_pair objects (default = synthetic)", {
  set.seed(42)
  orig <- data.frame(
    age = sample(c("young", "old"), 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE),
    income = sample(c("low", "medium", "high"), 50, replace = TRUE)
  )
  syn <- orig
  syn$income <- sample(syn$income)

  pair <- synth_pair(orig, syn,
                     key_vars = c("age", "gender"),
                     target_var = "income")

  result <- ldiversity(pair)
  expect_s3_class(result, "ldiversityRisk")
  expect_equal(result$sensitive_var, "income")
  expect_equal(result$n_records, 50)
})

test_that("ldiversity.synth_pair data='original' uses original data", {
  set.seed(42)
  orig <- data.frame(
    age = sample(c("young", "old"), 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE),
    income = sample(c("low", "medium", "high"), 50, replace = TRUE)
  )
  # Synthetic has only one sensitive value -> distinct_l = 1
  syn <- orig
  syn$income <- rep("low", 50)

  pair <- synth_pair(orig, syn,
                     key_vars = c("age", "gender"),
                     target_var = "income")

  result_syn <- ldiversity(pair, data = "synthetic")
  result_orig <- ldiversity(pair, data = "original")

  # Synthetic has only 1 distinct value per EC
  expect_equal(result_syn$distinct_l, 1)
  # Original should have higher diversity
  expect_true(result_orig$distinct_l > result_syn$distinct_l)
})

test_that("ldiversity.synth_pair errors without key_vars", {
  pair <- synth_pair(data.frame(a = 1:5), data.frame(a = 1:5))
  expect_error(ldiversity(pair), "key_vars")
})

test_that("ldiversity.synth_pair errors without target_var", {
  pair <- synth_pair(data.frame(a = 1:5, b = 1:5),
                     data.frame(a = 1:5, b = 1:5),
                     key_vars = "a")
  expect_error(ldiversity(pair), "target_var")
})

# --- S3 methods: print, summary, plot ---

test_that("print.ldiversity runs without error", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_output(print(result), "l-Diversity Assessment")
  expect_output(print(result), "Key variables")
  expect_output(print(result), "Sensitive variable")
})

test_that("summary.ldiversity returns correct class", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  s <- summary(result)
  expect_s3_class(s, "summary.ldiversityRisk")
  expect_true(is.data.frame(s$ec_summary))
  expect_true(is.data.frame(s$worst_ec))
})

test_that("print.summary.ldiversity runs without error", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_output(print(summary(result)), "Summary: l-Diversity Assessment")
})

test_that("plot.ldiversity which=1 runs without error", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_silent(plot(result, which = 1))
})

test_that("plot.ldiversity which=2 runs without error", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_silent(plot(result, which = 2))
})

test_that("plot.ldiversity which=1:2 runs without error", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_silent(plot(result, which = 1:2))
})

# --- Consistency checks ---

test_that("pct_violating_distinct is consistent with n_violating_distinct", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease", l = 3)
  expected_pct <- 100 * result$n_violating_distinct / result$n_records
  expect_equal(result$pct_violating_distinct, expected_pct)
})

test_that("pct_violating_entropy is consistent with n_violating_entropy", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease", l = 3)
  expected_pct <- 100 * result$n_violating_entropy / result$n_records
  expect_equal(result$pct_violating_entropy, expected_pct)
})

test_that("n_violating is zero when l = 1", {
  data <- make_test_data()
  result <- ldiversity(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease", l = 1)
  expect_equal(result$n_violating_distinct, 0)
  expect_true(result$satisfies_distinct_l)
})

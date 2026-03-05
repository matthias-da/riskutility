# Tests for delta_presence (delta-Presence Risk Assessment)

# --- Setup: shared test data ---

make_test_data <- function(n = 200, seed = 123) {
  set.seed(seed)
  X <- data.frame(
    age = sample(c("young", "middle", "old"), n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
  )
  Y <- data.frame(
    age = sample(c("young", "middle", "old"), n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
  )
  list(X = X, Y = Y)
}

# --- Class structure ---

test_that("delta_presence returns correct S3 class", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender", "region"))
  expect_s3_class(result, "delta_presence")
})

test_that("delta_presence result has expected fields", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender", "region"),
                           delta_min = 0.1, delta_max = 0.9)

  expect_named(result, c("membership_prob", "per_combination",
                          "delta_min", "delta_max",
                          "n_violations_lower", "n_violations_upper",
                          "pct_violations", "satisfies_delta",
                          "privacy_pass", "n_original", "n_synthetic",
                          "key_vars"))
  expect_equal(result$delta_min, 0.1)
  expect_equal(result$delta_max, 0.9)
  expect_equal(result$key_vars, c("age", "gender", "region"))
})

test_that("per_combination data.frame has correct columns", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"))
  expect_true(is.data.frame(result$per_combination))
  expect_named(result$per_combination,
               c("combination", "F_k", "f_k", "prob"))
})

test_that("membership_prob has same length as original data", {
  d <- make_test_data(n = 100)
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"))
  expect_length(result$membership_prob, nrow(d$X))
})

# --- Core computation ---

test_that("membership probability is f_k / F_k", {
  # Simple deterministic case
  X <- data.frame(key = c("A", "A", "A", "B", "B"))
  Y <- data.frame(key = c("A", "A", "B"))

  result <- delta_presence(X, Y, key_vars = "key")

  # A: F_k=3, f_k=2, prob = 2/3

  # B: F_k=2, f_k=1, prob = 1/2
  combo <- result$per_combination
  prob_A <- combo$prob[combo$combination == "A"]
  prob_B <- combo$prob[combo$combination == "B"]
  expect_equal(prob_A, 2/3)
  expect_equal(prob_B, 1/2)

  # Per-record probabilities
  expect_equal(result$membership_prob, c(2/3, 2/3, 2/3, 1/2, 1/2))
})

test_that("membership probability is capped at 1", {
  # f_k > F_k: synthetic overrepresents some combinations
  X <- data.frame(key = c("A", "B"))
  Y <- data.frame(key = c("A", "A", "A", "B"))

  result <- delta_presence(X, Y, key_vars = "key")

  combo <- result$per_combination
  prob_A <- combo$prob[combo$combination == "A"]
  expect_equal(prob_A, 1.0)  # 3/1 capped to 1
})

test_that("combinations in X not in Y have zero probability", {
  X <- data.frame(key = c("A", "A", "B", "C"))
  Y <- data.frame(key = c("A"))

  result <- delta_presence(X, Y, key_vars = "key")

  combo <- result$per_combination
  expect_equal(combo$prob[combo$combination == "B"], 0)
  expect_equal(combo$prob[combo$combination == "C"], 0)
  expect_equal(combo$f_k[combo$combination == "B"], 0)
})

test_that("identical data produces membership probability of 1 for all", {
  X <- data.frame(
    a = c("X", "X", "Y"),
    b = c(1, 2, 1)
  )
  Y <- X  # exact copy

  result <- delta_presence(X, Y, key_vars = c("a", "b"))
  expect_true(all(result$membership_prob == 1.0))
})

test_that("completely disjoint data produces all-zero probabilities", {
  X <- data.frame(key = c("A", "A", "B"))
  Y <- data.frame(key = c("C", "D"))

  expect_warning(
    result <- delta_presence(X, Y, key_vars = "key"),
    "QI combination.*not found in original"
  )
  expect_true(all(result$membership_prob == 0))
})

# --- Violation detection ---

test_that("delta_max violations detected correctly", {
  X <- data.frame(key = c("A", "B"))
  Y <- data.frame(key = c("A", "B"))

  result <- delta_presence(X, Y, key_vars = "key", delta_max = 0.5)

  # Both have prob = 1/1 = 1.0, which exceeds 0.5
  expect_equal(result$n_violations_upper, 2)
  expect_false(result$satisfies_delta)
  expect_false(result$privacy_pass)
})

test_that("delta_min violations detected correctly", {
  X <- data.frame(key = c("A", "A", "A", "A", "A", "B"))
  Y <- data.frame(key = c("A"))

  result <- delta_presence(X, Y, key_vars = "key", delta_min = 0.5)

  # A: prob = 1/5 = 0.2, below delta_min=0.5 -> 5 violations
  # B: prob = 0/1 = 0, below delta_min=0.5 -> 1 violation
  expect_equal(result$n_violations_lower, 6)
  expect_false(result$satisfies_delta)
})

test_that("no violations when bounds are 0 and 1", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"),
                           delta_min = 0.0, delta_max = 1.0)
  expect_true(result$satisfies_delta)
  expect_true(result$privacy_pass)
  expect_equal(result$n_violations_lower, 0)
  expect_equal(result$n_violations_upper, 0)
  expect_equal(result$pct_violations, 0)
})

test_that("satisfies_delta consistent with violation counts", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"),
                           delta_min = 0.3, delta_max = 0.7)
  expected_total <- result$n_violations_lower + result$n_violations_upper
  expect_equal(result$satisfies_delta,
               expected_total == 0)
  expect_equal(result$privacy_pass, result$satisfies_delta)
})

test_that("pct_violations is fraction of total records", {
  d <- make_test_data(n = 100)
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"),
                           delta_min = 0.3, delta_max = 0.8)
  n_total_viol <- result$n_violations_lower + result$n_violations_upper
  expect_equal(result$pct_violations, n_total_viol / result$n_original)
})

# --- Multiple key variables ---

test_that("delta_presence works with single key variable", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y, key_vars = "age")
  expect_s3_class(result, "delta_presence")
  # With 3 categories, should have at most 3 combinations
  expect_lte(nrow(result$per_combination), 3)
})

test_that("delta_presence works with multiple key variables", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y,
                           key_vars = c("age", "gender", "region"))
  expect_s3_class(result, "delta_presence")
  # 3 * 2 * 4 = 24 possible combinations
  expect_lte(nrow(result$per_combination), 24)
})

# --- Fabricated combinations warning ---

test_that("warning issued for fabricated combinations in Y", {
  X <- data.frame(key = c("A", "B"))
  Y <- data.frame(key = c("A", "C"))  # C not in X

  expect_warning(
    delta_presence(X, Y, key_vars = "key"),
    "QI combination.*not found in original"
  )
})

test_that("no warning when all Y combinations exist in X", {
  X <- data.frame(key = c("A", "B", "C"))
  Y <- data.frame(key = c("A", "B"))

  expect_no_warning(
    delta_presence(X, Y, key_vars = "key")
  )
})

# --- Edge cases ---

test_that("delta_presence works with numeric key variables", {
  X <- data.frame(a = c(1, 1, 2, 2, 3), b = c(10, 20, 10, 20, 10))
  Y <- data.frame(a = c(1, 2, 3), b = c(10, 20, 10))

  result <- delta_presence(X, Y, key_vars = c("a", "b"))
  expect_s3_class(result, "delta_presence")
  expect_true(all(result$membership_prob >= 0))
  expect_true(all(result$membership_prob <= 1))
})

test_that("delta_presence handles single record datasets", {
  X <- data.frame(key = "A")
  Y <- data.frame(key = "A")

  result <- delta_presence(X, Y, key_vars = "key")
  expect_equal(result$membership_prob, 1.0)
  expect_equal(result$n_original, 1)
  expect_equal(result$n_synthetic, 1)
})

test_that("delta_presence handles large number of unique combinations", {
  set.seed(42)
  X <- data.frame(a = seq_len(50), b = seq_len(50))
  Y <- data.frame(a = sample(50, 30), b = sample(50, 30))

  result <- delta_presence(X, Y, key_vars = c("a", "b"))
  expect_s3_class(result, "delta_presence")
  expect_length(result$membership_prob, 50)
})

test_that("n_original and n_synthetic are correct", {
  d <- make_test_data(n = 150)
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"))
  expect_equal(result$n_original, 150)
  expect_equal(result$n_synthetic, 150)
})

test_that("n_original reflects data after NA removal", {
  d <- make_test_data(n = 50)
  d$X$age[1:5] <- NA
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"),
                           na.rm = TRUE)
  expect_equal(result$n_original, 45)
})

# --- Error handling ---

test_that("delta_presence errors on non-data.frame X", {
  expect_error(delta_presence(1:10, data.frame(a = 1:5), key_vars = "a"),
               "X must be a data frame")
})

test_that("delta_presence errors on non-data.frame Y", {
  expect_error(delta_presence(data.frame(a = 1:5), 1:10, key_vars = "a"),
               "Y must be a data frame")
})

test_that("delta_presence errors on missing key variables in X", {
  X <- data.frame(a = 1:5)
  Y <- data.frame(a = 1:5, b = 1:5)
  expect_error(delta_presence(X, Y, key_vars = c("a", "b")),
               "Key variables missing in X")
})

test_that("delta_presence errors on missing key variables in Y", {
  X <- data.frame(a = 1:5, b = 1:5)
  Y <- data.frame(a = 1:5)
  expect_error(delta_presence(X, Y, key_vars = c("a", "b")),
               "Key variables missing in Y")
})

test_that("delta_presence errors on empty key_vars", {
  X <- data.frame(a = 1:5)
  Y <- data.frame(a = 1:5)
  expect_error(delta_presence(X, Y, key_vars = character(0)),
               "key_vars must be a non-empty character vector")
})

test_that("delta_presence errors on non-character key_vars", {
  X <- data.frame(a = 1:5)
  Y <- data.frame(a = 1:5)
  expect_error(delta_presence(X, Y, key_vars = 1),
               "key_vars must be a non-empty character vector")
})

test_that("delta_presence errors on invalid delta_min", {
  d <- make_test_data()
  expect_error(delta_presence(d$X, d$Y, key_vars = "age", delta_min = -0.1),
               "delta_min must be")
  expect_error(delta_presence(d$X, d$Y, key_vars = "age", delta_min = 1.5),
               "delta_min must be")
  expect_error(delta_presence(d$X, d$Y, key_vars = "age", delta_min = "abc"),
               "delta_min must be")
})

test_that("delta_presence errors on invalid delta_max", {
  d <- make_test_data()
  expect_error(delta_presence(d$X, d$Y, key_vars = "age", delta_max = -0.1),
               "delta_max must be")
  expect_error(delta_presence(d$X, d$Y, key_vars = "age", delta_max = 1.5),
               "delta_max must be")
})

test_that("delta_presence errors when delta_min > delta_max", {
  d <- make_test_data()
  expect_error(delta_presence(d$X, d$Y, key_vars = "age",
                              delta_min = 0.8, delta_max = 0.2),
               "delta_min must be less than or equal to delta_max")
})

test_that("delta_presence errors when no complete cases in X", {
  X <- data.frame(a = c(NA, NA))
  Y <- data.frame(a = c("A", "B"))
  expect_error(delta_presence(X, Y, key_vars = "a"),
               "No complete cases remaining in X")
})

test_that("delta_presence errors when no complete cases in Y", {
  X <- data.frame(a = c("A", "B"))
  Y <- data.frame(a = c(NA, NA))
  expect_error(delta_presence(X, Y, key_vars = "a"),
               "No complete cases remaining in Y")
})

# --- synth_pair method ---

test_that("delta_presence works with synth_pair objects", {
  set.seed(123)
  orig <- data.frame(
    age = sample(c("young", "old"), 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE),
    income = sample(c("low", "high"), 50, replace = TRUE)
  )
  syn <- data.frame(
    age = sample(c("young", "old"), 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE),
    income = sample(c("low", "high"), 50, replace = TRUE)
  )

  pair <- synth_pair(orig, syn,
                     key_vars = c("age", "gender"),
                     target_var = "income")

  result <- delta_presence(pair)
  expect_s3_class(result, "delta_presence")
  expect_equal(result$key_vars, c("age", "gender"))
})

test_that("delta_presence.synth_pair passes delta bounds", {
  set.seed(42)
  orig <- data.frame(
    a = sample(c("X", "Y"), 30, replace = TRUE),
    b = sample(c("M", "F"), 30, replace = TRUE)
  )
  syn <- data.frame(
    a = sample(c("X", "Y"), 30, replace = TRUE),
    b = sample(c("M", "F"), 30, replace = TRUE)
  )

  pair <- synth_pair(orig, syn, key_vars = c("a", "b"))
  result <- delta_presence(pair, delta_min = 0.2, delta_max = 0.8)
  expect_equal(result$delta_min, 0.2)
  expect_equal(result$delta_max, 0.8)
})

test_that("delta_presence.synth_pair errors without key_vars", {
  pair <- synth_pair(data.frame(a = 1:5), data.frame(a = 1:5))
  expect_error(delta_presence(pair), "key_vars")
})

# --- S3 methods: print, summary, plot ---

test_that("print.delta_presence runs without error", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"))
  expect_output(print(result), "delta-Presence Risk Assessment")
  expect_output(print(result), "Membership Probability")
})

test_that("print.delta_presence shows PASS for satisfied bounds", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"),
                           delta_min = 0.0, delta_max = 1.0)
  expect_output(print(result), "PASS")
})

test_that("print.delta_presence shows WARNING for violated bounds", {
  X <- data.frame(key = c("A", "B"))
  Y <- data.frame(key = c("A", "B"))
  result <- delta_presence(X, Y, key_vars = "key", delta_max = 0.5)
  expect_output(print(result), "WARNING")
})

test_that("summary.delta_presence returns correct class", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"))
  s <- summary(result)
  expect_s3_class(s, "summary.delta_presence")
  expect_true(is.data.frame(s$prob_summary))
  expect_true(is.data.frame(s$worst_combinations))
})

test_that("summary includes zero membership count", {
  X <- data.frame(key = c("A", "B", "C"))
  Y <- data.frame(key = c("A"))
  result <- delta_presence(X, Y, key_vars = "key")
  s <- summary(result)
  expect_equal(s$n_zero_membership, 2)  # B and C have f_k = 0
})

test_that("print.summary.delta_presence runs without error", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"))
  expect_output(print(summary(result)),
                "Summary: delta-Presence Risk Assessment")
})

test_that("plot.delta_presence which=1 runs without error", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"),
                           delta_min = 0.2, delta_max = 0.8)
  expect_silent(plot(result, which = 1))
})

test_that("plot.delta_presence which=2 runs without error", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"))
  expect_silent(plot(result, which = 2))
})

test_that("plot.delta_presence which=1:2 runs without error", {
  d <- make_test_data()
  result <- delta_presence(d$X, d$Y, key_vars = c("age", "gender"),
                           delta_min = 0.1, delta_max = 0.9)
  expect_silent(plot(result, which = 1:2))
})

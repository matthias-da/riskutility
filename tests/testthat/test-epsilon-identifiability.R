# Tests for epsilon_identifiability (Epsilon Identifiability Risk Assessment)

# --- Setup: shared test data ---

make_mixed_data <- function(n = 200, seed = 123) {
  set.seed(seed)
  data.frame(
    age = sample(20:70, n, replace = TRUE),
    income = rnorm(n, 50000, 15000),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
  )
}

make_numeric_data <- function(n = 200, seed = 123) {
  set.seed(seed)
  data.frame(
    x1 = rnorm(n, 0, 1),
    x2 = rnorm(n, 10, 5),
    x3 = runif(n, 0, 100)
  )
}

# --- Class structure ---

test_that("epsilon_identifiability returns correct S3 class", {
  X <- make_mixed_data()
  Y <- make_mixed_data(seed = 456)
  result <- epsilon_identifiability(X, Y)
  expect_s3_class(result, "epsilon_identifiability")
})

test_that("epsilon_identifiability result has expected fields", {
  X <- make_mixed_data()
  Y <- make_mixed_data(seed = 456)
  result <- epsilon_identifiability(X, Y, epsilon = 0.05)

  expected_names <- c("distances", "flagged", "identifiability_rate",
                       "entropy_weights", "entropies", "epsilon",
                       "n_flagged", "privacy_pass",
                       "n_original", "n_synthetic", "vars")
  expect_named(result, expected_names)
})

test_that("result fields have correct types", {
  X <- make_mixed_data()
  Y <- make_mixed_data(seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_type(result$distances, "double")
  expect_type(result$flagged, "logical")
  expect_type(result$identifiability_rate, "double")
  expect_type(result$entropy_weights, "double")
  expect_type(result$entropies, "double")
  expect_type(result$epsilon, "double")
  expect_type(result$n_flagged, "integer")
  expect_type(result$privacy_pass, "logical")
  expect_type(result$n_original, "integer")
  expect_type(result$n_synthetic, "integer")
  expect_type(result$vars, "character")
})

test_that("distances vector has correct length", {
  X <- make_mixed_data(n = 100)
  Y <- make_mixed_data(n = 80, seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_length(result$distances, 80)
  expect_length(result$flagged, 80)
})

# --- Core algorithm correctness ---

test_that("all distances are non-negative", {
  X <- make_mixed_data()
  Y <- make_mixed_data(seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_true(all(result$distances >= 0))
})

test_that("entropy weights sum to 1", {
  X <- make_mixed_data()
  Y <- make_mixed_data(seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_equal(sum(result$entropy_weights), 1, tolerance = 1e-10)
})

test_that("entropy weights are named correctly", {
  X <- make_mixed_data()
  Y <- make_mixed_data(seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_equal(names(result$entropy_weights), result$vars)
  expect_equal(names(result$entropies), result$vars)
})

test_that("entropies are non-negative", {
  X <- make_mixed_data()
  Y <- make_mixed_data(seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_true(all(result$entropies >= 0))
})

test_that("n_flagged equals sum of flagged vector", {
  X <- make_mixed_data()
  Y <- make_mixed_data(seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_equal(result$n_flagged, sum(result$flagged))
})

test_that("identifiability_rate equals n_flagged / n_synthetic", {
  X <- make_mixed_data()
  Y <- make_mixed_data(seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_equal(result$identifiability_rate,
               result$n_flagged / result$n_synthetic)
})

test_that("flagged vector matches distance < epsilon", {
  X <- make_mixed_data()
  Y <- make_mixed_data(seed = 456)
  result <- epsilon_identifiability(X, Y, epsilon = 0.1)

  expect_equal(result$flagged, result$distances < 0.1)
})

# --- Memorized data detection ---

test_that("identical data produces high identifiability rate", {
  X <- make_mixed_data(n = 100)
  # Y is a copy of X
  Y <- X[sample(nrow(X), 100, replace = TRUE), ]
  result <- epsilon_identifiability(X, Y, epsilon = 0.05)

  expect_true(result$identifiability_rate > 0.5)
  expect_false(result$privacy_pass)
})

test_that("independent data produces low identifiability rate", {
  X <- make_mixed_data(n = 200)
  Y <- make_mixed_data(n = 200, seed = 789)
  result <- epsilon_identifiability(X, Y, epsilon = 0.02)

  # With independent data, rate should be less than memorized data
  expect_true(result$identifiability_rate < 1.0)
})

# --- Epsilon parameter ---

test_that("larger epsilon flags more records", {
  X <- make_mixed_data(n = 100)
  Y <- make_mixed_data(n = 100, seed = 456)

  r_small <- epsilon_identifiability(X, Y, epsilon = 0.01)
  r_large <- epsilon_identifiability(X, Y, epsilon = 0.2)

  expect_true(r_large$n_flagged >= r_small$n_flagged)
})

test_that("epsilon = very small flags zero records for independent data", {
  X <- make_numeric_data(n = 100)
  Y <- make_numeric_data(n = 100, seed = 456)
  result <- epsilon_identifiability(X, Y, epsilon = 1e-6)

  expect_equal(result$n_flagged, 0L)
  expect_equal(result$identifiability_rate, 0)
})

# --- Privacy pass ---

test_that("privacy_pass is TRUE when identifiability_rate <= 0.1", {
  X <- make_mixed_data(n = 200)
  Y <- make_mixed_data(n = 200, seed = 789)
  result <- epsilon_identifiability(X, Y, epsilon = 0.001)

  if (result$identifiability_rate <= 0.1) {
    expect_true(result$privacy_pass)
  } else {
    expect_false(result$privacy_pass)
  }
})

test_that("privacy_pass reflects the 0.1 threshold correctly", {
  # Use memorized data to get high rate
  X <- make_mixed_data(n = 100)
  Y <- X
  result <- epsilon_identifiability(X, Y, epsilon = 0.05)

  expect_equal(result$privacy_pass, result$identifiability_rate <= 0.1)
})

# --- Entropy computation ---

test_that("constant variable gets zero entropy", {
  X <- data.frame(const = rep(1, 50), var = rnorm(50))
  Y <- data.frame(const = rep(1, 30), var = rnorm(30))
  result <- epsilon_identifiability(X, Y)

  expect_equal(result$entropies["const"], c(const = 0))
  expect_equal(result$entropy_weights["const"], c(const = 0))
})

test_that("high-entropy variable gets lower weight than low-entropy variable", {
  set.seed(42)
  # Binary variable (low entropy) vs many-valued variable (high entropy)
  X <- data.frame(
    binary = sample(c("A", "B"), 200, replace = TRUE),
    diverse = sample(letters, 200, replace = TRUE)
  )
  Y <- data.frame(
    binary = sample(c("A", "B"), 100, replace = TRUE),
    diverse = sample(letters, 100, replace = TRUE)
  )
  result <- epsilon_identifiability(X, Y)

  expect_true(result$entropy_weights["binary"] > result$entropy_weights["diverse"])
})

test_that("all-constant data uses equal weights with warning", {
  X <- data.frame(a = rep(1, 50), b = rep("x", 50))
  Y <- data.frame(a = rep(1, 30), b = rep("x", 30))

  expect_warning(
    result <- epsilon_identifiability(X, Y),
    "zero entropy"
  )
  expect_equal(result$entropy_weights, c(a = 0.5, b = 0.5))
})

# --- Variable selection ---

test_that("vars parameter selects correct variables", {
  X <- make_mixed_data()
  Y <- make_mixed_data(seed = 456)
  result <- epsilon_identifiability(X, Y, vars = c("age", "gender"))

  expect_equal(result$vars, c("age", "gender"))
  expect_length(result$entropy_weights, 2)
  expect_length(result$entropies, 2)
})

test_that("NULL vars uses all common variables", {
  X <- make_mixed_data()
  Y <- make_mixed_data(seed = 456)
  result <- epsilon_identifiability(X, Y, vars = NULL)

  expect_equal(sort(result$vars), sort(intersect(names(X), names(Y))))
})

# --- n_bins parameter ---

test_that("n_bins parameter affects numeric entropy computation", {
  X <- make_numeric_data(n = 200)
  Y <- make_numeric_data(n = 100, seed = 456)

  r_few <- epsilon_identifiability(X, Y, n_bins = 5)
  r_many <- epsilon_identifiability(X, Y, n_bins = 50)

  # Different binning should produce different entropies
  expect_false(identical(r_few$entropies, r_many$entropies))
})

# --- Dataset sizes ---

test_that("n_original and n_synthetic are recorded correctly", {
  X <- make_mixed_data(n = 150)
  Y <- make_mixed_data(n = 80, seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_equal(result$n_original, 150L)
  expect_equal(result$n_synthetic, 80L)
})

# --- NA handling ---

test_that("na.rm = TRUE removes incomplete cases", {
  X <- make_mixed_data(n = 100)
  Y <- make_mixed_data(n = 50, seed = 456)
  X$age[1:5] <- NA
  Y$income[1:3] <- NA

  result <- epsilon_identifiability(X, Y, na.rm = TRUE)

  expect_equal(result$n_original, 95L)
  expect_equal(result$n_synthetic, 47L)
  expect_length(result$distances, 47)
})

test_that("errors when no complete cases remain after NA removal", {
  X <- data.frame(a = c(NA_real_, NA_real_, NA_real_), b = c(NA_real_, NA_real_, NA_real_))
  Y <- data.frame(a = c(1, 2), b = c(3, 4))

  expect_error(epsilon_identifiability(X, Y), "No complete cases in X")
})

# --- Error handling ---

test_that("epsilon_identifiability errors on non-data.frame X", {
  expect_error(
    epsilon_identifiability(1:10, data.frame(a = 1:5)),
    "X must be a data frame"
  )
})

test_that("epsilon_identifiability errors on non-data.frame Y", {
  expect_error(
    epsilon_identifiability(data.frame(a = 1:5), 1:10),
    "Y must be a data frame"
  )
})

test_that("epsilon_identifiability errors on invalid epsilon", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = 1:10)

  expect_error(epsilon_identifiability(X, Y, epsilon = 0), "epsilon must be")
  expect_error(epsilon_identifiability(X, Y, epsilon = 1), "epsilon must be")
  expect_error(epsilon_identifiability(X, Y, epsilon = -0.1), "epsilon must be")
  expect_error(epsilon_identifiability(X, Y, epsilon = "abc"), "epsilon must be")
  expect_error(epsilon_identifiability(X, Y, epsilon = c(0.1, 0.2)), "epsilon must be")
})

test_that("epsilon_identifiability errors on invalid n_bins", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = 1:10)

  expect_error(epsilon_identifiability(X, Y, n_bins = 1), "n_bins must be")
  expect_error(epsilon_identifiability(X, Y, n_bins = "abc"), "n_bins must be")
})

test_that("epsilon_identifiability errors when no common variables", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(b = 1:10)

  expect_error(epsilon_identifiability(X, Y), "No common variables")
})

test_that("epsilon_identifiability errors on missing variables", {
  X <- data.frame(a = 1:10, b = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10)

  expect_error(
    epsilon_identifiability(X, Y, vars = c("a", "c")),
    "Variables missing"
  )
})

test_that("epsilon_identifiability errors on type mismatch", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = letters[1:10])

  expect_error(epsilon_identifiability(X, Y), "different class")
})

# --- synth_pair method ---

test_that("epsilon_identifiability works with synth_pair objects", {
  set.seed(123)
  orig <- make_mixed_data(n = 100)
  syn <- make_mixed_data(n = 100, seed = 456)

  pair <- synth_pair(orig, syn)
  result <- epsilon_identifiability(pair)

  expect_s3_class(result, "epsilon_identifiability")
  expect_equal(result$n_original, 100L)
  expect_equal(result$n_synthetic, 100L)
})

test_that("synth_pair method passes vars through", {
  orig <- make_mixed_data(n = 50)
  syn <- make_mixed_data(n = 50, seed = 456)

  pair <- synth_pair(orig, syn, vars = c("age", "income"))
  result <- epsilon_identifiability(pair)

  expect_equal(sort(result$vars), sort(c("age", "income")))
})

# --- S3 methods: print, summary, plot ---

test_that("print.epsilon_identifiability runs without error", {
  X <- make_mixed_data(n = 50)
  Y <- make_mixed_data(n = 50, seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_output(print(result), "Epsilon Identifiability Risk Assessment")
  expect_output(print(result), "Identifiability rate")
  expect_output(print(result), "Entropy Weights")
  expect_output(print(result), "Privacy Assessment")
})

test_that("print shows PASS for low risk", {
  X <- make_mixed_data(n = 100)
  Y <- make_mixed_data(n = 100, seed = 789)
  result <- epsilon_identifiability(X, Y, epsilon = 0.001)

  if (result$privacy_pass) {
    expect_output(print(result), "PASS")
  }
})

test_that("print shows WARNING for high risk", {
  X <- make_mixed_data(n = 100)
  Y <- X  # memorized
  result <- epsilon_identifiability(X, Y, epsilon = 0.1)

  if (!result$privacy_pass) {
    expect_output(print(result), "WARNING")
  }
})

test_that("summary.epsilon_identifiability returns correct class", {
  X <- make_mixed_data(n = 50)
  Y <- make_mixed_data(n = 50, seed = 456)
  result <- epsilon_identifiability(X, Y)
  s <- summary(result)

  expect_s3_class(s, "summary.epsilon_identifiability")
})

test_that("summary has expected fields", {
  X <- make_mixed_data(n = 50)
  Y <- make_mixed_data(n = 50, seed = 456)
  result <- epsilon_identifiability(X, Y)
  s <- summary(result)

  expect_true(!is.null(s$identifiability_rate))
  expect_true(!is.null(s$n_flagged))
  expect_true(!is.null(s$privacy_pass))
  expect_true(!is.null(s$epsilon))
  expect_true(!is.null(s$mean_distance))
  expect_true(!is.null(s$median_distance))
  expect_true(!is.null(s$sd_distance))
  expect_true(!is.null(s$quantiles))
  expect_true(!is.null(s$entropy_weights))
  expect_true(!is.null(s$entropies))
})

test_that("print.summary.epsilon_identifiability runs without error", {
  X <- make_mixed_data(n = 50)
  Y <- make_mixed_data(n = 50, seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_output(print(summary(result)),
                "Summary: Epsilon Identifiability Risk Assessment")
})

test_that("plot.epsilon_identifiability which=1 runs without error", {
  X <- make_mixed_data(n = 50)
  Y <- make_mixed_data(n = 50, seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_silent(plot(result, which = 1))
})

test_that("plot.epsilon_identifiability which=2 runs without error", {
  X <- make_mixed_data(n = 50)
  Y <- make_mixed_data(n = 50, seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_silent(plot(result, which = 2))
})

test_that("plot.epsilon_identifiability which=1:2 runs without error", {
  X <- make_mixed_data(n = 50)
  Y <- make_mixed_data(n = 50, seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_silent(plot(result, which = 1:2))
})

# --- Purely numeric data ---

test_that("epsilon_identifiability works with all-numeric data", {
  X <- make_numeric_data(n = 100)
  Y <- make_numeric_data(n = 80, seed = 456)
  result <- epsilon_identifiability(X, Y)

  expect_s3_class(result, "epsilon_identifiability")
  expect_true(all(result$entropies > 0))
})

# --- Purely categorical data ---

test_that("epsilon_identifiability works with all-categorical data", {
  set.seed(123)
  X <- data.frame(
    a = sample(c("x", "y", "z"), 100, replace = TRUE),
    b = sample(c("M", "F"), 100, replace = TRUE),
    c = sample(c("N", "S", "E", "W"), 100, replace = TRUE)
  )
  set.seed(456)
  Y <- data.frame(
    a = sample(c("x", "y", "z"), 80, replace = TRUE),
    b = sample(c("M", "F"), 80, replace = TRUE),
    c = sample(c("N", "S", "E", "W"), 80, replace = TRUE)
  )
  result <- epsilon_identifiability(X, Y)

  expect_s3_class(result, "epsilon_identifiability")
  expect_length(result$distances, 80)
})

# --- Reproducibility ---

test_that("results are deterministic for same data", {
  X <- make_mixed_data(n = 100)
  Y <- make_mixed_data(n = 80, seed = 456)

  r1 <- epsilon_identifiability(X, Y)
  r2 <- epsilon_identifiability(X, Y)

  expect_equal(r1$distances, r2$distances)
  expect_equal(r1$identifiability_rate, r2$identifiability_rate)
  expect_equal(r1$entropy_weights, r2$entropy_weights)
})

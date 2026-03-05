# Tests for hitting_rate (Hitting Rate)

# --- Class structure and fields ---

test_that("hitting_rate returns correct S3 class structure", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    age = rnorm(n, 40, 10),
    income = rnorm(n, 50000, 15000)
  )
  Y <- data.frame(
    age = rnorm(80, 40, 10),
    income = rnorm(80, 50000, 15000)
  )

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_s3_class(result, "hitting_rate")
  expect_true("rate" %in% names(result))
  expect_true("min_distances" %in% names(result))
  expect_true("hits" %in% names(result))
  expect_true("n_hits" %in% names(result))
  expect_true("threshold" %in% names(result))
  expect_true("privacy_pass" %in% names(result))
  expect_true("n_original" %in% names(result))
  expect_true("n_synthetic" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("vars" %in% names(result))
  expect_true("rate_at_zero" %in% names(result))
})

test_that("hitting_rate privacy_pass is logical", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_type(result$privacy_pass, "logical")
})

test_that("hitting_rate hits is logical vector of correct length", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_type(result$hits, "logical")
  expect_equal(length(result$hits), 40)
  expect_equal(length(result$min_distances), 40)
})

# --- Rate values and consistency ---

test_that("rate equals n_hits / n_synthetic", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_equal(result$rate, result$n_hits / result$n_synthetic)
})

test_that("n_hits equals sum of hits", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_equal(result$n_hits, sum(result$hits))
})

test_that("hits correspond to min_distances <= threshold", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  result <- hitting_rate(X, Y, threshold = 0.1, method = "euclidean")

  expect_equal(result$hits, result$min_distances <= result$threshold)
})

test_that("rate is between 0 and 1", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_true(result$rate >= 0 && result$rate <= 1)
})

test_that("rate_at_zero is between 0 and rate", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_true(result$rate_at_zero >= 0)
  expect_true(result$rate_at_zero <= result$rate)
})

# --- Min distances ---

test_that("min_distances are non-negative", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_true(all(result$min_distances >= 0))
})

test_that("min_distances has correct length", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_equal(length(result$min_distances), 40)
})

# --- Threshold behavior ---

test_that("higher threshold yields higher or equal rate", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  r1 <- hitting_rate(X, Y, threshold = 0.01, method = "euclidean")
  r2 <- hitting_rate(X, Y, threshold = 0.05, method = "euclidean")
  r3 <- hitting_rate(X, Y, threshold = 0.2, method = "euclidean")

  expect_true(r1$rate <= r2$rate)
  expect_true(r2$rate <= r3$rate)
})

test_that("threshold = 0 gives only exact matches", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, threshold = 0, method = "euclidean")

  # For random continuous data, no exact matches expected
  expect_equal(result$rate, 0)
  expect_equal(result$n_hits, 0)
  expect_equal(result$rate, result$rate_at_zero)
})

test_that("very large threshold hits all records", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(30), x2 = rnorm(30))
  Y <- data.frame(x1 = rnorm(20), x2 = rnorm(20))

  result <- hitting_rate(X, Y, threshold = 1e6, method = "euclidean")

  expect_equal(result$rate, 1)
  expect_equal(result$n_hits, 20)
})

test_that("default threshold is 0.05", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_equal(result$threshold, 0.05)
})

# --- Memorization detection ---

test_that("copied data has high hitting rate", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  # Y is exact copy of X
  Y <- X[sample(nrow(X), n, replace = TRUE), ]

  result <- hitting_rate(X, Y, threshold = 0.05, method = "euclidean")

  # All records should be hits (exact copies have distance 0)
  expect_equal(result$rate, 1)
  expect_equal(result$n_hits, n)
  expect_equal(result$rate_at_zero, 1)
  expect_false(result$privacy_pass)
})

test_that("independent random data has low hitting rate with small threshold", {
  set.seed(42)
  n <- 200
  X <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )
  Y <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )

  result <- hitting_rate(X, Y, threshold = 0.01, method = "euclidean")

  # Random 3D data should have very few near hits at tau=0.01
  expect_true(result$rate < 0.5)
})

# --- Distance methods ---

test_that("hitting_rate works with gower method", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    age = rnorm(n, 40, 10),
    gender = sample(c("M", "F"), n, replace = TRUE),
    income = rnorm(n, 50000, 15000)
  )
  Y <- data.frame(
    age = rnorm(40, 40, 10),
    gender = sample(c("M", "F"), 40, replace = TRUE),
    income = rnorm(40, 50000, 15000)
  )

  result <- hitting_rate(X, Y, method = "gower")

  expect_s3_class(result, "hitting_rate")
  expect_equal(result$method, "gower")
  expect_true(result$rate >= 0 && result$rate <= 1)
})

test_that("hitting_rate works with euclidean method", {
  set.seed(123)
  n <- 60
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_s3_class(result, "hitting_rate")
  expect_equal(result$method, "euclidean")
})

test_that("euclidean errors with non-numeric data", {
  X <- data.frame(x = c("a", "b", "c", "d", "e"))
  Y <- data.frame(x = c("a", "b", "c"))

  expect_error(hitting_rate(X, Y, method = "euclidean"), "numeric")
})

# --- Input validation ---

test_that("hitting_rate errors when X is not a data frame", {
  Y <- data.frame(x1 = rnorm(10))
  expect_error(hitting_rate(1:10, Y), "X must be a data frame")
})

test_that("hitting_rate errors when Y is not a data frame", {
  X <- data.frame(x1 = rnorm(10))
  expect_error(hitting_rate(X, 1:10), "Y must be a data frame")
})

test_that("hitting_rate errors with negative threshold", {
  X <- data.frame(x1 = rnorm(10))
  Y <- data.frame(x1 = rnorm(10))
  expect_error(hitting_rate(X, Y, threshold = -0.1), "non-negative")
})

test_that("hitting_rate errors with non-numeric threshold", {
  X <- data.frame(x1 = rnorm(10))
  Y <- data.frame(x1 = rnorm(10))
  expect_error(hitting_rate(X, Y, threshold = "abc"), "non-negative")
})

test_that("hitting_rate errors when no common variables", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(b = 1:10)

  expect_error(hitting_rate(X, Y), "No common variables")
})

test_that("hitting_rate errors with mismatched variable types", {
  X <- data.frame(x = 1:10)
  Y <- data.frame(x = letters[1:10])

  expect_error(hitting_rate(X, Y), "different class")
})

test_that("hitting_rate errors with missing variables", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10)

  expect_error(hitting_rate(X, Y, vars = c("a", "b")), "Variables missing in X")
})

# --- Variable selection ---

test_that("hitting_rate uses specified vars", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))

  result <- hitting_rate(X, Y, vars = c("x1", "x2"), method = "euclidean")

  expect_equal(result$vars, c("x1", "x2"))
})

test_that("hitting_rate auto-selects common vars", {
  set.seed(123)
  X <- data.frame(a = rnorm(60), b = rnorm(60), c_only = rnorm(60))
  Y <- data.frame(a = rnorm(40), b = rnorm(40), d_only = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_equal(sort(result$vars), c("a", "b"))
})

# --- NA handling ---

test_that("hitting_rate handles NA values with na.rm = TRUE", {
  set.seed(123)
  X <- data.frame(x1 = c(rnorm(18), NA, NA), x2 = rnorm(20))
  Y <- data.frame(x1 = c(rnorm(9), NA), x2 = rnorm(10))

  result <- hitting_rate(X, Y, method = "euclidean", na.rm = TRUE)

  expect_s3_class(result, "hitting_rate")
  # NA records removed: X has 18 complete, Y has 9 complete
  expect_equal(result$n_original, 18)
  expect_equal(result$n_synthetic, 9)
})

# --- Dataset sizes stored correctly ---

test_that("hitting_rate stores correct dataset sizes", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
  Y <- data.frame(x1 = rnorm(80), x2 = rnorm(80))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_equal(result$n_original, 100)
  expect_equal(result$n_synthetic, 80)
})

# --- synth_pair method ---

test_that("hitting_rate.synth_pair works", {
  set.seed(123)
  n <- 60
  original <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  synthetic <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  pair <- synth_pair(original = original, synthetic = synthetic)
  result <- hitting_rate(pair, method = "euclidean")

  expect_s3_class(result, "hitting_rate")
  expect_true(!is.null(result$rate))
})

test_that("hitting_rate.synth_pair uses vars from pair", {
  set.seed(123)
  original <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  synthetic <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))

  pair <- synth_pair(original = original, synthetic = synthetic,
                     vars = c("x1", "x2"))
  result <- hitting_rate(pair, method = "euclidean")

  expect_equal(result$vars, c("x1", "x2"))
})

# --- print, summary, plot ---

test_that("hitting_rate print method works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_output(print(result), "Hitting Rate Privacy Metric")
  expect_output(print(result), "Hitting Rate Results")
  expect_output(print(result), "Privacy Assessment")
  expect_output(print(result), "Min Distance Distribution")
})

test_that("hitting_rate summary returns correct class", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")
  s <- summary(result)

  expect_s3_class(s, "summary.hitting_rate")
  expect_true(!is.null(s$rate))
  expect_true(!is.null(s$threshold))
  expect_true(!is.null(s$privacy_pass))
  expect_true(!is.null(s$quantiles))
  expect_true(!is.null(s$rates_at_thresholds))
  expect_true(!is.null(s$mean_min_distance))
  expect_true(!is.null(s$sd_min_distance))
})

test_that("hitting_rate print.summary works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")
  s <- summary(result)

  expect_output(print(s), "Summary: Hitting Rate Privacy Metric")
  expect_output(print(s), "Key Metrics")
  expect_output(print(s), "Min Distance Distribution")
  expect_output(print(s), "Hitting Rate at Standard Thresholds")
})

test_that("hitting_rate plot method works without error for which = 1", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_no_error(plot(result, which = 1))
})

test_that("hitting_rate plot method works without error for which = 2", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_no_error(plot(result, which = 2))
})

test_that("hitting_rate plot method works without error for which = 1:2", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- hitting_rate(X, Y, method = "euclidean")

  expect_no_error(plot(result, which = 1:2))
})

# --- privacy_pass logic ---

test_that("privacy_pass is TRUE when rate <= 0.1", {
  set.seed(42)
  n <- 200
  # High-dimensional data: very few near hits expected
  X <- data.frame(
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
    x4 = rnorm(n), x5 = rnorm(n)
  )
  Y <- data.frame(
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
    x4 = rnorm(n), x5 = rnorm(n)
  )

  result <- hitting_rate(X, Y, threshold = 0.01, method = "euclidean")

  # With small threshold and high dimensions, rate should be low
  if (result$rate <= 0.1) {
    expect_true(result$privacy_pass)
  } else {
    expect_false(result$privacy_pass)
  }
})

test_that("privacy_pass is FALSE for memorized data", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  Y <- X  # exact copy

  result <- hitting_rate(X, Y, threshold = 0.05, method = "euclidean")

  expect_false(result$privacy_pass)
})

# --- Gower distances bounded in [0,1] ---

test_that("gower min_distances are bounded in [0,1]", {
  set.seed(123)
  X <- data.frame(
    age = rnorm(60, 40, 10),
    gender = sample(c("M", "F"), 60, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), 60, replace = TRUE)
  )
  Y <- data.frame(
    age = rnorm(40, 40, 10),
    gender = sample(c("M", "F"), 40, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), 40, replace = TRUE)
  )

  result <- hitting_rate(X, Y, method = "gower")

  expect_true(all(result$min_distances >= 0))
  expect_true(all(result$min_distances <= 1))
})

# --- Internal helper ---

test_that(".compute_rates_at_thresholds returns correct format", {
  distances <- c(0, 0.005, 0.01, 0.03, 0.07, 0.12, 0.25, 0.5)
  rates <- riskutility:::.compute_rates_at_thresholds(distances)

  expect_type(rates, "double")
  expect_equal(length(rates), 8)
  expect_true(all(rates >= 0 & rates <= 1))
  # Monotonically non-decreasing
  expect_true(all(diff(rates) >= 0))
  # At tau=0, only exact 0 counts
  expect_equal(unname(rates[1]), 1 / 8)  # one distance == 0
})

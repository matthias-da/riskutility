# Tests for linkability (Linkability Risk)

# --- Class structure and fields ---

test_that("linkability returns correct S3 class structure", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    age = sample(20:70, n, replace = TRUE),
    income = rnorm(n, 50000, 15000),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
  )
  Y <- data.frame(
    age = sample(20:70, n, replace = TRUE),
    income = rnorm(n, 50000, 15000),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
  )

  result <- linkability(X, Y, n_attacks = 100, seed = 42)

  expect_s3_class(result, "linkability")
  expect_true("risk" %in% names(result))
  expect_true("risk_ci" %in% names(result))
  expect_true("risk_attack" %in% names(result))
  expect_true("risk_attack_ci" %in% names(result))
  expect_true("risk_control" %in% names(result))
  expect_true("risk_control_ci" %in% names(result))
  expect_true("n_attacks" %in% names(result))
  expect_true("n_success" %in% names(result))
  expect_true("n_control_success" %in% names(result))
  expect_true("links" %in% names(result))
  expect_true("links_control" %in% names(result))
  expect_true("privacy_pass" %in% names(result))
  expect_true("n_original" %in% names(result))
  expect_true("n_synthetic" %in% names(result))
  expect_true("n_holdout" %in% names(result))
  expect_true("aux_cols" %in% names(result))
  expect_true("secret_cols" %in% names(result))
  expect_true("n_neighbors" %in% names(result))
  expect_true("vars" %in% names(result))
  expect_true("confidence_level" %in% names(result))
})

test_that("linkability privacy_pass is logical", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- linkability(X, Y, n_attacks = 50, seed = 42)

  expect_type(result$privacy_pass, "logical")
})

test_that("linkability risk is between 0 and 1", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))

  result <- linkability(X, Y, n_attacks = 100, seed = 42)

  expect_true(result$risk >= 0 && result$risk <= 1)
  expect_true(result$risk_attack >= 0 && result$risk_attack <= 1)
  expect_true(result$risk_control >= 0 && result$risk_control <= 1)
})

test_that("links vector has correct length", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))

  result <- linkability(X, Y, n_attacks = 30, seed = 42)

  expect_equal(length(result$links), result$n_attacks)
  expect_equal(length(result$links_control), result$n_attacks)
  expect_type(result$links, "logical")
})

# --- Risk values for random vs memorized data ---

test_that("attack rate is higher for copied/memorized data", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    x1 = seq_len(n),
    x2 = seq_len(n) * 2,
    x3 = rnorm(n),
    x4 = rnorm(n)
  )
  # Y is exact copy - maximum memorization
  Y_copy <- X

  result_copy <- linkability(X, Y_copy, n_attacks = 100, seed = 42)

  # Independent data
  Y_rand <- data.frame(
    x1 = sample(seq_len(n)),
    x2 = sample(seq_len(n) * 2),
    x3 = rnorm(n),
    x4 = rnorm(n)
  )
  result_rand <- linkability(X, Y_rand, n_attacks = 100, seed = 42)

  # Copied data should have higher attack rate

  expect_true(result_copy$risk_attack >= result_rand$risk_attack)
})

test_that("attack rate exceeds control rate for memorized data", {
  set.seed(123)
  n <- 80
  train <- data.frame(
    x1 = seq_len(n),
    x2 = seq_len(n) * 3,
    x3 = rnorm(n),
    x4 = rnorm(n)
  )
  holdout <- data.frame(
    x1 = (n + 1):(2 * n),
    x2 = ((n + 1):(2 * n)) * 3,
    x3 = rnorm(n),
    x4 = rnorm(n)
  )
  # Synthetic is copy of training
  Y <- train

  result <- linkability(train, Y, holdout = holdout,
                         n_attacks = 80, seed = 42)

  expect_true(result$n_success >= result$n_control_success)
})

# --- Column splitting ---

test_that("random column split produces non-overlapping groups", {
  set.seed(123)
  X <- data.frame(a = 1:20, b = 1:20, c = 1:20, d = 1:20)
  Y <- data.frame(a = 1:10, b = 1:10, c = 1:10, d = 1:10)

  result <- linkability(X, Y, n_attacks = 10, seed = 42)

  expect_true(length(result$aux_cols) > 0)
  expect_true(length(result$secret_cols) > 0)
  expect_equal(length(intersect(result$aux_cols, result$secret_cols)), 0)
  expect_equal(sort(c(result$aux_cols, result$secret_cols)),
               sort(result$vars))
})

test_that("user-specified aux_cols are respected", {
  set.seed(123)
  X <- data.frame(a = rnorm(40), b = rnorm(40), c = rnorm(40), d = rnorm(40))
  Y <- data.frame(a = rnorm(20), b = rnorm(20), c = rnorm(20), d = rnorm(20))

  result <- linkability(X, Y, aux_cols = c("a", "b"),
                         n_attacks = 20, seed = 42)

  expect_equal(sort(result$aux_cols), c("a", "b"))
  expect_equal(sort(result$secret_cols), c("c", "d"))
})

test_that("error when aux_cols covers all variables", {
  X <- data.frame(a = 1:10, b = 1:10)
  Y <- data.frame(a = 1:5, b = 1:5)

  expect_error(linkability(X, Y, aux_cols = c("a", "b")),
               "secret_cols is empty")
})

test_that("error when aux_cols not in vars", {
  X <- data.frame(a = 1:10, b = 1:10)
  Y <- data.frame(a = 1:5, b = 1:5)

  expect_error(linkability(X, Y, aux_cols = "z"),
               "aux_cols not in vars")
})

# --- n_neighbors ---

test_that("n_neighbors = 1 works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  Y <- data.frame(x1 = rnorm(30), x2 = rnorm(30), x3 = rnorm(30))

  result <- linkability(X, Y, n_neighbors = 1, n_attacks = 30, seed = 42)

  expect_s3_class(result, "linkability")
  expect_equal(result$n_neighbors, 1)
})

test_that("n_neighbors > 1 increases success rate", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    x1 = seq_len(n),
    x2 = rnorm(n),
    x3 = rnorm(n),
    x4 = rnorm(n)
  )
  Y <- X  # exact copy

  r1 <- linkability(X, Y, n_neighbors = 1, n_attacks = 50, seed = 42)
  r5 <- linkability(X, Y, n_neighbors = 5, n_attacks = 50, seed = 42)

  # More neighbors = more lenient matching = at least as many successes
  expect_true(r5$n_success >= r1$n_success)
})

# --- Holdout handling ---

test_that("linkability works with explicit holdout", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))
  holdout <- data.frame(x1 = rnorm(30), x2 = rnorm(30), x3 = rnorm(30))

  result <- linkability(X, Y, holdout = holdout, n_attacks = 30, seed = 42)

  expect_s3_class(result, "linkability")
  expect_equal(result$n_original, nrow(X))
  expect_equal(result$n_holdout, nrow(holdout))
})

test_that("linkability works with holdout_fraction", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))

  result <- linkability(X, Y, holdout_fraction = 0.3,
                         n_attacks = 30, seed = 42)

  expect_s3_class(result, "linkability")
  expect_equal(result$n_holdout, floor(n * 0.3))
  expect_equal(result$n_original, n - floor(n * 0.3))
})

test_that("linkability seed produces reproducible results", {
  set.seed(99)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))

  r1 <- linkability(X, Y, n_attacks = 50, seed = 42)
  r2 <- linkability(X, Y, n_attacks = 50, seed = 42)

  expect_equal(r1$risk, r2$risk)
  expect_equal(r1$n_success, r2$n_success)
  expect_equal(r1$links, r2$links)
  expect_equal(r1$aux_cols, r2$aux_cols)
})

# --- synth_pair method ---

test_that("linkability.synth_pair works", {
  set.seed(123)
  original <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  synthetic <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))

  pair <- synth_pair(original = original, synthetic = synthetic)
  result <- linkability(pair, n_attacks = 30, seed = 42)

  expect_s3_class(result, "linkability")
  expect_true(!is.null(result$risk))
})

test_that("linkability.synth_pair uses holdout when available", {
  set.seed(123)
  original <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  synthetic <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))
  holdout <- data.frame(x1 = rnorm(30), x2 = rnorm(30), x3 = rnorm(30))

  pair <- synth_pair(original = original, synthetic = synthetic,
                     holdout = holdout)
  result <- linkability(pair, n_attacks = 30, seed = 42)

  expect_equal(result$n_original, 60)
  expect_equal(result$n_holdout, 30)
})

# --- Input validation ---

test_that("linkability input validation works", {
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20), x3 = rnorm(20))
  Y <- data.frame(x1 = rnorm(10), x2 = rnorm(10), x3 = rnorm(10))

  expect_error(linkability(1:10, Y), "X must be a data frame")
  expect_error(linkability(X, 1:10), "Y must be a data frame")
  expect_error(linkability(X, Y, holdout = "not a df"),
               "holdout must be a data frame")
})

test_that("linkability requires at least 2 variables", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = 1:5)

  expect_error(linkability(X, Y), "At least 2 variables")
})

test_that("linkability errors when no common variables", {
  X <- data.frame(a = 1:10, b = 1:10)
  Y <- data.frame(c = 1:10, d = 1:10)

  expect_error(linkability(X, Y), "At least 2 variables")
})

test_that("linkability errors with mismatched variable types", {
  X <- data.frame(x = 1:10, y = 1:10)
  Y <- data.frame(x = letters[1:10], y = 1:10)

  expect_error(linkability(X, Y), "different class")
})

test_that("linkability errors with missing variables", {
  X <- data.frame(a = 1:10, b = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10, c = 1:10)

  expect_error(linkability(X, Y, vars = c("a", "b", "c")),
               "Variables missing in X")
})

test_that("linkability validates n_attacks", {
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  Y <- data.frame(x1 = rnorm(10), x2 = rnorm(10))

  expect_error(linkability(X, Y, n_attacks = -1), "positive integer")
})

test_that("linkability validates n_neighbors", {
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  Y <- data.frame(x1 = rnorm(10), x2 = rnorm(10))

  expect_error(linkability(X, Y, n_neighbors = 0), "positive integer")
})

test_that("linkability validates confidence_level", {
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  Y <- data.frame(x1 = rnorm(10), x2 = rnorm(10))

  expect_error(linkability(X, Y, confidence_level = 1.5),
               "between 0 and 1")
})

test_that("linkability handles NA values", {
  set.seed(123)
  X <- data.frame(x1 = c(rnorm(18), NA, NA), x2 = rnorm(20), x3 = rnorm(20))
  Y <- data.frame(x1 = c(rnorm(9), NA), x2 = rnorm(10), x3 = rnorm(10))

  result <- linkability(X, Y, n_attacks = 20, seed = 42, na.rm = TRUE)

  expect_s3_class(result, "linkability")
})

# --- n_attacks capping ---

test_that("n_attacks is capped at nrow(Y)", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))
  Y <- data.frame(x1 = rnorm(10), x2 = rnorm(10), x3 = rnorm(10))

  result <- linkability(X, Y, n_attacks = 1000, seed = 42)

  expect_equal(result$n_attacks, 10)
  expect_equal(length(result$links), 10)
})

# --- Variable selection ---

test_that("linkability uses specified vars", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60), x4 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40), x4 = rnorm(40))

  result <- linkability(X, Y, vars = c("x1", "x2", "x3"),
                         n_attacks = 20, seed = 42)

  expect_equal(sort(result$vars), c("x1", "x2", "x3"))
})

test_that("linkability auto-selects common vars", {
  set.seed(123)
  X <- data.frame(a = rnorm(60), b = rnorm(60), c_only = rnorm(60))
  Y <- data.frame(a = rnorm(40), b = rnorm(40), d_only = rnorm(40))

  result <- linkability(X, Y, n_attacks = 20, seed = 42)

  expect_equal(sort(result$vars), c("a", "b"))
})

# --- Mixed data types ---

test_that("linkability handles mixed data types", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    age = sample(20:70, n, replace = TRUE),
    income = rnorm(n, 50000, 15000),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE))
  )
  Y <- data.frame(
    age = sample(20:70, n, replace = TRUE),
    income = rnorm(n, 50000, 15000),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE))
  )

  result <- linkability(X, Y, n_attacks = 50, seed = 42)

  expect_s3_class(result, "linkability")
})

# --- Dataset sizes stored correctly ---

test_that("linkability stores correct dataset sizes", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
  Y <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))

  result <- linkability(X, Y, holdout_fraction = 0.4,
                         n_attacks = 50, seed = 42)

  expect_equal(result$n_synthetic, 80)
  expect_equal(result$n_holdout, 40)
  expect_equal(result$n_original, 60)
})

# --- print, summary, plot ---

test_that("linkability print method works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))

  result <- linkability(X, Y, n_attacks = 30, seed = 42)

  expect_output(print(result), "Linkability Risk Assessment")
  expect_output(print(result), "Attack Results")
  expect_output(print(result), "Risk Score")
  expect_output(print(result), "Privacy Assessment")
})

test_that("linkability summary returns correct class", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))

  result <- linkability(X, Y, n_attacks = 30, seed = 42)
  s <- summary(result)

  expect_s3_class(s, "summary.linkability")
  expect_true(!is.null(s$risk))
  expect_true(!is.null(s$aux_cols))
  expect_true(!is.null(s$secret_cols))
})

test_that("linkability print.summary works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))

  result <- linkability(X, Y, n_attacks = 30, seed = 42)
  s <- summary(result)

  expect_output(print(s), "Summary: Linkability Risk")
  expect_output(print(s), "Key Metrics")
  expect_output(print(s), "Link Results")
})

test_that("linkability plot method works without error", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))

  result <- linkability(X, Y, n_attacks = 30, seed = 42)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

# --- Consistency checks ---

test_that("n_success matches sum of links", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))

  result <- linkability(X, Y, n_attacks = 50, seed = 42)

  expect_equal(result$n_success, sum(result$links))
  expect_equal(result$n_control_success, sum(result$links_control))
})

test_that("risk_attack equals n_success / n_attacks", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))

  result <- linkability(X, Y, n_attacks = 50, seed = 42)

  expect_equal(result$risk_attack, result$n_success / result$n_attacks)
})

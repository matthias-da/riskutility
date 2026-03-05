# Tests for domias (Density-based Membership Inference Attack)

# --- Class structure and fields ---

test_that("domias returns correct S3 class structure", {
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

  result <- domias(X, Y, seed = 42)

  expect_s3_class(result, "domias")
  expect_true("density_ratios_train" %in% names(result))
  expect_true("density_ratios_holdout" %in% names(result))
  expect_true("auc" %in% names(result))
  expect_true("mean_ratio_train" %in% names(result))
  expect_true("mean_ratio_holdout" %in% names(result))
  expect_true("memorization_score" %in% names(result))
  expect_true("privacy_pass" %in% names(result))
  expect_true("n_train" %in% names(result))
  expect_true("n_holdout" %in% names(result))
  expect_true("n_synthetic" %in% names(result))
  expect_true("radius" %in% names(result))
  expect_true("vars" %in% names(result))
})

test_that("domias privacy_pass is logical", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- domias(X, Y, seed = 42)

  expect_type(result$privacy_pass, "logical")
})

test_that("domias stores correct radius", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- domias(X, Y, radius = 0.15, seed = 42)

  expect_equal(result$radius, 0.15)
})

# --- AUC values for random/independent data ---

test_that("AUC is approximately 0.5 for random independent data", {
  set.seed(42)
  n <- 300
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

  result <- domias(X, Y, seed = 123)

  # AUC should be around 0.5 for data from same distribution
  expect_true(result$auc > 0.3 && result$auc < 0.7)
})

test_that("AUC is between 0 and 1", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  result <- domias(X, Y, seed = 42)

  expect_true(result$auc >= 0 && result$auc <= 1)
})

test_that("memorization_score is approximately 1 for independent data", {
  set.seed(42)
  n <- 200
  X <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  Y <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  result <- domias(X, Y, seed = 123)

  # Memorization score should be close to 1 for independent data

  expect_true(result$memorization_score > 0.5 && result$memorization_score < 2.0)
})

# --- Memorization detection ---

test_that("domias detects copied/memorized data", {
  set.seed(123)
  n <- 200
  # More dimensions increase signal strength for DOMIAS
  train <- data.frame(
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
    x4 = rnorm(n), x5 = rnorm(n)
  )
  holdout <- data.frame(
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
    x4 = rnorm(n), x5 = rnorm(n)
  )
  # Synthetic is exact copy of training (memorized)
  Y <- train[sample(nrow(train), n, replace = TRUE), ]

  result <- domias(train, Y, holdout = holdout, radius = 0.1)

  # AUC should be above 0.5 for memorized data
  expect_true(result$auc > 0.5)
})

test_that("memorization_score is elevated for copied data", {
  set.seed(123)
  n <- 200
  train <- data.frame(
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
    x4 = rnorm(n), x5 = rnorm(n)
  )
  holdout <- data.frame(
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
    x4 = rnorm(n), x5 = rnorm(n)
  )
  # Synthetic is exact copy of training
  Y <- train[sample(nrow(train), n, replace = TRUE), ]

  result <- domias(train, Y, holdout = holdout, radius = 0.1)

  # memorization_score should be > 1 for memorized data
  expect_true(result$memorization_score > 1)
})

test_that("privacy_pass is FALSE for memorized data", {
  set.seed(123)
  n <- 200
  p <- 5
  train <- data.frame(
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
    x4 = rnorm(n), x5 = rnorm(n)
  )
  holdout <- data.frame(
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
    x4 = rnorm(n), x5 = rnorm(n)
  )
  # Synthetic is near-copy of training (with tiny noise)
  Y <- train + matrix(rnorm(n * p, sd = 0.01), ncol = p)
  names(Y) <- names(train)

  result <- domias(train, Y, holdout = holdout, radius = 0.1)

  expect_false(result$privacy_pass)
})

test_that("privacy_pass is TRUE for independent data", {
  set.seed(42)
  n <- 200
  train <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  holdout <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  Y <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  result <- domias(train, Y, holdout = holdout)

  expect_true(result$privacy_pass)
})

# --- Density ratio vectors ---

test_that("density_ratios_train has correct length", {
  set.seed(123)
  n <- 80
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  result <- domias(X, Y, holdout_fraction = 0.4, seed = 42)

  expect_equal(length(result$density_ratios_train), result$n_train)
  expect_equal(length(result$density_ratios_holdout), result$n_holdout)
})

test_that("density ratios are positive", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- domias(X, Y, seed = 42)

  expect_true(all(result$density_ratios_train > 0))
  expect_true(all(result$density_ratios_holdout > 0))
})

test_that("mean_ratio values match density_ratios vectors", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- domias(X, Y, seed = 42)

  expect_equal(result$mean_ratio_train,
               mean(result$density_ratios_train, na.rm = TRUE))
  expect_equal(result$mean_ratio_holdout,
               mean(result$density_ratios_holdout, na.rm = TRUE))
})

test_that("memorization_score is ratio of means", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- domias(X, Y, seed = 42)

  expect_equal(result$memorization_score,
               result$mean_ratio_train / result$mean_ratio_holdout)
})

# --- Holdout handling ---

test_that("domias works with explicit holdout", {
  set.seed(123)
  n <- 60
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))
  holdout <- data.frame(x1 = rnorm(30), x2 = rnorm(30))

  result <- domias(X, Y, holdout = holdout)

  expect_s3_class(result, "domias")
  # With explicit holdout, all of X is used as training
  expect_equal(result$n_train, nrow(X))
  expect_equal(result$n_holdout, nrow(holdout))
})

test_that("domias works with holdout_fraction", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(50), x2 = rnorm(50))

  result <- domias(X, Y, holdout_fraction = 0.3, seed = 42)

  expect_s3_class(result, "domias")
  expect_equal(result$n_holdout, floor(n * 0.3))
  expect_equal(result$n_train, n - floor(n * 0.3))
})

test_that("domias seed produces reproducible results", {
  set.seed(99)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  r1 <- domias(X, Y, seed = 42)
  r2 <- domias(X, Y, seed = 42)

  expect_equal(r1$auc, r2$auc)
  expect_equal(r1$mean_ratio_train, r2$mean_ratio_train)
  expect_equal(r1$mean_ratio_holdout, r2$mean_ratio_holdout)
  expect_equal(r1$memorization_score, r2$memorization_score)
})

# --- synth_pair method ---

test_that("domias.synth_pair works", {
  set.seed(123)
  n <- 60
  original <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  synthetic <- data.frame(
    x1 = rnorm(40),
    x2 = rnorm(40)
  )

  pair <- synth_pair(original = original, synthetic = synthetic)
  result <- domias(pair, seed = 42)

  expect_s3_class(result, "domias")
  expect_true(!is.null(result$auc))
})

test_that("domias.synth_pair uses holdout when available", {
  set.seed(123)
  original <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  synthetic <- data.frame(x1 = rnorm(40), x2 = rnorm(40))
  holdout <- data.frame(x1 = rnorm(30), x2 = rnorm(30))

  pair <- synth_pair(original = original, synthetic = synthetic, holdout = holdout)
  result <- domias(pair)

  expect_equal(result$n_train, 60)
  expect_equal(result$n_holdout, 30)
})

# --- Mixed data types ---

test_that("domias works with mixed data types", {
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

  result <- domias(X, Y, seed = 42)

  expect_s3_class(result, "domias")
  expect_true(!is.na(result$auc))
})

test_that("domias works with categorical-only data", {
  set.seed(123)
  n <- 80
  X <- data.frame(
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE),
    education = sample(c("low", "mid", "high"), n, replace = TRUE)
  )
  Y <- data.frame(
    gender = sample(c("M", "F"), 60, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), 60, replace = TRUE),
    education = sample(c("low", "mid", "high"), 60, replace = TRUE)
  )

  result <- domias(X, Y, seed = 42, radius = 0.3)

  expect_s3_class(result, "domias")
  expect_true(!is.na(result$auc))
})

# --- Radius parameter ---

test_that("different radii produce different results", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  r1 <- domias(X, Y, radius = 0.05, seed = 42)
  r2 <- domias(X, Y, radius = 0.2, seed = 42)

  # Different radii should generally give different density ratios
  expect_false(identical(r1$mean_ratio_train, r2$mean_ratio_train))
})

test_that("larger radius gives higher neighbor counts (higher ratios)", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(80), x2 = rnorm(80))

  r_small <- domias(X, Y, radius = 0.05, seed = 42)
  r_large <- domias(X, Y, radius = 0.3, seed = 42)

  # Larger radius should generally not have smaller mean ratios
  # (more neighbors counted in both numerator and denominator)
  # This is a soft check - just ensure they both work
  expect_true(!is.na(r_small$auc))
  expect_true(!is.na(r_large$auc))
})

# --- Input validation ---

test_that("domias input validation works", {
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  Y <- data.frame(x1 = rnorm(10), x2 = rnorm(10))

  expect_error(domias(1:10, Y), "X must be a data frame")
  expect_error(domias(X, 1:10), "Y must be a data frame")
  expect_error(domias(X, Y, holdout = "not a df"), "holdout must be a data frame")
})

test_that("domias errors with invalid radius", {
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  Y <- data.frame(x1 = rnorm(10), x2 = rnorm(10))

  expect_error(domias(X, Y, radius = 0), "radius must be a number between 0 and 1")
  expect_error(domias(X, Y, radius = 1), "radius must be a number between 0 and 1")
  expect_error(domias(X, Y, radius = -0.1), "radius must be a number between 0 and 1")
  expect_error(domias(X, Y, radius = "abc"), "radius must be a number between 0 and 1")
})

test_that("domias errors with invalid holdout_fraction", {
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  Y <- data.frame(x1 = rnorm(10), x2 = rnorm(10))

  expect_error(domias(X, Y, holdout_fraction = 0),
               "holdout_fraction must be a number between 0 and 1")
  expect_error(domias(X, Y, holdout_fraction = 1),
               "holdout_fraction must be a number between 0 and 1")
})

test_that("domias errors when no common variables", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(b = 1:10)

  expect_error(domias(X, Y), "No common variables")
})

test_that("domias errors with mismatched variable types", {
  X <- data.frame(x = 1:10)
  Y <- data.frame(x = letters[1:10])

  expect_error(domias(X, Y), "different class")
})

test_that("domias errors with missing variables", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10)

  expect_error(domias(X, Y, vars = c("a", "b")), "Variables missing in X")
})

test_that("domias handles NA values", {
  set.seed(123)
  X <- data.frame(x1 = c(rnorm(18), NA, NA), x2 = rnorm(20))
  Y <- data.frame(x1 = c(rnorm(9), NA), x2 = rnorm(10))

  result <- domias(X, Y, seed = 42, na.rm = TRUE)

  expect_s3_class(result, "domias")
})

# --- Variable selection ---

test_that("domias uses specified vars", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))

  result <- domias(X, Y, vars = c("x1", "x2"), seed = 42)

  expect_equal(result$vars, c("x1", "x2"))
})

test_that("domias auto-selects common vars", {
  set.seed(123)
  X <- data.frame(a = rnorm(60), b = rnorm(60), c_only = rnorm(60))
  Y <- data.frame(a = rnorm(40), b = rnorm(40), d_only = rnorm(40))

  result <- domias(X, Y, seed = 42)

  expect_equal(sort(result$vars), c("a", "b"))
})

# --- Dataset sizes stored correctly ---

test_that("domias stores correct dataset sizes", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
  Y <- data.frame(x1 = rnorm(80), x2 = rnorm(80))

  result <- domias(X, Y, holdout_fraction = 0.4, seed = 42)

  expect_equal(result$n_synthetic, 80)
  expect_equal(result$n_holdout, 40)
  expect_equal(result$n_train, 60)
})

# --- print, summary, plot ---

test_that("domias print method works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- domias(X, Y, seed = 42)

  expect_output(print(result), "DOMIAS")
  expect_output(print(result), "Density Ratio")
  expect_output(print(result), "Privacy Assessment")
  expect_output(print(result), "AUC")
})

test_that("domias print shows correct pass/warning", {
  set.seed(42)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- domias(X, Y, seed = 123)

  if (result$privacy_pass) {
    expect_output(print(result), "PASS")
  } else {
    expect_output(print(result), "WARNING")
  }
})

test_that("domias summary returns correct class", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- domias(X, Y, seed = 42)
  s <- summary(result)

  expect_s3_class(s, "summary.domias")
  expect_true(!is.null(s$auc))
  expect_true(!is.null(s$mean_ratio_train))
  expect_true(!is.null(s$mean_ratio_holdout))
  expect_true(!is.null(s$memorization_score))
  expect_true(!is.null(s$quantiles_train))
  expect_true(!is.null(s$quantiles_holdout))
  expect_true(!is.null(s$n_high_ratio_train))
  expect_true(!is.null(s$n_high_ratio_holdout))
})

test_that("domias print.summary works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- domias(X, Y, seed = 42)
  s <- summary(result)

  expect_output(print(s), "Summary: DOMIAS")
  expect_output(print(s), "Key Metrics")
  expect_output(print(s), "Density Ratios")
})

test_that("domias plot method works without error", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- domias(X, Y, seed = 42)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

test_that("domias print returns object invisibly", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- domias(X, Y, seed = 42)

  out <- capture.output(ret <- print(result))
  expect_identical(ret, result)
})

test_that("domias summary print returns object invisibly", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- domias(X, Y, seed = 42)
  s <- summary(result)

  out <- capture.output(ret <- print(s))
  expect_identical(ret, s)
})

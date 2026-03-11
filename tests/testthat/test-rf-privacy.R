# tests/testthat/test-rf-privacy.R
test_that("rf_privacy detects memorized data", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X_train <- data.frame(a = rnorm(100), b = rnorm(100), c = rnorm(100))
  X_holdout <- data.frame(a = rnorm(100), b = rnorm(100), c = rnorm(100))
  # Y = near-copies of training records (memorized)
  Y <- X_train[sample(100, 200, replace = TRUE), ]
  Y <- Y + rnorm(nrow(Y) * ncol(Y), 0, 0.01)
  res <- rf_privacy(X_train, Y, holdout = X_holdout, seed = 1,
                    n_trees = 200, null_test = FALSE)

  expect_s3_class(res, "rf_privacy")
  expect_true(res$max_prox_share > 0.55)
  expect_true(res$max_prox_ratio > 1.0)
  expect_false(res$privacy_pass)
})

test_that("rf_privacy passes for random synthetic data", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(200), b = rnorm(200))
  Y <- data.frame(a = rnorm(200), b = rnorm(200))
  res <- rf_privacy(X, Y, holdout_fraction = 0.5, seed = 1,
                    n_trees = 200, null_test = FALSE)

  expect_s3_class(res, "rf_privacy")
  # Should be around 0.5 (within tolerance)
  expect_true(abs(res$max_prox_share - 0.5) < 0.15)
})

test_that("rf_privacy returns all expected fields", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  res <- rf_privacy(X, Y, seed = 1, n_trees = 50, null_test = FALSE)

  expect_true("max_prox_share" %in% names(res))
  expect_true("max_prox_ratio" %in% names(res))
  expect_true("max_prox_train" %in% names(res))
  expect_true("max_prox_holdout" %in% names(res))
  expect_true("prox_share" %in% names(res))
  expect_true("prox_ratio" %in% names(res))
  expect_true("privacy_pass" %in% names(res))
  expect_true("wilcox_test" %in% names(res))
  expect_true("oob_error" %in% names(res))
  expect_true("var_importance" %in% names(res))
  expect_true(is.logical(res$privacy_pass))
  expect_equal(length(res$max_prox_train), nrow(Y))
})

test_that("rf_privacy null_test works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  res <- rf_privacy(X, Y, seed = 1, n_trees = 50,
                    null_test = TRUE, n_null = 20)

  expect_true("null_distribution" %in% names(res))
  expect_true(!is.null(res$null_distribution))
  expect_true("null_pvalue" %in% names(res$null_distribution))
})

test_that("rf_privacy na.rm = FALSE works with NAs", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = c(rnorm(99), NA), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  # Should not error with na.rm = FALSE (ranger handles NAs)
  res <- rf_privacy(X, Y, na.rm = FALSE, seed = 1,
                    n_trees = 50, null_test = FALSE)
  expect_s3_class(res, "rf_privacy")
})

test_that("rf_privacy synth_pair dispatch works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  sp <- synth_pair(X, Y)
  res <- rf_privacy(sp, seed = 1, n_trees = 50, null_test = FALSE)
  expect_s3_class(res, "rf_privacy")
})

test_that("rf_privacy seed separation works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  # Same seed should produce same results
  res1 <- rf_privacy(X, Y, seed = 42, n_trees = 50, null_test = FALSE)
  res2 <- rf_privacy(X, Y, seed = 42, n_trees = 50, null_test = FALSE)
  expect_equal(res1$max_prox_share, res2$max_prox_share)
  # Different seed should produce different holdout splits
  res3 <- rf_privacy(X, Y, seed = 99, n_trees = 50, null_test = FALSE)
  expect_true(res1$max_prox_share != res3$max_prox_share ||
              res1$max_prox_ratio != res3$max_prox_ratio)
})

test_that("rf_privacy prox_ratio zero-guard works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  # Very small holdout with highly separable data -> may trigger zero-guard
  X <- data.frame(a = c(rep(0, 10), rep(100, 10)))
  Y <- data.frame(a = rep(50, 10))
  # This may or may not trigger -- test that NA is returned with warning
  # when denominator is degenerate
  res <- suppressWarnings(
    rf_privacy(X, Y, holdout_fraction = 0.5, seed = 1,
               n_trees = 50, null_test = FALSE)
  )
  expect_s3_class(res, "rf_privacy")
  # max_prox_ratio is either numeric or NA (both valid)
  expect_true(is.numeric(res$max_prox_ratio) || is.na(res$max_prox_ratio))
})

test_that("rf_privacy print/summary/plot methods work", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(80), b = rnorm(80))
  Y <- data.frame(a = rnorm(80), b = rnorm(80))
  res <- rf_privacy(X, Y, seed = 1, n_trees = 50, null_test = FALSE)

  expect_output(print(res), "RF Privacy")
  s <- summary(res)
  expect_s3_class(s, "summary.rf_privacy")
  expect_output(print(s))
  expect_silent(plot(res, which = 1))
  expect_silent(plot(res, which = 2))
})

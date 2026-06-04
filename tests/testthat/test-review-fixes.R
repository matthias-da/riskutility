# Regression tests for the correctness fixes from the 2026-06 full-suite review.

test_that("gower() returns distance 0 / utility 1 for identical data", {
  X <- data.frame(
    a = c(1, 2, 3, 4),
    b = factor(c("x", "y", "x", "y")),
    c = c(10, 20, 30, 40)
  )
  res <- gower(X, X)
  expect_s3_class(res, "gower")
  expect_equal(res$gower_distance, 0, tolerance = 1e-8)
  expect_equal(res$utility_score, 1, tolerance = 1e-8)
})

test_that("gower() stays bounded in [0,1] and requires equal nrow", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- gower(X, Y)
  expect_gte(res$gower_distance, 0)
  expect_lte(res$gower_distance, 1)
  expect_error(gower(X, Y[1:10, ]), "same number of rows")
})

test_that("KLDiv()/JSDiv() are finite with zero-probability entries", {
  p <- c(0, 0.5, 0.5)
  q <- c(0.2, 0.3, 0.5)
  expect_true(is.finite(KLDiv(p, q)))           # 0 * log(0/q) handled as 0
  expect_equal(KLDiv(c(0.5, 0.5), c(0.5, 0.5)), 0, tolerance = 1e-12)
  expect_true(is.finite(JSDiv(p, q)))
  expect_gte(JSDiv(p, q), 0)
})

test_that("nndr() flags exact copies (NNDR = 0) instead of dropping them as NA", {
  set.seed(42)
  base <- data.frame(a = rnorm(30), b = rnorm(30))
  dup  <- base[1, , drop = FALSE]
  original <- rbind(base, dup, dup)             # row 1 now appears 3x (d1 = d2 = 0)
  holdout  <- data.frame(a = rnorm(30), b = rnorm(30))
  synth    <- rbind(data.frame(a = rnorm(20), b = rnorm(20)), dup)  # copies the record
  res <- nndr(original, synth, holdout = holdout, method = "euclidean")
  expect_true(is.finite(res$mean_nndr_train))
  expect_equal(min(res$nndr_train, na.rm = TRUE), 0, tolerance = 1e-12)
})

test_that("plot.denpca runs for both 'which' values (was a crash)", {
  skip_if_not_installed("ggplot2")
  set.seed(7)
  X <- data.frame(v1 = rnorm(80), v2 = rnorm(80, 1),
                  v3 = rnorm(80, 2), v4 = rnorm(80, 3))
  Y <- data.frame(v1 = rnorm(80, .2), v2 = rnorm(80, 1.1),
                  v3 = rnorm(80, 2.1), v4 = rnorm(80, 3.1))
  d <- densitydiff_pca(X, Y, bayesspace = FALSE)
  expect_s3_class(d, "denpca")
  expect_false(is.null(d$kl))
  expect_no_error(suppressWarnings(plot(d, which = 1)))
  expect_no_error(suppressWarnings(plot(d, which = 2)))
})

test_that("plot.denratio runs (smoke test)", {
  skip_if_not_installed("ggplot2")
  set.seed(8)
  d <- densitydiff_1d_num(rnorm(150, 40, 10), rnorm(150, 42, 11))
  expect_s3_class(d, "denratio")
  expect_no_error(suppressWarnings(plot(d)))
})

test_that("densitydiff_kl_num trivariate is finite on a shared grid", {
  skip_if_not_installed("misc3d")
  set.seed(9)
  X <- data.frame(a = rnorm(100), b = rnorm(100), c = rnorm(100))
  Y <- data.frame(a = rnorm(100, .1), b = rnorm(100, .1), c = rnorm(100, .1))
  val <- densitydiff_kl_num(X, Y)
  expect_true(is.finite(val))
})

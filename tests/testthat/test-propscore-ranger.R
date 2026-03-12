# tests/testthat/test-propscore-ranger.R
test_that("propscore(method = 'ranger') returns correct class", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  res <- propscore(X, Y, method = "ranger", proximity = "none")

  expect_s3_class(res, "propscore")
  expect_equal(res$method, "ranger")
  expect_true("oob_error" %in% names(res))
  expect_true("var_importance" %in% names(res))
})

test_that("propscore(method = 'ranger') identical data gives OOB ~0.5", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(200), b = rnorm(200))
  Y <- X  # identical
  res <- propscore(X, Y, method = "ranger", proximity = "none")

  expect_true(res$oob_error > 0.35)  # near 0.5
})

test_that("propscore(method = 'ranger') different data gives OOB < 0.3", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(200), b = rnorm(200))
  Y <- data.frame(a = rnorm(200) + 3, b = rnorm(200) + 3)
  res <- propscore(X, Y, method = "ranger", proximity = "none")

  expect_true(res$oob_error < 0.3)
})

test_that("propscore proximity = 'summary' returns structural metrics", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  res <- propscore(X, Y, method = "ranger", proximity = "summary")

  expect_true("within_orig_prox" %in% names(res))
  expect_true("within_synth_prox" %in% names(res))
  expect_true("cross_prox" %in% names(res))
  expect_true("structure_ratio" %in% names(res))
  expect_null(res$proximity_matrix)
})

test_that("propscore proximity = 'full' stores matrix", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  res <- propscore(X, Y, method = "ranger", proximity = "full")

  expect_true(!is.null(res$proximity_matrix))
  expect_equal(nrow(res$proximity_matrix), 60)
  expect_equal(ncol(res$proximity_matrix), 60)
})

test_that("propscore match.arg rejects invalid method", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = 1:10)
  expect_error(propscore(X, Y, method = "invalid"), "arg")
})

test_that("propscore existing methods still work after match.arg", {
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50))
  Y <- data.frame(a = rnorm(50), b = rnorm(50))

  res_rf <- propscore(X, Y, method = "rf")
  expect_s3_class(res_rf, "propscore")

  res_lr <- propscore(X, Y, method = "logreg")
  expect_s3_class(res_lr, "propscore")
})

test_that("propscore emits message for ranger params with non-ranger method", {
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50))
  Y <- data.frame(a = rnorm(50), b = rnorm(50))
  expect_message(
    propscore(X, Y, method = "rf", proximity = "summary"),
    "ranger"
  )
})

test_that("propscore(method = 'ranger') synth_pair dispatch works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50))
  Y <- data.frame(a = rnorm(50), b = rnorm(50))
  sp <- synth_pair(X, Y)
  res <- propscore(sp, method = "ranger", proximity = "summary")
  expect_s3_class(res, "propscore")
  expect_true("structure_ratio" %in% names(res))
})

test_that("propscore plot which = 3 and 4 work for ranger", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50))
  Y <- data.frame(a = rnorm(50), b = rnorm(50))
  res <- propscore(X, Y, method = "ranger", proximity = "summary")

  expect_silent(plot(res, which = 3))
  expect_silent(plot(res, which = 4))
})

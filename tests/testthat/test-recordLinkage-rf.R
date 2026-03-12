# tests/testthat/test-recordLinkage-rf.R
test_that("recordLinkage(method = 'rf') independent matching works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = sample(letters[1:5], 50, TRUE))
  Y <- data.frame(a = rnorm(50), b = sample(letters[1:5], 50, TRUE))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                       n_trees = 100, matching = "independent")

  expect_s3_class(res, "recordLinkageRisk")
  expect_true("var_importance" %in% names(res))
  expect_true(all(res$per_record$risk >= 0 & res$per_record$risk <= 1))
})

test_that("recordLinkage(method = 'rf') bijective matching works", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("clue")
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                       n_trees = 100, matching = "bijective")

  expect_s3_class(res, "recordLinkageRisk")
  # Bijective risk is binary
  expect_true(all(res$per_record$risk %in% c(0, 1)))
})

test_that("recordLinkage(method = 'rf') with blocking works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rep(c("A", "B"), 50))
  Y <- data.frame(a = rnorm(100), b = rep(c("A", "B"), 50))
  res <- recordLinkage(X, Y, key = "a", block = "b", method = "rf",
                       n_trees = 50)

  expect_s3_class(res, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'rf') rf_global = TRUE works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rep(c("A", "B"), 50))
  Y <- data.frame(a = rnorm(100), b = rep(c("A", "B"), 50))
  res <- recordLinkage(X, Y, key = "a", block = "b", method = "rf",
                       n_trees = 50, rf_global = TRUE)

  expect_s3_class(res, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'rf') direction = 'reverse' works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  res_fwd <- recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                           n_trees = 50, direction = "forward")
  res_rev <- recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                           n_trees = 50, direction = "reverse")

  expect_s3_class(res_fwd, "recordLinkageRisk")
  expect_s3_class(res_rev, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'rf') emits message for weights", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  expect_message(
    recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                  weights = c(a = 1, b = 2), n_trees = 50),
    "importance"
  )
})

test_that("recordLinkage(method = 'rf') small block fallback", {
  skip_if_not_installed("ranger")
  set.seed(1)
  # Block C has only 1 record per class -> should fall back
  X <- data.frame(a = rnorm(31), b = c(rep("A", 15), rep("B", 15), "C"))
  Y <- data.frame(a = rnorm(31), b = c(rep("A", 15), rep("B", 15), "C"))
  expect_message(
    res <- recordLinkage(X, Y, key = "a", block = "b",
                         method = "rf", n_trees = 50),
    "blocks"
  )
  expect_s3_class(res, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'rf') truth = 'row' works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- X[sample(20), ]  # permuted copy
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                       truth = "row", n_trees = 50)
  expect_s3_class(res, "recordLinkageRisk")
  expect_true("true_in_set" %in% names(res$per_record))
})

test_that("recordLinkage(method = 'rf') var_importance plot works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                       n_trees = 50)
  # which = 3 is existing var importance plot
  expect_silent(plot(res, which = 3))
})

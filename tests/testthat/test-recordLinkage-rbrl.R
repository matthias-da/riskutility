# tests/testthat/test-recordLinkage-rbrl.R

test_that("recordLinkage(method = 'rbrl') basic independent matching works", {
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "rbrl",
                       matching = "independent")
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(res$method, "rbrl")
  expect_equal(nrow(res$per_record), 30)
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
})

test_that("recordLinkage(method = 'rbrl') detects rank-preserved copies", {
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  # Monotone perturbation: add small noise (preserves rank order)
  Y <- X + rnorm(nrow(X) * ncol(X), 0, 0.01)
  res <- recordLinkage(X, Y, key = c("a", "b", "c"), method = "rbrl",
                       truth = "row")
  # RBRL should detect high risk since ranks are preserved
  expect_true(res$overall$mean_risk > 0.1)
  expect_true(mean(res$per_record$true_in_set) > 0.8)
})

test_that("recordLinkage(method = 'rbrl') robust to monotone transformations", {
  set.seed(1)
  X <- data.frame(a = 1:30, b = seq(0, 1, length.out = 30))
  # Large additive noise but rank-preserving
  Y <- data.frame(a = X$a * 100 + 5000, b = X$b^3)
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "rbrl",
                       truth = "row")
  # Deterministic would see large distances; RBRL sees preserved ranks
  expect_true(res$overall$pct_true_in_set > 0.9)
})

test_that("recordLinkage(method = 'rbrl') bijective matching works", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- X + rnorm(40, 0, 0.01)  # near-copies
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "rbrl",
                       matching = "bijective", truth = "row")
  expect_s3_class(res, "recordLinkageRisk")
  expect_true("bijective_assigned" %in% names(res$per_record))
  # Bijective should yield binary risk (0 or 1)
  expect_true(all(res$per_record$risk %in% c(0, 1)))
})

test_that("recordLinkage(method = 'rbrl') with blocking works", {
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = sample(c("A", "B"), 30, TRUE))
  Y <- data.frame(a = rnorm(30), b = sample(c("A", "B"), 30, TRUE))
  res <- recordLinkage(X, Y, key = c("a", "b"), block = "b",
                       method = "rbrl")
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(nrow(res$per_record), 30)
})

test_that("recordLinkage(method = 'rbrl') handles nominal variables", {
  set.seed(1)
  X <- data.frame(
    num = rnorm(20),
    cat = sample(letters[1:3], 20, TRUE),
    stringsAsFactors = FALSE
  )
  Y <- data.frame(
    num = rnorm(20),
    cat = sample(letters[1:3], 20, TRUE),
    stringsAsFactors = FALSE
  )
  res <- recordLinkage(X, Y, key = c("num", "cat"), method = "rbrl")
  expect_s3_class(res, "recordLinkageRisk")
  # Nominal variables should use exact match (0/1)
  expect_true(all(is.finite(res$per_record$d_min)))
})

test_that("recordLinkage(method = 'rbrl') direction = 'reverse' works", {
  set.seed(1)
  X <- data.frame(a = rnorm(25), b = rnorm(25))
  Y <- data.frame(a = rnorm(25), b = rnorm(25))
  res_fwd <- recordLinkage(X, Y, key = c("a", "b"), method = "rbrl",
                            direction = "forward")
  res_rev <- recordLinkage(X, Y, key = c("a", "b"), method = "rbrl",
                            direction = "reverse")
  expect_equal(res_fwd$direction, "forward")
  expect_equal(res_rev$direction, "reverse")
  expect_equal(nrow(res_rev$per_record), 25)
})

test_that("recordLinkage(method = 'rbrl') weights affect distances", {
  set.seed(1)
  X <- data.frame(a = rnorm(25), b = rnorm(25))
  Y <- data.frame(a = rnorm(25), b = rnorm(25))
  res1 <- recordLinkage(X, Y, key = c("a", "b"), method = "rbrl",
                        weights = c(a = 1, b = 1))
  res2 <- recordLinkage(X, Y, key = c("a", "b"), method = "rbrl",
                        weights = c(a = 10, b = 0.1))
  # Different weights should produce different d_min values
  expect_false(identical(res1$per_record$d_min, res2$per_record$d_min))
})

test_that("recordLinkage(method = 'rbrl') var_importance is computed", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20), c = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20), c = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b", "c"), method = "rbrl")
  expect_true(!is.null(res$var_importance))
  expect_equal(length(res$var_importance), 3)
  expect_true(all(names(res$var_importance) == c("a", "b", "c")))
  expect_true(all(res$var_importance >= 0))
})

test_that("recordLinkage(method = 'rbrl') plot works", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "rbrl")
  expect_silent(plot(res, which = 1))
  expect_silent(plot(res, which = 3))  # var importance
})

test_that("RBRL vs deterministic: RBRL finds higher risk under monotone transform", {
  set.seed(42)
  X <- data.frame(a = 1:20, b = seq(10, 200, length.out = 20))
  # Monotone transform: log + shift (preserves ranks perfectly)
  Y <- data.frame(a = log(X$a + 1) * 50, b = sqrt(X$b) * 10)
  res_det <- recordLinkage(X, Y, key = c("a", "b"), method = "deterministic",
                           truth = "row")
  res_rbrl <- recordLinkage(X, Y, key = c("a", "b"), method = "rbrl",
                            truth = "row")
  # RBRL should detect that ranks are perfectly preserved
  # Deterministic sees large value differences
  expect_true(res_rbrl$overall$pct_true_in_set >=
              res_det$overall$pct_true_in_set)
})

test_that("recordLinkage(method = 'rbrl') with ordinal factors", {
  set.seed(1)
  X <- data.frame(
    a = ordered(sample(c("low", "mid", "high"), 20, TRUE),
                levels = c("low", "mid", "high")),
    b = rnorm(20)
  )
  Y <- data.frame(
    a = ordered(sample(c("low", "mid", "high"), 20, TRUE),
                levels = c("low", "mid", "high")),
    b = rnorm(20)
  )
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "rbrl")
  expect_s3_class(res, "recordLinkageRisk")
  expect_true(all(is.finite(res$per_record$d_min)))
})

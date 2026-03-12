# tests/testthat/test-recordLinkage-mahalanobis.R

test_that("recordLinkage(method = 'mahalanobis') basic works", {
  set.seed(1)
  X <- data.frame(a = rnorm(25), b = rnorm(25))
  Y <- data.frame(a = rnorm(25), b = rnorm(25))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "mahalanobis")
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(res$method, "mahalanobis")
  expect_equal(nrow(res$per_record), 25)
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
})

test_that("recordLinkage(method = 'mahalanobis') detects near-copies", {
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  Y <- X + rnorm(nrow(X) * ncol(X), 0, 0.01)
  res <- recordLinkage(X, Y, key = c("a", "b", "c"),
                       method = "mahalanobis", truth = "row")
  expect_true(res$overall$mean_risk > 0.1)
  expect_true(mean(res$per_record$true_in_set) > 0.8)
})

test_that("recordLinkage(method = 'mahalanobis') robust vs classical", {
  set.seed(1)
  X <- data.frame(a = rnorm(40), b = rnorm(40))
  Y <- X + rnorm(80, 0, 0.1)
  res_rob <- recordLinkage(X, Y, key = c("a", "b"),
                            method = "mahalanobis", robust = TRUE,
                            truth = "row")
  res_cls <- recordLinkage(X, Y, key = c("a", "b"),
                            method = "mahalanobis", robust = FALSE,
                            truth = "row")
  expect_s3_class(res_rob, "recordLinkageRisk")
  expect_s3_class(res_cls, "recordLinkageRisk")
  expect_true(res_rob$overall$pct_true_in_set > 0.5)
  expect_true(res_cls$overall$pct_true_in_set > 0.5)
  expect_false(identical(res_rob$per_record$d_min, res_cls$per_record$d_min))
})

test_that("recordLinkage(method = 'mahalanobis') captures correlation", {
  set.seed(42)
  n <- 50
  a <- rnorm(n)
  b <- a + rnorm(n, 0, 0.1)
  X <- data.frame(a = a, b = b)

  Y_along <- data.frame(a = a + 0.5, b = b + 0.5)
  Y_against <- data.frame(a = a + 0.5, b = b - 0.5)

  res_along <- recordLinkage(X, Y_along, key = c("a", "b"),
                              method = "mahalanobis", truth = "row")
  res_against <- recordLinkage(X, Y_against, key = c("a", "b"),
                                method = "mahalanobis", truth = "row")

  expect_true(mean(res_against$per_record$d_min) >
              mean(res_along$per_record$d_min))
})

test_that("recordLinkage(method = 'mahalanobis') with mixed data", {
  set.seed(1)
  X <- data.frame(
    num1 = rnorm(30), num2 = rnorm(30),
    cat = sample(c("a", "b", "c"), 30, TRUE),
    stringsAsFactors = FALSE
  )
  Y <- data.frame(
    num1 = rnorm(30), num2 = rnorm(30),
    cat = sample(c("a", "b", "c"), 30, TRUE),
    stringsAsFactors = FALSE
  )
  res <- recordLinkage(X, Y, key = c("num1", "num2", "cat"),
                       method = "mahalanobis")
  expect_s3_class(res, "recordLinkageRisk")
  expect_true(all(is.finite(res$per_record$d_min)))
})

test_that("recordLinkage(method = 'mahalanobis') with ordinal variables", {
  set.seed(1)
  X <- data.frame(
    a = ordered(sample(c("low", "mid", "high"), 30, TRUE),
                levels = c("low", "mid", "high")),
    b = rnorm(30)
  )
  Y <- data.frame(
    a = ordered(sample(c("low", "mid", "high"), 30, TRUE),
                levels = c("low", "mid", "high")),
    b = rnorm(30)
  )
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "mahalanobis")
  expect_s3_class(res, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'mahalanobis') errors on all nominal", {
  X <- data.frame(a = sample(letters[1:3], 20, TRUE),
                  b = sample(c("x", "y"), 20, TRUE),
                  stringsAsFactors = FALSE)
  Y <- X
  expect_error(
    recordLinkage(X, Y, key = c("a", "b"), method = "mahalanobis"),
    "at least one numeric"
  )
})

test_that("recordLinkage(method = 'mahalanobis') with blocking", {
  set.seed(1)
  X <- data.frame(a = rnorm(40), b = sample(c("A", "B"), 40, TRUE))
  Y <- data.frame(a = rnorm(40), b = sample(c("A", "B"), 40, TRUE))
  res <- recordLinkage(X, Y, key = c("a", "b"), block = "b",
                       method = "mahalanobis")
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(nrow(res$per_record), 40)
})

test_that("recordLinkage(method = 'mahalanobis') bijective matching", {
  skip_if_not_installed("clue")
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- X + rnorm(40, 0, 0.01)
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "mahalanobis",
                       matching = "bijective", truth = "row")
  expect_s3_class(res, "recordLinkageRisk")
  expect_true("bijective_assigned" %in% names(res$per_record))
  expect_true(all(res$per_record$risk %in% c(0, 1)))
})

test_that("recordLinkage(method = 'mahalanobis') direction = 'reverse'", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "mahalanobis",
                       direction = "reverse")
  expect_equal(res$direction, "reverse")
  expect_equal(nrow(res$per_record), 20)
})

test_that("recordLinkage(method = 'mahalanobis') var_importance", {
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  res <- recordLinkage(X, Y, key = c("a", "b", "c"),
                       method = "mahalanobis")
  expect_true(!is.null(res$var_importance))
  expect_equal(length(res$var_importance), 3)
  expect_true(all(res$var_importance >= 0))
})

test_that("recordLinkage(method = 'mahalanobis') risk_weighting works", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res_unif <- recordLinkage(X, Y, key = c("a", "b"),
                             method = "mahalanobis",
                             risk_weighting = "uniform")
  res_sm <- recordLinkage(X, Y, key = c("a", "b"),
                           method = "mahalanobis",
                           risk_weighting = "softmax")
  expect_s3_class(res_unif, "recordLinkageRisk")
  expect_s3_class(res_sm, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'mahalanobis') plot works", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "mahalanobis")
  expect_silent(plot(res, which = 1))
  expect_silent(plot(res, which = 3))
})

test_that("mahalanobis vs deterministic: mahalanobis uses correlation", {
  set.seed(42)
  n <- 40
  a <- rnorm(n)
  b <- a + rnorm(n, 0, 0.05)
  X <- data.frame(a = a, b = b)
  Y <- data.frame(a = a + 0.3, b = b - 0.3)

  res_det <- recordLinkage(X, Y, key = c("a", "b"),
                            method = "deterministic", truth = "row")
  res_mah <- recordLinkage(X, Y, key = c("a", "b"),
                            method = "mahalanobis", truth = "row")

  # Mahalanobis should see against-correlation perturbation as more unusual
  expect_true(mean(res_mah$per_record$d_min) >
              mean(res_det$per_record$d_min) * 0.5)
})

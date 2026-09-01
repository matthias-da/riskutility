# Test: permutation importance in compare_feature_importance()
#
# The permutation backend is implemented internally (.permutation_importance);
# it must not require the vip package (archived from CRAN 2026-07-08) and must
# accept the caret-style metric names the documentation lists ("Accuracy",
# "RMSE", ...), plus a custom metric function.

make_class_data <- function(seed) {
  set.seed(seed)
  n <- 200
  inf <- rnorm(n)
  noise <- rnorm(n)
  outcome <- factor(ifelse(inf + rnorm(n, 0, 0.3) > 0, "A", "B"))
  data.frame(inf = inf, noise = noise, outcome = outcome)
}

make_reg_data <- function(seed) {
  set.seed(seed)
  n <- 200
  inf <- rnorm(n)
  noise <- rnorm(n)
  outcome <- 2 * inf + rnorm(n, 0, 0.5)
  data.frame(inf = inf, noise = noise, outcome = outcome)
}

test_that("permutation importance works with caret-style metric names, without vip", {
  skip_if_not_installed("caret")
  skip_if_not_installed("rpart")

  X <- make_class_data(1)
  Y <- make_class_data(2)
  pred_fn <- function(object, newdata) predict(object, newdata, type = "raw")

  set.seed(42)
  res <- compare_feature_importance(X, Y, outcome ~ .,
                                    method = "rpart",
                                    importance_type = "permutation",
                                    pred_wrapper = pred_fn,
                                    metric = "Accuracy", nsim = 3)

  expect_s3_class(res, "compare_feature_importance")
  imp <- setNames(res$comparison$importance_X, res$comparison$feature)
  # permuting the informative feature must hurt accuracy far more than noise
  expect_gt(imp["inf"], imp["noise"])
  expect_gt(imp["inf"], 0)
})

test_that("permutation importance handles error metrics with the correct sign", {
  skip_if_not_installed("caret")

  X <- make_reg_data(3)
  Y <- make_reg_data(4)
  pred_fn <- function(object, newdata) predict(object, newdata)

  set.seed(42)
  res <- compare_feature_importance(X, Y, outcome ~ .,
                                    method = "lm",
                                    importance_type = "permutation",
                                    pred_wrapper = pred_fn,
                                    metric = "RMSE", nsim = 3)

  imp <- setNames(res$comparison$importance_X, res$comparison$feature)
  # RMSE is smaller-is-better: permuting the informative feature must
  # INCREASE the error, giving a large positive importance
  expect_gt(imp["inf"], imp["noise"])
  expect_gt(imp["inf"], 0)
})

test_that("permutation importance accepts a custom metric function", {
  skip_if_not_installed("caret")
  skip_if_not_installed("rpart")

  X <- make_class_data(5)
  Y <- make_class_data(6)
  pred_fn <- function(object, newdata) predict(object, newdata, type = "raw")
  acc <- function(obs, pred) mean(pred == obs)

  set.seed(42)
  res <- compare_feature_importance(X, Y, outcome ~ .,
                                    method = "rpart",
                                    importance_type = "permutation",
                                    pred_wrapper = pred_fn,
                                    metric = acc, nsim = 3)

  imp <- setNames(res$comparison$importance_X, res$comparison$feature)
  expect_gt(imp["inf"], imp["noise"])
})

test_that("unknown metric names error informatively", {
  skip_if_not_installed("caret")
  skip_if_not_installed("rpart")

  X <- make_class_data(7)
  Y <- make_class_data(8)
  pred_fn <- function(object, newdata) predict(object, newdata, type = "raw")

  expect_error(
    compare_feature_importance(X, Y, outcome ~ .,
                               method = "rpart",
                               importance_type = "permutation",
                               pred_wrapper = pred_fn,
                               metric = "NotAMetric", nsim = 2),
    "metric"
  )
})

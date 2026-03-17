test_that("tstr returns correct S3 class with expected fields", {
  skip_if_not_installed("caret")
  set.seed(1)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- 2 * X$x1 + 0.5 * X$x2 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- 2 * Y$x1 + 0.5 * Y$x2 + rnorm(n, sd = 0.5)

  res <- tstr(X, Y, target_var = "y", seed = 1)

  expect_s3_class(res, "tstr")
  expect_true("tstr_ratio" %in% names(res))
  expect_true("utility_score" %in% names(res))
  expect_true("tstr_performance" %in% names(res))
  expect_true("trtr_performance" %in% names(res))
  expect_true("metric" %in% names(res))
  expect_true("model" %in% names(res))
  expect_true("target_var" %in% names(res))
  expect_true("n_train_orig" %in% names(res))
  expect_true("n_test" %in% names(res))
  expect_true("n_synth" %in% names(res))
  expect_true("test_fraction" %in% names(res))
  expect_true("per_target" %in% names(res))
})

test_that("tstr ratio is near 1 for identical distributions (regression)", {
  skip_if_not_installed("caret")
  set.seed(2)
  n <- 300
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- 3 * X$x1 - X$x2 + rnorm(n, sd = 0.3)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- 3 * Y$x1 - Y$x2 + rnorm(n, sd = 0.3)

  res <- tstr(X, Y, target_var = "y", seed = 2)

  expect_equal(res$metric, "R2")
  # Ratio should be close to 1 for identical generating process

  expect_true(res$tstr_ratio > 0.7)
  expect_true(res$tstr_ratio < 1.3)
})

test_that("tstr detects poor synthetic data (low ratio)", {
  skip_if_not_installed("caret")
  set.seed(3)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- 5 * X$x1 + 3 * X$x2 + rnorm(n, sd = 0.5)

  # Poor synthetic: relationship is destroyed (random noise for y)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- rnorm(n, sd = 10)

  res <- tstr(X, Y, target_var = "y", seed = 3)

  # Poor synthetic should have much lower ratio
  expect_true(res$tstr_ratio < 0.5)
})

test_that("utility_score is in [0,1]", {
  skip_if_not_installed("caret")
  set.seed(4)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- 2 * X$x1 + rnorm(n, sd = 0.5)

  # Good synthetic
  Y_good <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y_good$y <- 2 * Y_good$x1 + rnorm(n, sd = 0.5)
  res1 <- tstr(X, Y_good, target_var = "y", seed = 4)
  expect_true(res1$utility_score >= 0 && res1$utility_score <= 1)

  # Poor synthetic
  Y_poor <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y_poor$y <- rnorm(n, sd = 10)
  res2 <- tstr(X, Y_poor, target_var = "y", seed = 4)
  expect_true(res2$utility_score >= 0 && res2$utility_score <= 1)
})

test_that("classification target uses AUC", {
  skip_if_not_installed("caret")
  set.seed(5)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- factor(ifelse(X$x1 + X$x2 + rnorm(n, sd = 0.5) > 0, "A", "B"))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- factor(ifelse(Y$x1 + Y$x2 + rnorm(n, sd = 0.5) > 0, "A", "B"))

  res <- tstr(X, Y, target_var = "y", seed = 5)

  expect_equal(res$metric, "AUC")
  expect_true(res$tstr_performance >= 0 && res$tstr_performance <= 1)
  expect_true(res$trtr_performance >= 0 && res$trtr_performance <= 1)
})

test_that("multiple target variables are handled", {
  skip_if_not_installed("caret")
  set.seed(6)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y1 <- 2 * X$x1 + rnorm(n, sd = 0.5)
  X$y2 <- -X$x2 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y1 <- 2 * Y$x1 + rnorm(n, sd = 0.5)
  Y$y2 <- -Y$x2 + rnorm(n, sd = 0.5)

  res <- tstr(X, Y, target_var = c("y1", "y2"), seed = 6)

  expect_s3_class(res, "tstr")
  expect_equal(length(res$per_target), 2)
  expect_true("y1" %in% names(res$per_target))
  expect_true("y2" %in% names(res$per_target))
  # Aggregate ratio is average of per-target ratios
  avg_ratio <- mean(c(res$per_target$y1$tstr_ratio, res$per_target$y2$tstr_ratio))
  expect_equal(res$tstr_ratio, avg_ratio)
})

test_that("test_fraction parameter controls split", {
  skip_if_not_installed("caret")
  set.seed(7)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- X$x1 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- Y$x1 + rnorm(n, sd = 0.5)

  res <- tstr(X, Y, target_var = "y", test_fraction = 0.2, seed = 7)

  # With n=200 and test_fraction=0.2, expect ~40 test, ~160 train
  expect_equal(res$n_test, 40)
  expect_equal(res$n_train_orig, 160)
  expect_equal(res$test_fraction, 0.2)
})

test_that("pre-split test_data overrides test_fraction", {
  skip_if_not_installed("caret")
  set.seed(8)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- X$x1 + rnorm(n, sd = 0.5)

  # Pre-split test data
  test_df <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  test_df$y <- test_df$x1 + rnorm(50, sd = 0.5)

  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- Y$x1 + rnorm(n, sd = 0.5)

  res <- tstr(X, Y, target_var = "y", test_data = test_df, seed = 8)

  # Full X used for training, pre-split test_data used for test
  expect_equal(res$n_train_orig, 200)
  expect_equal(res$n_test, 50)
})

test_that("synth_pair dispatch works", {
  skip_if_not_installed("caret")
  set.seed(9)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- X$x1 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- Y$x1 + rnorm(n, sd = 0.5)

  sp <- list(original = X, synthetic = Y)
  class(sp) <- "synth_pair"

  res <- tstr(sp, target_var = "y", seed = 9)
  expect_s3_class(res, "tstr")
  expect_equal(res$n_synth, 200)
})

test_that("print method works without error", {
  skip_if_not_installed("caret")
  set.seed(10)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- X$x1 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- Y$x1 + rnorm(n, sd = 0.5)
  res <- tstr(X, Y, target_var = "y", seed = 10)

  expect_output(print(res), "Train on Synthetic")
  expect_output(print(res), "TRTR")
  expect_output(print(res), "TSTR")
  expect_output(print(res), "Utility")
})

test_that("print returns invisible x", {
  skip_if_not_installed("caret")
  set.seed(11)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- X$x1 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- Y$x1 + rnorm(n, sd = 0.5)
  res <- tstr(X, Y, target_var = "y", seed = 11)
  out <- capture.output(ret <- print(res))
  expect_identical(ret, res)
})

test_that("summary method returns summary.tstr object", {
  skip_if_not_installed("caret")
  set.seed(12)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- X$x1 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- Y$x1 + rnorm(n, sd = 0.5)
  res <- tstr(X, Y, target_var = "y", seed = 12)

  s <- summary(res)
  expect_s3_class(s, "summary.tstr")
  expect_true("tstr_ratio" %in% names(s))
  expect_true("utility_score" %in% names(s))
  expect_true("performance_gap" %in% names(s))
})

test_that("print.summary.tstr works", {
  skip_if_not_installed("caret")
  set.seed(13)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- X$x1 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- Y$x1 + rnorm(n, sd = 0.5)
  res <- tstr(X, Y, target_var = "y", seed = 13)

  s <- summary(res)
  expect_output(print(s), "Summary.*TSTR")
  expect_output(print(s), "Performance Comparison")
})

test_that("plot method works without error", {
  skip_if_not_installed("caret")
  set.seed(14)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- X$x1 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- Y$x1 + rnorm(n, sd = 0.5)
  res <- tstr(X, Y, target_var = "y", seed = 14)

  expect_silent(plot(res, which = 1))
})

test_that("plot method works with multiple targets", {
  skip_if_not_installed("caret")
  set.seed(15)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y1 <- X$x1 + rnorm(n, sd = 0.5)
  X$y2 <- X$x2 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y1 <- Y$x1 + rnorm(n, sd = 0.5)
  Y$y2 <- Y$x2 + rnorm(n, sd = 0.5)
  res <- tstr(X, Y, target_var = c("y1", "y2"), seed = 15)

  expect_silent(plot(res, which = 1))
})

test_that("error without target_var", {
  skip_if_not_installed("caret")
  X <- data.frame(x = rnorm(100), y = rnorm(100))
  Y <- data.frame(x = rnorm(100), y = rnorm(100))

  expect_error(tstr(X, Y), "target_var.*required")
})

test_that("error with missing target variable in data", {
  skip_if_not_installed("caret")
  X <- data.frame(x = rnorm(100), y = rnorm(100))
  Y <- data.frame(x = rnorm(100), z = rnorm(100))

  expect_error(tstr(X, Y, target_var = "y"), "missing in Y")
})

test_that("NA handling removes rows", {
  skip_if_not_installed("caret")
  set.seed(16)
  n <- 200
  X <- data.frame(x1 = c(rnorm(n - 5), rep(NA, 5)), x2 = rnorm(n))
  X$y <- ifelse(is.na(X$x1), NA, X$x1 + rnorm(n, sd = 0.5))
  Y <- data.frame(x1 = c(rnorm(n - 3), rep(NA, 3)), x2 = rnorm(n))
  Y$y <- ifelse(is.na(Y$x1), NA, Y$x1 + rnorm(n, sd = 0.5))

  res <- tstr(X, Y, target_var = "y", na.rm = TRUE, seed = 16)

  expect_s3_class(res, "tstr")
  # After NA removal, we should have fewer total rows
  expect_true(res$n_train_orig + res$n_test <= 195)
  expect_true(res$n_synth <= 197)
})

test_that("model parameter selects correct model", {
  skip_if_not_installed("caret")
  set.seed(17)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- X$x1 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- Y$x1 + rnorm(n, sd = 0.5)

  res <- tstr(X, Y, target_var = "y", model = "glm", seed = 17)
  expect_equal(res$model, "glm")
})

test_that("rpart model works", {
  skip_if_not_installed("caret")
  skip_if_not_installed("rpart")
  set.seed(18)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- X$x1 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- Y$x1 + rnorm(n, sd = 0.5)

  res <- tstr(X, Y, target_var = "y", model = "rpart", seed = 18)
  expect_s3_class(res, "tstr")
  expect_equal(res$model, "rpart")
})

test_that("seed parameter gives reproducible results", {
  skip_if_not_installed("caret")
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y <- X$x1 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y <- Y$x1 + rnorm(n, sd = 0.5)

  res1 <- tstr(X, Y, target_var = "y", seed = 42)
  res2 <- tstr(X, Y, target_var = "y", seed = 42)
  expect_equal(res1$tstr_ratio, res2$tstr_ratio)
  expect_equal(res1$utility_score, res2$utility_score)
})

test_that("R2 metric is used for regression", {
  skip_if_not_installed("caret")
  set.seed(20)
  n <- 200
  X <- data.frame(x1 = rnorm(n))
  X$y <- X$x1 + rnorm(n, sd = 0.3)
  Y <- data.frame(x1 = rnorm(n))
  Y$y <- Y$x1 + rnorm(n, sd = 0.3)

  res <- tstr(X, Y, target_var = "y", seed = 20)
  expect_equal(res$metric, "R2")
})

test_that("per_target has correct structure for single target", {
  skip_if_not_installed("caret")
  set.seed(21)
  n <- 200
  X <- data.frame(x1 = rnorm(n))
  X$y <- X$x1 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n))
  Y$y <- Y$x1 + rnorm(n, sd = 0.5)

  res <- tstr(X, Y, target_var = "y", seed = 21)
  expect_equal(length(res$per_target), 1)
  expect_true("tstr_performance" %in% names(res$per_target$y))
  expect_true("trtr_performance" %in% names(res$per_target$y))
  expect_true("tstr_ratio" %in% names(res$per_target$y))
  expect_true("metric" %in% names(res$per_target$y))
})

test_that("print shows per-target info for multiple targets", {
  skip_if_not_installed("caret")
  set.seed(22)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  X$y1 <- X$x1 + rnorm(n, sd = 0.5)
  X$y2 <- X$x2 + rnorm(n, sd = 0.5)
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y$y1 <- Y$x1 + rnorm(n, sd = 0.5)
  Y$y2 <- Y$x2 + rnorm(n, sd = 0.5)
  res <- tstr(X, Y, target_var = c("y1", "y2"), seed = 22)

  expect_output(print(res), "Per-Target")
  expect_output(print(res), "y1")
  expect_output(print(res), "y2")
})

test_that("classification with character target works", {
  skip_if_not_installed("caret")
  set.seed(23)
  n <- 200
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), stringsAsFactors = FALSE)
  X$y <- ifelse(X$x1 + rnorm(n, sd = 2) > 0, "pos", "neg")
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), stringsAsFactors = FALSE)
  Y$y <- ifelse(Y$x1 + rnorm(n, sd = 2) > 0, "pos", "neg")

  res <- suppressWarnings(tstr(X, Y, target_var = "y", seed = 23))
  expect_equal(res$metric, "AUC")
  expect_s3_class(res, "tstr")
})

# Tests for pMSE (Propensity Score Mean Squared Error)

library(testthat)

test_that("pMSE returns correct S3 class structure", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    income = rnorm(n, 50000, 10000)
  )
  Y <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    income = rnorm(n, 50000, 10000)
  )

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  expect_s3_class(result, "pMSE")
  expect_true("pMSE" %in% names(result))
  expect_true("pMSE_null" %in% names(result))
  expect_true("pMSE_ratio" %in% names(result))
  expect_true("S_pMSE" %in% names(result))
  expect_true("SPECKS" %in% names(result))
  expect_true("ks_pvalue" %in% names(result))
  expect_true("PO50" %in% names(result))
  expect_true("propensity_original" %in% names(result))
  expect_true("propensity_synthetic" %in% names(result))
  expect_true("c" %in% names(result))
  expect_true("k" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("utility_interpretation" %in% names(result))
})

test_that("pMSE values are non-negative", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  expect_true(result$pMSE >= 0)
  expect_true(result$pMSE_null >= 0)
  expect_true(result$SPECKS >= 0 && result$SPECKS <= 1)
})

test_that("pMSE propensity scores are in [0, 1]", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  expect_true(all(result$propensity_original >= 0 & result$propensity_original <= 1))
  expect_true(all(result$propensity_synthetic >= 0 & result$propensity_synthetic <= 1))
})

test_that("pMSE propensity scores have correct lengths", {
  set.seed(123)
  n_x <- 80
  n_y <- 120
  X <- data.frame(x1 = rnorm(n_x))
  Y <- data.frame(x1 = rnorm(n_y))

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  expect_equal(length(result$propensity_original), n_x)
  expect_equal(length(result$propensity_synthetic), n_y)
  expect_equal(result$n_original, n_x)
  expect_equal(result$n_synthetic, n_y)
  expect_equal(result$N, n_x + n_y)
})

test_that("pMSE c proportion is correctly computed", {
  set.seed(123)
  n_x <- 80
  n_y <- 120
  X <- data.frame(x1 = rnorm(n_x))
  Y <- data.frame(x1 = rnorm(n_y))

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  expect_equal(result$c, n_y / (n_x + n_y))
})

test_that("pMSE ratio is close to 1 for similar data (logit)", {
  set.seed(123)
  n <- 200
  X <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE))
  )
  Y <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE))
  )

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  # For data from same distribution, ratio should be moderate (roughly around 1)
  expect_true(result$pMSE_ratio < 10)
})

test_that("pMSE PO50 is near 50% for similar data", {
  set.seed(123)
  n <- 200
  X <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE))
  )
  Y <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE))
  )

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  # PO50 should be roughly around 50% for indistinguishable data
  expect_true(result$PO50 >= 30 && result$PO50 <= 70)
})

test_that("pMSE detects different distributions", {
  set.seed(123)
  n <- 200
  X <- data.frame(
    age = sample(18:30, n, replace = TRUE),
    income = rnorm(n, 30000, 5000)
  )
  Y <- data.frame(
    age = sample(60:80, n, replace = TRUE),
    income = rnorm(n, 80000, 5000)
  )

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  # pMSE should be noticeably larger for very different data
  expect_true(result$pMSE > 0.01)
  expect_true(result$SPECKS > 0.1)
})

test_that("pMSE works with logit method and maxorder = 1", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = factor(sample(c("A", "B"), n, replace = TRUE))
  )
  Y <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = factor(sample(c("A", "B"), n, replace = TRUE))
  )

  result <- pMSE(X, Y, method = "logit", maxorder = 1, seed = 42)

  expect_s3_class(result, "pMSE")
  expect_equal(result$method, "logit")
  expect_equal(result$maxorder, 1)
  # With interactions, k should be larger than main effects only
  result_main <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)
  expect_true(result$k >= result_main$k)
})

test_that("pMSE works with cart method", {
  skip_if_not_installed("rpart")
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- pMSE(X, Y, method = "cart", nperms = 5, seed = 42)

  expect_s3_class(result, "pMSE")
  expect_equal(result$method, "cart")
  expect_true(is.na(result$k))  # k is NA for CART
  # With nperms > 0, null should be estimated
  expect_false(is.na(result$pMSE_null))
})

test_that("pMSE cart with nperms = 0 gives NA null", {
  skip_if_not_installed("rpart")
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- pMSE(X, Y, method = "cart", nperms = 0, seed = 42)

  expect_true(is.na(result$pMSE_null))
  expect_true(is.na(result$pMSE_ratio))
  expect_true(is.na(result$S_pMSE))
})

test_that("pMSE works with a single variable", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n))

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  expect_s3_class(result, "pMSE")
  expect_equal(length(result$vars), 1)
})

test_that("pMSE works with explicit vars parameter", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- pMSE(X, Y, vars = c("x1", "x2"), method = "logit", maxorder = 0, seed = 42)

  expect_equal(length(result$vars), 2)
  expect_true("x1" %in% result$vars)
  expect_true("x2" %in% result$vars)
})

test_that("pMSE input validation: X must be data frame", {
  Y <- data.frame(x1 = 1:10)
  expect_error(pMSE(1:10, Y), "X must be a data frame")
})

test_that("pMSE input validation: Y must be data frame", {
  X <- data.frame(x1 = 1:10)
  expect_error(pMSE(X, 1:10), "Y must be a data frame")
})

test_that("pMSE errors when no common variables", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(b = 1:10)

  expect_error(pMSE(X, Y), "No common variables")
})

test_that("pMSE errors when specified variables missing in X", {
  X <- data.frame(x1 = 1:10)
  Y <- data.frame(x1 = 1:10, x2 = rnorm(10))

  expect_error(pMSE(X, Y, vars = c("x1", "x2")), "Variables missing in X")
})

test_that("pMSE errors when specified variables missing in Y", {
  X <- data.frame(x1 = 1:10, x2 = rnorm(10))
  Y <- data.frame(x1 = 1:10)

  expect_error(pMSE(X, Y, vars = c("x1", "x2")), "Variables missing in Y")
})

test_that("pMSE handles NA values via removal", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = c(rnorm(n - 3), NA, NA, NA), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = c(rnorm(n - 1), NA))

  expect_message(
    result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42, na.rm = TRUE),
    "Removing"
  )

  expect_s3_class(result, "pMSE")
  # Number of original + synthetic should be less than 2n
  expect_true(result$N < 2 * n)
})

test_that("pMSE print method works", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  expect_output(print(result), "Propensity Score Mean Squared Error")
  expect_output(print(result), "pMSE")
  expect_output(print(result), "SPECKS")
  expect_output(print(result), "PO50")
  expect_output(print(result), "Interpretation")
})

test_that("pMSE summary method returns correct class", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)
  s <- summary(result)

  expect_s3_class(s, "summary.pMSE")
  expect_true(!is.null(s$pMSE))
  expect_true(!is.null(s$pMSE_null))
  expect_true(!is.null(s$pMSE_ratio))
  expect_true(!is.null(s$SPECKS))
  expect_true(!is.null(s$PO50))
  expect_true(!is.null(s$prop_original_summary))
  expect_true(!is.null(s$prop_synthetic_summary))
  expect_true(!is.null(s$method))
  expect_true(!is.null(s$vars))
  expect_output(print(s), "Summary: Propensity Score Mean Squared Error")
})

test_that("pMSE plot method works for all which values", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 3))
  expect_no_error(plot(result, which = 4))
  expect_no_error(plot(result, which = 1:4))
})

test_that("pMSE works via synth_pair dispatch", {
  set.seed(123)
  n <- 100
  original <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    income = rnorm(n, 50000, 10000)
  )
  synthetic <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    income = rnorm(n, 50000, 10000)
  )

  pair <- synth_pair(original, synthetic)
  result <- pMSE(pair, method = "logit", maxorder = 0, seed = 42)

  expect_s3_class(result, "pMSE")
  expect_equal(result$n_original, n)
  expect_equal(result$n_synthetic, n)
})

test_that("pMSE logit analytical null distribution is correctly computed", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  # Verify null formula: E[pMSE] = (k-1) * (1-c)^2 * c / N
  expected_null <- (result$k - 1) * (1 - result$c)^2 * result$c / result$N
  expect_equal(result$pMSE_null, expected_null, tolerance = 1e-10)
})

test_that("pMSE seed parameter ensures reproducibility", {
  set.seed(1)
  n <- 100
  X <- data.frame(
    x1 = rnorm(n),
    x2 = factor(sample(c("A", "B"), n, replace = TRUE))
  )
  Y <- data.frame(
    x1 = rnorm(n),
    x2 = factor(sample(c("A", "B"), n, replace = TRUE))
  )

  result1 <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)
  result2 <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  expect_equal(result1$pMSE, result2$pMSE)
  expect_equal(result1$SPECKS, result2$SPECKS)
})

test_that("pMSE interpretation string is non-empty", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n))
  Y <- data.frame(x1 = rnorm(n))

  result <- pMSE(X, Y, method = "logit", maxorder = 0, seed = 42)

  expect_true(nchar(result$utility_interpretation) > 0)
  # Should start with a quality label
  expect_true(grepl("^(GOOD|MODERATE|FAIR|POOR)", result$utility_interpretation))
})

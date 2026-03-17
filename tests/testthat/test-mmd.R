test_that("mmd returns correct S3 class with expected fields", {
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  res <- mmd(X, Y, seed = 1)

  expect_s3_class(res, "mmd")
  expect_true("mmd2" %in% names(res))
  expect_true("sigma" %in% names(res))
  expect_true("kernel" %in% names(res))
  expect_true("method" %in% names(res))
  expect_true("n_X" %in% names(res))
  expect_true("n_Y" %in% names(res))
  expect_true("n_vars" %in% names(res))
  expect_true("vars" %in% names(res))
  expect_true("standardized" %in% names(res))
  expect_true("utility_score" %in% names(res))
  expect_true("perm_pvalue" %in% names(res))
  expect_true("perm_null" %in% names(res))
})

test_that("mmd2 is non-negative (or near zero)", {
  set.seed(2)
  X <- data.frame(x = rnorm(200), y = rnorm(200))
  Y <- data.frame(x = rnorm(200), y = rnorm(200))
  res <- mmd(X, Y, seed = 2)
  # Unbiased estimator can be slightly negative, but should be close to zero

  expect_true(res$mmd2 > -0.1)
})

test_that("mmd2 is near zero for identical distributions (large sample)", {
  set.seed(3)
  n <- 1000
  X <- data.frame(x = rnorm(n), y = rnorm(n))
  Y <- data.frame(x = rnorm(n), y = rnorm(n))
  res <- mmd(X, Y, seed = 3)
  # Should be close to zero
  expect_true(abs(res$mmd2) < 0.05)
  expect_true(res$utility_score > 0.8)
})

test_that("mmd detects different distributions", {
  set.seed(4)
  X <- data.frame(x = rnorm(300, 0, 1), y = rnorm(300, 0, 1))
  Y <- data.frame(x = rnorm(300, 5, 1), y = rnorm(300, 5, 1))
  res <- mmd(X, Y, seed = 4)
  expect_true(res$mmd2 > 0.01)
  # Utility score should be noticeably below 1 for shifted distributions
  expect_true(res$utility_score < 0.8)
})

test_that("utility_score is in [0,1]", {
  set.seed(5)
  X <- data.frame(x = rnorm(100), y = rnorm(100))
  # Similar
  Y1 <- data.frame(x = rnorm(100), y = rnorm(100))
  res1 <- mmd(X, Y1, seed = 5)
  expect_true(res1$utility_score >= 0 && res1$utility_score <= 1)

  # Very different
  Y2 <- data.frame(x = rnorm(100, 10, 1), y = rnorm(100, 10, 1))
  res2 <- mmd(X, Y2, seed = 5)
  expect_true(res2$utility_score >= 0 && res2$utility_score <= 1)
})

test_that("mmd auto-detects numeric variables and skips categoricals", {
  set.seed(6)
  X <- data.frame(
    num1 = rnorm(100),
    num2 = rnorm(100),
    cat1 = sample(letters[1:3], 100, replace = TRUE),
    stringsAsFactors = FALSE
  )
  Y <- data.frame(
    num1 = rnorm(100),
    num2 = rnorm(100),
    cat1 = sample(letters[1:3], 100, replace = TRUE),
    stringsAsFactors = FALSE
  )
  res <- mmd(X, Y, seed = 6)
  expect_equal(res$vars, c("num1", "num2"))
  expect_equal(res$n_vars, 2)
})

test_that("vars parameter selects specific variables", {
  set.seed(7)
  X <- data.frame(a = rnorm(100), b = rnorm(100), c = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100), c = rnorm(100))
  res <- mmd(X, Y, vars = c("a", "c"), seed = 7)
  expect_equal(res$vars, c("a", "c"))
  expect_equal(res$n_vars, 2)
})

test_that("synth_pair dispatch works", {
  set.seed(8)
  X <- data.frame(x = rnorm(100), y = rnorm(100))
  Y <- data.frame(x = rnorm(100), y = rnorm(100))
  sp <- list(
    original = X,
    synthetic = Y,
    num_vars = c("x", "y")
  )
  class(sp) <- "synth_pair"
  res <- mmd(sp, seed = 8)
  expect_s3_class(res, "mmd")
  expect_equal(res$n_vars, 2)
})

test_that("permutation test returns valid p-value", {
  set.seed(9)
  X <- data.frame(x = rnorm(100), y = rnorm(100))
  Y <- data.frame(x = rnorm(100, 3, 1), y = rnorm(100, 3, 1))
  res <- mmd(X, Y, n_perm = 100, seed = 9)
  expect_true(!is.null(res$perm_pvalue))
  expect_true(res$perm_pvalue >= 0 && res$perm_pvalue <= 1)
  expect_length(res$perm_null, 100)
  # With different distributions, p should be small
  expect_true(res$perm_pvalue < 0.1)
})

test_that("permutation p-value is large for identical distributions", {
  set.seed(10)
  X <- data.frame(x = rnorm(200), y = rnorm(200))
  Y <- data.frame(x = rnorm(200), y = rnorm(200))
  res <- mmd(X, Y, n_perm = 200, seed = 10)
  # Should not reject
  expect_true(res$perm_pvalue > 0.01)
})

test_that("RFF method works and approximates exact", {
  set.seed(11)
  X <- data.frame(x = rnorm(200), y = rnorm(200))
  Y <- data.frame(x = rnorm(200, 2, 1), y = rnorm(200, 2, 1))

  res_exact <- mmd(X, Y, method = "exact", seed = 11)
  res_rff <- mmd(X, Y, method = "rff", n_features = 2000, seed = 11)

  expect_s3_class(res_rff, "mmd")
  expect_equal(res_rff$method, "rff")
  expect_equal(res_rff$n_features, 2000L)

  # Both should detect the difference
  expect_true(res_exact$mmd2 > 0.01)
  expect_true(res_rff$mmd2 > 0.01)

  # RFF should be a reasonable approximation (within factor of 5)
  expect_true(res_rff$mmd2 / res_exact$mmd2 > 0.2)
  expect_true(res_rff$mmd2 / res_exact$mmd2 < 5)
})

test_that("exact method emits message for n > 5000", {
  set.seed(12)
  X <- data.frame(x = rnorm(5001))
  Y <- data.frame(x = rnorm(100))
  expect_message(mmd(X, Y, seed = 12), "rff")
})

test_that("gaussian kernel is default", {
  set.seed(13)
  X <- data.frame(x = rnorm(100))
  Y <- data.frame(x = rnorm(100))
  res <- mmd(X, Y, seed = 13)
  expect_equal(res$kernel, "gaussian")
})

test_that("rational_quadratic kernel works", {
  set.seed(14)
  X <- data.frame(x = rnorm(100), y = rnorm(100))
  Y <- data.frame(x = rnorm(100, 2), y = rnorm(100, 2))
  res <- mmd(X, Y, kernel = "rational_quadratic", seed = 14)
  expect_s3_class(res, "mmd")
  expect_equal(res$kernel, "rational_quadratic")
  expect_true(res$mmd2 > 0)
})

test_that("standardize parameter controls standardization", {
  set.seed(15)
  X <- data.frame(x = rnorm(100, 0, 1), y = rnorm(100, 0, 1000))
  Y <- data.frame(x = rnorm(100, 0, 1), y = rnorm(100, 0, 1000))
  res_std <- mmd(X, Y, standardize = TRUE, seed = 15)
  res_raw <- mmd(X, Y, standardize = FALSE, seed = 15)

  expect_true(res_std$standardized)
  expect_false(res_raw$standardized)

  # The sigma values should differ due to standardization

  expect_true(res_std$sigma != res_raw$sigma)
})

test_that("print method works without error", {
  set.seed(16)
  X <- data.frame(x = rnorm(100))
  Y <- data.frame(x = rnorm(100))
  res <- mmd(X, Y, seed = 16)
  expect_output(print(res), "Maximum Mean Discrepancy")
  expect_output(print(res), "MMD\\^2")
  expect_output(print(res), "Utility")
})

test_that("print method shows permutation info when available", {
  set.seed(17)
  X <- data.frame(x = rnorm(100))
  Y <- data.frame(x = rnorm(100, 3))
  res <- mmd(X, Y, n_perm = 50, seed = 17)
  expect_output(print(res), "Perm")
})

test_that("summary method returns summary.mmd object", {
  set.seed(18)
  X <- data.frame(x = rnorm(100), y = rnorm(100))
  Y <- data.frame(x = rnorm(100), y = rnorm(100))
  res <- mmd(X, Y, seed = 18)
  s <- summary(res)
  expect_s3_class(s, "summary.mmd")
  expect_true("mmd2" %in% names(s))
  expect_true("utility_score" %in% names(s))
})

test_that("print.summary.mmd works", {
  set.seed(19)
  X <- data.frame(x = rnorm(100))
  Y <- data.frame(x = rnorm(100))
  res <- mmd(X, Y, n_perm = 50, seed = 19)
  s <- summary(res)
  expect_output(print(s), "Summary.*MMD")
  expect_output(print(s), "Permutation Test")
})

test_that("plot method works with permutation results", {
  set.seed(20)
  X <- data.frame(x = rnorm(100))
  Y <- data.frame(x = rnorm(100, 3))
  res <- mmd(X, Y, n_perm = 50, seed = 20)
  # Should not error
  expect_silent(plot(res, which = 1))
})

test_that("plot method gives message when no permutation results", {
  set.seed(21)
  X <- data.frame(x = rnorm(100))
  Y <- data.frame(x = rnorm(100))
  res <- mmd(X, Y, seed = 21)
  expect_message(plot(res, which = 1), "No permutation results")
})

test_that("NA handling with na.rm = TRUE", {
  set.seed(22)
  X <- data.frame(x = c(rnorm(98), NA, NA), y = rnorm(100))
  Y <- data.frame(x = rnorm(100), y = c(rnorm(99), NA))
  res <- mmd(X, Y, na.rm = TRUE, seed = 22)
  expect_s3_class(res, "mmd")
  expect_equal(res$n_X, 98)
  expect_equal(res$n_Y, 99)
})

test_that("error when no numeric variables", {
  X <- data.frame(a = letters[1:10], b = LETTERS[1:10], stringsAsFactors = FALSE)
  Y <- data.frame(a = letters[1:10], b = LETTERS[1:10], stringsAsFactors = FALSE)
  expect_error(mmd(X, Y), "No numeric variables")
})

test_that("error when specified vars are non-numeric", {
  X <- data.frame(a = letters[1:10], b = 1:10, stringsAsFactors = FALSE)
  Y <- data.frame(a = letters[1:10], b = 1:10, stringsAsFactors = FALSE)
  expect_error(mmd(X, Y, vars = c("a")), "Non-numeric variables")
})

test_that("error when specified vars missing from data", {
  X <- data.frame(a = rnorm(10))
  Y <- data.frame(a = rnorm(10))
  expect_error(mmd(X, Y, vars = c("a", "nonexistent")), "missing")
})

test_that("sigma is positive", {
  set.seed(23)
  X <- data.frame(x = rnorm(100))
  Y <- data.frame(x = rnorm(100))
  res <- mmd(X, Y, seed = 23)
  expect_true(res$sigma > 0)
})

test_that("seed parameter gives reproducible results", {
  X <- data.frame(x = rnorm(100), y = rnorm(100))
  Y <- data.frame(x = rnorm(100), y = rnorm(100))
  res1 <- mmd(X, Y, seed = 42)
  res2 <- mmd(X, Y, seed = 42)
  expect_equal(res1$mmd2, res2$mmd2)
  expect_equal(res1$sigma, res2$sigma)
})

test_that("RFF with permutation test works", {
  set.seed(24)
  X <- data.frame(x = rnorm(100))
  Y <- data.frame(x = rnorm(100, 3))
  res <- mmd(X, Y, method = "rff", n_perm = 50, seed = 24)
  expect_true(!is.null(res$perm_pvalue))
  expect_true(res$perm_pvalue < 0.1)
})

test_that("mmd works with single variable", {
  set.seed(25)
  X <- data.frame(x = rnorm(200))
  Y <- data.frame(x = rnorm(200))
  res <- mmd(X, Y, seed = 25)
  expect_s3_class(res, "mmd")
  expect_equal(res$n_vars, 1)
})

test_that("mmd works with many variables", {
  set.seed(26)
  X <- as.data.frame(matrix(rnorm(100 * 10), ncol = 10))
  Y <- as.data.frame(matrix(rnorm(100 * 10), ncol = 10))
  res <- mmd(X, Y, seed = 26)
  expect_s3_class(res, "mmd")
  expect_equal(res$n_vars, 10)
})

test_that("mmd works with unequal sample sizes", {
  set.seed(27)
  X <- data.frame(x = rnorm(50), y = rnorm(50))
  Y <- data.frame(x = rnorm(200), y = rnorm(200))
  res <- mmd(X, Y, seed = 27)
  expect_equal(res$n_X, 50)
  expect_equal(res$n_Y, 200)
})

test_that("rational_quadratic detects shifted distributions", {
  set.seed(28)
  X <- data.frame(x = rnorm(200), y = rnorm(200))
  Y <- data.frame(x = rnorm(200, 5), y = rnorm(200, 5))

  res_gauss <- mmd(X, Y, kernel = "gaussian", seed = 28)
  res_rq <- mmd(X, Y, kernel = "rational_quadratic", seed = 28)

  expect_true(res_gauss$mmd2 > 0.01)
  expect_true(res_rq$mmd2 > 0.01)
})

test_that("print returns invisible x", {
  set.seed(29)
  X <- data.frame(x = rnorm(50))
  Y <- data.frame(x = rnorm(50))
  res <- mmd(X, Y, seed = 29)
  out <- capture.output(ret <- print(res))
  expect_identical(ret, res)
})

test_that("mmd method field is correct", {
  set.seed(30)
  X <- data.frame(x = rnorm(100))
  Y <- data.frame(x = rnorm(100))
  res_exact <- mmd(X, Y, method = "exact", seed = 30)
  res_rff <- mmd(X, Y, method = "rff", seed = 30)
  expect_equal(res_exact$method, "exact")
  expect_equal(res_rff$method, "rff")
  expect_null(res_exact$n_features)
  expect_equal(res_rff$n_features, 500L)
})

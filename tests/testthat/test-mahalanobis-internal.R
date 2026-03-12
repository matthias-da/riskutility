# tests/testthat/test-mahalanobis-internal.R

# --- .mahal_prepare() tests ---

test_that(".mahal_prepare returns correct structure with robust = TRUE", {
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50))
  type <- c(a = "numeric", b = "numeric", c = "numeric")
  res <- .mahal_prepare(X, names(X), type, robust = TRUE)

  expect_type(res, "list")
  expect_true(is.matrix(res$cov_inv))
  expect_equal(dim(res$cov_inv), c(3, 3))
  expect_true(is.numeric(res$center))
  expect_equal(length(res$center), 3)
  expect_equal(res$numeric_keys, c("a", "b", "c"))
  expect_equal(res$nominal_keys, character(0))
  expect_true(res$robust)
  expect_true(is.numeric(res$chi_sq_threshold))
})

test_that(".mahal_prepare with robust = FALSE uses classical covariance", {
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50))
  type <- c(a = "numeric", b = "numeric")
  res <- .mahal_prepare(X, names(X), type, robust = FALSE)

  expect_false(res$robust)
  expect_true(is.matrix(res$cov_inv))
  expect_equal(dim(res$cov_inv), c(2, 2))
})

test_that(".mahal_prepare separates numeric and nominal keys", {
  set.seed(1)
  X <- data.frame(
    num1 = rnorm(40), num2 = rnorm(40),
    cat1 = sample(letters[1:3], 40, TRUE),
    stringsAsFactors = FALSE
  )
  type <- c(num1 = "numeric", num2 = "numeric", cat1 = "nominal")
  res <- .mahal_prepare(X, names(X), type, robust = TRUE)

  expect_equal(res$numeric_keys, c("num1", "num2"))
  expect_equal(res$nominal_keys, "cat1")
  expect_equal(dim(res$cov_inv), c(2, 2))
  expect_true(res$alpha > 0 && res$alpha < 1)
})

test_that(".mahal_prepare handles ordinal as numeric", {
  set.seed(1)
  X <- data.frame(
    a = ordered(sample(c("low", "mid", "high"), 40, TRUE),
                levels = c("low", "mid", "high")),
    b = rnorm(40)
  )
  type <- c(a = "ordinal", b = "numeric")
  res <- .mahal_prepare(X, names(X), type, robust = TRUE)

  expect_equal(res$numeric_keys, c("a", "b"))
  expect_equal(dim(res$cov_inv), c(2, 2))
})

test_that(".mahal_prepare errors on purely nominal data", {
  X <- data.frame(
    a = sample(letters[1:3], 30, TRUE),
    b = sample(c("x", "y"), 30, TRUE),
    stringsAsFactors = FALSE
  )
  type <- c(a = "nominal", b = "nominal")
  expect_error(
    .mahal_prepare(X, names(X), type, robust = TRUE),
    "at least one numeric"
  )
})

test_that(".mahal_prepare falls back when p >= n for MCD", {
  set.seed(1)
  # 5 observations, 6 variables -> MCD can't work
  X <- data.frame(matrix(rnorm(30), nrow = 5, ncol = 6))
  type <- setNames(rep("numeric", 6), names(X))
  # Two warnings: fallback + singular covariance
  suppressWarnings(
    expect_warning(
      res <- .mahal_prepare(X, names(X), type, robust = TRUE),
      "Falling back"
    )
  )
  expect_true(is.matrix(res$cov_inv))
  expect_false(res$robust)
})

test_that(".mahal_prepare handles near-singular covariance", {
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50))
  X$c <- X$a + X$b  # perfectly collinear
  type <- c(a = "numeric", b = "numeric", c = "numeric")
  # MCD and our code both warn about singularity
  suppressWarnings(
    res <- .mahal_prepare(X, names(X), type, robust = TRUE)
  )
  expect_true(is.matrix(res$cov_inv))
  expect_true(all(is.finite(res$cov_inv)))
})

test_that(".mahal_prepare computes chi-squared threshold", {
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50))
  type <- c(a = "numeric", b = "numeric", c = "numeric")
  res <- .mahal_prepare(X, names(X), type, robust = TRUE)

  expected <- sqrt(qchisq(0.975, df = 3))
  expect_equal(res$chi_sq_threshold, expected)
})

# --- .mahal_dist() tests ---

test_that(".mahal_dist computes correct distances for pure numeric", {
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  type <- c(a = "numeric", b = "numeric")
  prep <- .mahal_prepare(X, names(X), type, robust = FALSE)

  x_row <- X[1, , drop = FALSE]
  candidates <- X[1:5, , drop = FALSE]
  di <- .mahal_dist(x_row, candidates, prep, type)

  expect_equal(length(di), 5)
  expect_true(all(di >= 0))
  expect_true(di[1] < 1e-10)  # distance to self should be ~0
})

test_that(".mahal_dist with mixed data combines Mahalanobis and nominal", {
  set.seed(1)
  X <- data.frame(
    num = rnorm(30),
    cat = sample(c("a", "b"), 30, TRUE),
    stringsAsFactors = FALSE
  )
  type <- c(num = "numeric", cat = "nominal")
  prep <- .mahal_prepare(X, names(X), type, robust = FALSE)

  x_row <- X[1, , drop = FALSE]
  candidates <- X[1:5, , drop = FALSE]
  di <- .mahal_dist(x_row, candidates, prep, type)

  expect_equal(length(di), 5)
  expect_true(all(di >= 0))
  expect_true(di[1] < 1e-10)  # self-distance = 0
})

test_that(".mahal_dist accounts for correlation structure", {
  set.seed(42)
  n <- 100
  a <- rnorm(n)
  b <- a + rnorm(n, 0, 0.1)  # b ~= a (highly correlated)
  X <- data.frame(a = a, b = b)
  type <- c(a = "numeric", b = "numeric")
  prep <- .mahal_prepare(X, names(X), type, robust = FALSE)

  x_base <- data.frame(a = 0, b = 0)
  cand_along <- data.frame(a = 1, b = 1)     # along correlation
  cand_against <- data.frame(a = 1, b = -1)  # against correlation

  d_along <- .mahal_dist(x_base, cand_along, prep, type)
  d_against <- .mahal_dist(x_base, cand_against, prep, type)

  # Against-correlation should be larger (more unusual)
  expect_true(d_against > d_along)
})

test_that(".mahal_dist normalizes to reasonable range", {
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50))
  type <- c(a = "numeric", b = "numeric")
  prep <- .mahal_prepare(X, names(X), type, robust = FALSE)

  all_dists <- numeric(0)
  for (i in 1:10) {
    di <- .mahal_dist(X[i, , drop = FALSE], X, prep, type)
    all_dists <- c(all_dists, di)
  }
  expect_true(all(all_dists >= 0))
  # Most in-distribution distances should be <= 1 after chi-sq normalization
  expect_true(median(all_dists) < 2)
})

test_that(".mahal_dist handles single numeric key", {
  set.seed(1)
  X <- data.frame(a = rnorm(30))
  type <- c(a = "numeric")
  prep <- .mahal_prepare(X, names(X), type, robust = FALSE)

  di <- .mahal_dist(X[1, , drop = FALSE], X[1:5, , drop = FALSE], prep, type)
  expect_equal(length(di), 5)
  expect_true(di[1] < 1e-10)
})

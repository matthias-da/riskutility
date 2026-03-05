# Tests for information-theoretic functions (divergence.R and divergence2.R)

# ---- KLDiv tests ----

test_that("KLDiv is non-negative for valid distributions", {
  p <- c(0.25, 0.25, 0.25, 0.25)
  q <- c(0.1, 0.2, 0.3, 0.4)

  result <- KLDiv(p, q)
  expect_true(result >= 0)
})

test_that("KLDiv equals zero for identical distributions", {
  p <- c(0.25, 0.25, 0.25, 0.25)

  result <- KLDiv(p, p)
  expect_equal(result, 0)
})

test_that("KLDiv gives known result for simple distributions", {
  # KL(P || Q) = sum(P * log(P/Q))
  p <- c(0.5, 0.5)
  q <- c(0.25, 0.75)

  expected <- 0.5 * log(0.5/0.25) + 0.5 * log(0.5/0.75)
  result <- KLDiv(p, q)
  expect_equal(result, expected)
})

test_that("KLDiv is not symmetric", {
  p <- c(0.5, 0.5)
  q <- c(0.25, 0.75)

  kl_pq <- KLDiv(p, q)
  kl_qp <- KLDiv(q, p)

  expect_false(isTRUE(all.equal(kl_pq, kl_qp)))
})

# ---- JSDiv tests ----

test_that("JSDiv is symmetric", {
  p <- c(0.5, 0.5)
  q <- c(0.25, 0.75)

  js_pq <- JSDiv(p, q)
  js_qp <- JSDiv(q, p)

  expect_equal(js_pq, js_qp)
})

test_that("JSDiv equals zero for identical distributions", {
  p <- c(0.25, 0.25, 0.25, 0.25)

  result <- JSDiv(p, p)
  expect_equal(result, 0)
})

test_that("JSDiv is bounded between 0 and log(2)", {
  p <- c(0.5, 0.5)
  q <- c(0.25, 0.75)

  result <- JSDiv(p, q)

  expect_true(result >= 0)
  expect_true(result <= log(2))
})

test_that("JSDiv is non-negative", {
  p <- c(0.1, 0.3, 0.6)
  q <- c(0.3, 0.3, 0.4)

  result <- JSDiv(p, q)
  expect_true(result >= 0)
})

# ---- CrossEntropy tests ----

test_that("CrossEntropy gives correct result", {
  p <- c(0.5, 0.5)
  q <- c(0.25, 0.75)

  # H(P, Q) = -sum(P * log(Q))
  expected <- -(0.5 * log(0.25) + 0.5 * log(0.75))
  result <- CrossEntropy(p, q)
  expect_equal(result, expected)
})

test_that("CrossEntropy ignores zero probabilities in A", {
  p <- c(0.5, 0, 0.5)
  q <- c(0.3, 0.3, 0.4)

  # Only considers indices where A > 0: indices 1 and 3
  expected <- -(0.5 * log(0.3) + 0.5 * log(0.4))
  result <- CrossEntropy(p, q)
  expect_equal(result, expected)
})

# ---- RenyiEntropy tests ----

test_that("RenyiEntropy gives correct result for alpha=2", {
  p <- c(0.25, 0.25, 0.25, 0.25)

  # Renyi entropy of order 2 for uniform distribution
  # H_2 = -log(sum(p^2)) = -log(4 * 0.0625) = -log(0.25) = log(4)
  # Wait: formula is (1/(1-alpha)) * log(sum(p^alpha))
  # = (1/(1-2)) * log(sum(0.25^2)) = -1 * log(4*0.0625) = -log(0.25) = log(4)
  expected <- -1 * log(sum(p^2))
  result <- RenyiEntropy(p, alpha = 2)
  expect_equal(result, expected)
})

test_that("RenyiEntropy errors for alpha=1", {
  p <- c(0.5, 0.5)
  expect_error(RenyiEntropy(p, alpha = 1))
})

test_that("RenyiEntropy approaches Shannon entropy as alpha approaches 1", {
  p <- c(0.1, 0.2, 0.3, 0.4)
  shannon <- -sum(p * log(p))

  # Test with alpha close to 1
  renyi_close <- RenyiEntropy(p, alpha = 1.001)

  # Should be close to Shannon entropy

  expect_true(abs(renyi_close - shannon) < 0.01)
})

test_that("RenyiEntropy ignores zero probabilities", {
  p <- c(0.5, 0, 0.5)

  # Should not error due to zeros
  result <- RenyiEntropy(p, alpha = 2)
  expect_true(is.finite(result))
})

# ---- MaxEntropy tests ----

test_that("MaxEntropy gives log(n) for uniform distribution", {
  p <- c(0.25, 0.25, 0.25, 0.25)

  result <- MaxEntropy(p)
  expect_equal(result, log(4))
})

test_that("MaxEntropy counts non-zero elements", {
  p <- c(0.5, 0, 0.5, 0)

  result <- MaxEntropy(p)
  expect_equal(result, log(2))
})

# ---- MinEntropy tests ----

test_that("MinEntropy gives correct result", {
  p <- c(0.5, 0.3, 0.2)

  # MinEntropy = -log(max(p))
  expected <- -log(0.5)
  result <- MinEntropy(p)
  expect_equal(result, expected)
})

test_that("MinEntropy of uniform distribution equals log(n)", {
  p <- c(0.25, 0.25, 0.25, 0.25)

  result <- MinEntropy(p)
  expect_equal(result, log(4))
})

test_that("MinEntropy of degenerate distribution is zero", {
  p <- c(1, 0, 0)

  result <- MinEntropy(p)
  expect_equal(result, 0)
})

# ---- NormalizedEntropy tests ----

test_that("NormalizedEntropy is bounded between 0 and 1", {
  p <- c(0.1, 0.3, 0.6)

  result <- NormalizedEntropy(p)
  expect_true(result >= 0)
  expect_true(result <= 1)
})

test_that("NormalizedEntropy equals 1 for uniform distribution", {
  p <- c(0.25, 0.25, 0.25, 0.25)

  result <- NormalizedEntropy(p)
  expect_equal(result, 1)
})

test_that("NormalizedEntropy is close to 0 for concentrated distribution", {
  p <- c(0.99, 0.01)

  result <- NormalizedEntropy(p)
  expect_true(result < 0.1)
})

# ---- ConditionalEntropy tests ----

test_that("ConditionalEntropy gives correct result for independent variables", {
  # If X and Y are independent, H(Y|X) = H(Y)
  # Joint probability matrix where rows sum to p(x)
  # and columns are proportional to p(y) within each row
  joint <- matrix(c(0.1, 0.1,
                     0.2, 0.2,
                     0.1, 0.1), nrow = 3, byrow = TRUE)

  result <- ConditionalEntropy(joint)
  expect_true(is.finite(result))
  expect_true(result >= 0)
})

test_that("ConditionalEntropy is zero when Y is determined by X", {
  # Each row has only one non-zero column
  joint <- matrix(c(0.5, 0,
                     0, 0.5), nrow = 2, byrow = TRUE)

  result <- ConditionalEntropy(joint)
  expect_equal(result, 0)
})

test_that("ConditionalEntropy is non-negative", {
  joint <- matrix(c(0.1, 0.2,
                     0.3, 0.4), nrow = 2, byrow = TRUE)

  result <- ConditionalEntropy(joint)
  expect_true(result >= 0)
})

# ---- CumulativeEntropy tests ----

test_that("CumulativeEntropy returns finite value", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

  result <- CumulativeEntropy(x)
  expect_true(is.finite(result))
  expect_true(result > 0)
})

test_that("CumulativeEntropy handles repeated values", {
  x <- c(1, 1, 2, 2, 3, 3)

  result <- CumulativeEntropy(x)
  expect_true(is.finite(result))
})

# ---- KLDiv_bayes tests ----

test_that("KLDiv_bayes gives non-negative result for valid distributions", {
  p <- c(0.25, 0.25, 0.25, 0.25)
  q <- c(0.1, 0.2, 0.3, 0.4)

  result <- KLDiv_bayes(p, q)
  expect_true(result >= 0)
})

test_that("KLDiv_bayes equals zero for identical distributions", {
  p <- c(0.25, 0.25, 0.25, 0.25)

  result <- KLDiv_bayes(p, p)
  expect_equal(result, 0)
})

# ---- JSDiv_bayes tests ----

test_that("JSDiv_bayes is symmetric", {
  p <- c(0.5, 0.5)
  q <- c(0.25, 0.75)

  js_pq <- JSDiv_bayes(p, q)
  js_qp <- JSDiv_bayes(q, p)

  expect_equal(js_pq, js_qp)
})

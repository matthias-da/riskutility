# tests/testthat/test-contingency-fidelity.R

test_that("contingency_fidelity returns correct S3 class", {
set.seed(1)
n <- 200
X <- data.frame(a = sample(c("x", "y"), n, TRUE),
                b = sample(c("p", "q"), n, TRUE),
                c = sample(c("1", "2", "3"), n, TRUE),
                stringsAsFactors = TRUE)
Y <- X  # identical copy
res <- contingency_fidelity(X, Y)
expect_s3_class(res, "contingency_fidelity")
expect_true(all(c("mean_tv", "utility_score", "pairwise",
                   "n_vars", "vars", "n_X", "n_Y") %in% names(res)))
})

test_that("near-zero TV for identical distributions", {
  set.seed(2)
  n <- 500
  X <- data.frame(a = sample(letters[1:3], n, TRUE),
                  b = sample(LETTERS[1:4], n, TRUE),
                  stringsAsFactors = TRUE)
  res <- contingency_fidelity(X, X)
  expect_equal(res$mean_tv, 0, tolerance = 1e-10)
  expect_equal(res$utility_score, 1, tolerance = 1e-10)
})

test_that("detects different distributions (TV > 0)", {
  set.seed(3)
  n <- 500
  X <- data.frame(a = sample(c("M", "F"), n, TRUE),
                  b = sample(c("low", "high"), n, TRUE),
                  stringsAsFactors = TRUE)
  # Synthetic has completely different joint distribution
  Y <- data.frame(a = sample(c("M", "F"), n, TRUE, prob = c(0.9, 0.1)),
                  b = sample(c("low", "high"), n, TRUE, prob = c(0.1, 0.9)),
                  stringsAsFactors = TRUE)
  res <- contingency_fidelity(X, Y)
  expect_gt(res$mean_tv, 0.05)
  expect_lt(res$utility_score, 1)
})

test_that("correct pairwise dimensions: 3 vars -> 3 pairs", {
  set.seed(4)
  n <- 100
  X <- data.frame(a = sample(c("x", "y"), n, TRUE),
                  b = sample(c("p", "q"), n, TRUE),
                  c = sample(c("1", "2"), n, TRUE),
                  stringsAsFactors = TRUE)
  res <- contingency_fidelity(X, X)
  expect_equal(nrow(res$pairwise), 3)
  expect_equal(res$n_vars, 3)
})

test_that("correct pairwise dimensions: 4 vars -> 6 pairs", {
  set.seed(5)
  n <- 100
  X <- data.frame(a = sample(c("x", "y"), n, TRUE),
                  b = sample(c("p", "q"), n, TRUE),
                  c = sample(c("1", "2"), n, TRUE),
                  d = sample(c("A", "B", "C"), n, TRUE),
                  stringsAsFactors = TRUE)
  res <- contingency_fidelity(X, X)
  expect_equal(nrow(res$pairwise), 6)
  expect_equal(res$n_vars, 4)
})

test_that("skips numeric variables with message", {
  set.seed(6)
  n <- 100
  X <- data.frame(a = sample(c("x", "y"), n, TRUE),
                  b = sample(c("p", "q"), n, TRUE),
                  c = sample(c("1", "2"), n, TRUE),
                  num1 = rnorm(n),
                  stringsAsFactors = TRUE)
  Y <- X
  expect_message(contingency_fidelity(X, Y), "Skipping non-categorical")
})

test_that("utility_score is in [0, 1]", {
  set.seed(7)
  n <- 200
  X <- data.frame(a = sample(c("M", "F"), n, TRUE),
                  b = sample(c("Y", "N"), n, TRUE),
                  stringsAsFactors = TRUE)
  Y <- data.frame(a = sample(c("M", "F"), n, TRUE, prob = c(0.8, 0.2)),
                  b = sample(c("Y", "N"), n, TRUE, prob = c(0.3, 0.7)),
                  stringsAsFactors = TRUE)
  res <- contingency_fidelity(X, Y)
  expect_gte(res$utility_score, 0)
  expect_lte(res$utility_score, 1)
})

test_that("TV distance hand-computed for simple 2x2 table", {
  # Original: a=A,b=X appears 50%; a=A,b=Y 0%; a=B,b=X 0%; a=B,b=Y 50%
  # Synthetic: a=A,b=X 25%; a=A,b=Y 25%; a=B,b=X 25%; a=B,b=Y 25%
  # TV = 0.5 * (|0.5-0.25| + |0-0.25| + |0-0.25| + |0.5-0.25|)
  #    = 0.5 * (0.25 + 0.25 + 0.25 + 0.25) = 0.5
  X <- data.frame(a = factor(c("A", "A", "B", "B"), levels = c("A", "B")),
                  b = factor(c("X", "X", "Y", "Y"), levels = c("X", "Y")))
  Y <- data.frame(a = factor(c("A", "A", "B", "B"), levels = c("A", "B")),
                  b = factor(c("X", "Y", "X", "Y"), levels = c("X", "Y")))
  res <- contingency_fidelity(X, Y)
  expect_equal(res$mean_tv, 0.5, tolerance = 1e-10)
  expect_equal(res$utility_score, 0.5, tolerance = 1e-10)
})

test_that("TV distance hand-computed: identical uniform 2x2", {
  # Both uniform: TV should be 0
  X <- data.frame(a = factor(c("A", "A", "B", "B")),
                  b = factor(c("X", "Y", "X", "Y")))
  Y <- data.frame(a = factor(c("A", "A", "B", "B")),
                  b = factor(c("X", "Y", "X", "Y")))
  res <- contingency_fidelity(X, Y)
  expect_equal(res$mean_tv, 0)
})

test_that("handles mismatched factor levels (level in orig not in synth)", {
  # Original has level "C" that synthetic doesn't have
  X <- data.frame(a = factor(c("A", "B", "C", "A", "B", "C"),
                             levels = c("A", "B", "C")),
                  b = factor(c("X", "Y", "X", "Y", "X", "Y")))
  Y <- data.frame(a = factor(c("A", "B", "A", "B", "A", "B"),
                             levels = c("A", "B")),
                  b = factor(c("X", "Y", "X", "Y", "X", "Y")))
  # Should not error: missing level "C" in Y gets proportion 0
  res <- contingency_fidelity(X, Y)
  expect_s3_class(res, "contingency_fidelity")
  expect_gt(res$mean_tv, 0)
})

test_that("handles level in synth not in orig", {
  X <- data.frame(a = factor(c("A", "A", "B", "B")),
                  b = factor(c("X", "Y", "X", "Y")))
  Y <- data.frame(a = factor(c("A", "A", "C", "C")),
                  b = factor(c("X", "Y", "X", "Y")))
  res <- contingency_fidelity(X, Y)
  expect_s3_class(res, "contingency_fidelity")
  expect_gt(res$mean_tv, 0)
})

test_that("synth_pair dispatch works", {
  set.seed(8)
  n <- 100
  X <- data.frame(a = sample(c("x", "y"), n, TRUE),
                  b = sample(c("p", "q"), n, TRUE),
                  c = sample(c("1", "2"), n, TRUE),
                  stringsAsFactors = TRUE)
  sp <- list(original = X, synthetic = X,
             cat_vars = c("a", "b", "c"))
  class(sp) <- "synth_pair"
  res <- contingency_fidelity(sp)
  expect_s3_class(res, "contingency_fidelity")
  expect_equal(res$mean_tv, 0, tolerance = 1e-10)
})

test_that("vars parameter selects specific variables", {
  set.seed(9)
  n <- 100
  X <- data.frame(a = sample(c("x", "y"), n, TRUE),
                  b = sample(c("p", "q"), n, TRUE),
                  c = sample(c("1", "2"), n, TRUE),
                  stringsAsFactors = TRUE)
  res <- contingency_fidelity(X, X, vars = c("a", "b"))
  expect_equal(res$n_vars, 2)
  expect_equal(nrow(res$pairwise), 1)
  expect_equal(res$vars, c("a", "b"))
})

test_that("error with fewer than 2 categorical variables", {
  X <- data.frame(a = factor(c("x", "y", "x")))
  Y <- data.frame(a = factor(c("x", "y", "x")))
  expect_error(contingency_fidelity(X, Y), "At least 2 categorical")
})

test_that("error with only numeric variables", {
  X <- data.frame(a = 1:10, b = rnorm(10))
  Y <- data.frame(a = 1:10, b = rnorm(10))
  expect_error(contingency_fidelity(X, Y), "At least 2 categorical")
})

test_that("print method works without error", {
  set.seed(10)
  n <- 50
  X <- data.frame(a = sample(c("x", "y"), n, TRUE),
                  b = sample(c("p", "q"), n, TRUE),
                  stringsAsFactors = TRUE)
  res <- contingency_fidelity(X, X)
  expect_output(print(res), "Contingency Fidelity")
  expect_output(print(res), "Utility score")
})

test_that("summary method returns summary object", {
  set.seed(11)
  n <- 100
  X <- data.frame(a = sample(c("x", "y"), n, TRUE),
                  b = sample(c("p", "q"), n, TRUE),
                  c = sample(c("1", "2"), n, TRUE),
                  stringsAsFactors = TRUE)
  res <- contingency_fidelity(X, X)
  s <- summary(res)
  expect_s3_class(s, "summary.contingency_fidelity")
  expect_true(all(c("mean_tv", "utility_score", "pairwise",
                     "max_tv", "min_tv", "sd_tv", "n_pairs") %in% names(s)))
  expect_output(print(s), "Summary: Contingency Fidelity")
})

test_that("plot method produces ggplot without error", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("reshape2")
  set.seed(12)
  n <- 100
  X <- data.frame(a = sample(c("x", "y"), n, TRUE),
                  b = sample(c("p", "q"), n, TRUE),
                  c = sample(c("1", "2"), n, TRUE),
                  stringsAsFactors = TRUE)
  Y <- data.frame(a = sample(c("x", "y"), n, TRUE),
                  b = sample(c("p", "q"), n, TRUE, prob = c(0.8, 0.2)),
                  c = sample(c("1", "2"), n, TRUE),
                  stringsAsFactors = TRUE)
  res <- contingency_fidelity(X, Y)
  # Should not error
  expect_invisible(plot(res, which = 1))
})

test_that("NA handling removes incomplete rows", {
  X <- data.frame(a = factor(c("A", "B", NA, "A", "B")),
                  b = factor(c("X", "Y", "X", NA, "Y")))
  Y <- data.frame(a = factor(c("A", "B", "A", "B", "A")),
                  b = factor(c("X", "Y", "X", "Y", "X")))
  res <- contingency_fidelity(X, Y, na.rm = TRUE)
  expect_equal(res$n_X, 3)  # rows 1, 2, 5 survive
  expect_s3_class(res, "contingency_fidelity")
})

test_that("character columns are handled like factors", {
  X <- data.frame(a = c("x", "y", "x", "y"),
                  b = c("p", "q", "p", "q"),
                  stringsAsFactors = FALSE)
  Y <- data.frame(a = c("x", "y", "x", "y"),
                  b = c("p", "q", "p", "q"),
                  stringsAsFactors = FALSE)
  res <- contingency_fidelity(X, Y)
  expect_s3_class(res, "contingency_fidelity")
  expect_equal(res$mean_tv, 0)
})

test_that("TV distance is symmetric", {
  set.seed(13)
  n <- 200
  X <- data.frame(a = sample(c("M", "F"), n, TRUE),
                  b = sample(c("Y", "N"), n, TRUE),
                  stringsAsFactors = TRUE)
  Y <- data.frame(a = sample(c("M", "F"), n, TRUE, prob = c(0.7, 0.3)),
                  b = sample(c("Y", "N"), n, TRUE, prob = c(0.4, 0.6)),
                  stringsAsFactors = TRUE)
  res_xy <- contingency_fidelity(X, Y)
  res_yx <- contingency_fidelity(Y, X)
  expect_equal(res_xy$mean_tv, res_yx$mean_tv, tolerance = 1e-10)
})

test_that("completely disjoint distributions give TV = 1", {
  # X has only (A, X) and (B, Y); Y has only (A, Y) and (B, X)
  X <- data.frame(a = factor(c("A", "A", "B", "B"), levels = c("A", "B")),
                  b = factor(c("X", "X", "Y", "Y"), levels = c("X", "Y")))
  Y <- data.frame(a = factor(c("A", "A", "B", "B"), levels = c("A", "B")),
                  b = factor(c("Y", "Y", "X", "X"), levels = c("X", "Y")))
  res <- contingency_fidelity(X, Y)
  expect_equal(res$mean_tv, 1, tolerance = 1e-10)
  expect_equal(res$utility_score, 0, tolerance = 1e-10)
})

test_that("print shows interpretation bands correctly", {
  # Excellent band (identical data)
  set.seed(14)
  n <- 100
  X <- data.frame(a = sample(c("x", "y"), n, TRUE),
                  b = sample(c("p", "q"), n, TRUE),
                  stringsAsFactors = TRUE)
  res <- contingency_fidelity(X, X)
  expect_output(print(res), "EXCELLENT")
})

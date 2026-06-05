# Unit tests for the shared .true_match_risk() back-end (used by all 8 methods,
# and the new common path for rf/embedding).

library(testthat)

test_that(".true_match_risk: maximize, nearest != true -> risk 0", {
  # similarity scores; candidate 30 is the closest (0.9) but the TRUE match is
  # candidate 20 (0.2). Under nearest, the true match is not the nearest -> 0.
  r <- riskutility:::.true_match_risk(
    scores = c(0.1, 0.2, 0.9), cand = c(10L, 20L, 30L),
    true_pos = 20L, maximize = TRUE,
    strategy = "nearest", risk_weighting = "uniform")
  expect_equal(r$risk, 0)
  expect_false(r$true_in_set)
})

test_that(".true_match_risk: topk uniform gives 1/k when true in set", {
  r <- riskutility:::.true_match_risk(
    scores = c(0.1, 0.2, 0.9), cand = c(10L, 20L, 30L),
    true_pos = 20L, maximize = TRUE,
    strategy = "topk", k = 3, risk_weighting = "uniform")
  expect_equal(r$cand_n, 3L)
  expect_true(r$true_in_set)
  expect_equal(r$risk, 1 / 3)
})

test_that(".true_match_risk: distance, true is the unique nearest -> risk 1", {
  r <- riskutility:::.true_match_risk(
    scores = c(0.1, 0.2, 0.9), cand = c(10L, 20L, 30L),
    true_pos = 10L, maximize = FALSE,
    strategy = "nearest", risk_weighting = "uniform")
  expect_equal(r$cand_n, 1L)
  expect_true(r$true_in_set)
  expect_equal(r$risk, 1)
})

test_that(".true_match_risk: empty candidate set -> risk 0", {
  r <- riskutility:::.true_match_risk(
    scores = numeric(0), cand = integer(0), true_pos = NA_integer_,
    maximize = TRUE, strategy = "nearest", risk_weighting = "uniform")
  expect_equal(r$risk, 0)
  expect_equal(r$cand_n, 0L)
  expect_false(r$true_in_set)
})

test_that(".true_match_risk: softmax weights sum semantics (true gets its share)", {
  # maximize: higher score = closer. With topk k=3 (all), softmax over pseudo-d.
  r <- riskutility:::.true_match_risk(
    scores = c(2, 2, 2), cand = c(1L, 2L, 3L), true_pos = 2L,
    maximize = TRUE, strategy = "topk", k = 3, risk_weighting = "softmax")
  # all equal -> uniform softmax -> 1/3 each
  expect_equal(r$risk, 1 / 3)
})

test_that("probabilistic: user m/u outside (0,1) are clamped with a warning", {
  set.seed(5)
  x <- data.frame(a = round(rnorm(8, 10, 3)))
  res <- NULL
  ws <- capture_warnings(
    res <- recordLinkage(x, x, key = "a", method = "probabilistic",
                         truth = "row", m_probs = c(a = 1), u_probs = c(a = 0))
  )
  expect_true(any(grepl("clamped", ws)))      # both m_probs and u_probs warn
  expect_true(all(is.finite(res$per_record$risk)))
})

test_that("probabilistic: unnamed m_probs is rejected", {
  x <- data.frame(a = 1:5, b = factor(letters[1:5]))
  expect_error(
    recordLinkage(x, x, key = c("a", "b"), method = "probabilistic",
                  truth = "row", m_probs = c(0.9, 0.9),
                  u_probs = c(a = 0.3, b = 0.3)),
    "named"
  )
})

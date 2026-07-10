library(testthat)

test_that("probabilistic: user m/u outside (0,1) are clamped with a warning", {
  set.seed(5)
  x <- data.frame(a = round(rnorm(8, 10, 3)))
  res <- NULL
  ws <- capture_warnings(
    res <- recordLinkage(x, x, key = "a", method = "probabilistic",
                         truth = "row",
                         control = rl_control(m_probs = c(a = 1),
                                              u_probs = c(a = 0)))
  )
  expect_true(any(grepl("clamped", ws)))
  expect_true(all(is.finite(res$per_record$risk)))
})

test_that("probabilistic: unnamed m_probs is rejected", {
  x <- data.frame(a = 1:5, b = factor(letters[1:5]))
  expect_error(
    recordLinkage(x, x, key = c("a", "b"), method = "probabilistic",
                  truth = "row",
                  control = rl_control(m_probs = c(0.9, 0.9),
                                       u_probs = c(a = 0.3, b = 0.3))),
    "named"
  )
})

test_that("probabilistic: direction = 'reverse' runs without warning", {
  x <- data.frame(a = 1:5, b = factor(letters[1:5]))
  res <- expect_no_warning(
    recordLinkage(x, x, key = c("a", "b"), method = "probabilistic",
                  direction = "anon_to_original")
  )
  expect_equal(res$direction, "anon_to_original")
})

test_that("probabilistic: na_anon != 'ignore' in rl_control() runs without warning", {
  x <- data.frame(a = 1:5, b = factor(letters[1:5]))
  expect_no_warning(
    recordLinkage(x, x, key = c("a", "b"), method = "probabilistic",
                  control = rl_control(na_anon = "match"))
  )
})

test_that("probabilistic: cand_n equals total block size, not positive-LR count", {
  set.seed(1)
  x <- data.frame(
    region = factor(rep(c("A", "B"), each = 5)),
    age    = c(20, 25, 30, 35, 40, 45, 50, 55, 60, 65)
  )
  res <- recordLinkage(x, x, key = c("region", "age"),
                       method = "probabilistic",
                       block = "region")
  # Each block has 5 records; cand_n should be 5 for every record
  expect_true(all(res$per_record$cand_n == 5L))
})

test_that(".choose_guess_set: maximize picks candidates at the maximum (nearest)", {
  d <- c(0.1, 0.9, 0.9, 0.3)
  cand <- c(10, 20, 30, 40)
  out <- riskutility:::.choose_guess_set(d, cand, strategy = "nearest",
                                         maximize = TRUE)
  expect_equal(sort(out), c(20, 30))
})

test_that(".choose_guess_set: maximize with topk picks the k largest, ties expand", {
  d <- c(5, 1, 5, 3, 2)
  cand <- 1:5
  out <- riskutility:::.choose_guess_set(d, cand, strategy = "topk", k = 2,
                                         maximize = TRUE)
  # Two largest values (5, 5) belong to cand 1 and 3
  expect_equal(sort(out), c(1, 3))
})

test_that(".choose_guess_set: maximize = FALSE (default) still picks the minimum", {
  d <- c(0.1, 0.9, 0.9, 0.3)
  cand <- c(10, 20, 30, 40)
  out <- riskutility:::.choose_guess_set(d, cand, strategy = "nearest")
  expect_equal(out, 10)
})

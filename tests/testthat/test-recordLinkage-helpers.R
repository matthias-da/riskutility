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

test_that("probabilistic: direction != 'forward' warns and is ignored", {
  x <- data.frame(a = 1:5, b = factor(letters[1:5]))
  expect_warning(
    recordLinkage(x, x, key = c("a", "b"), method = "probabilistic",
                  direction = "reverse"),
    "direction"
  )
})

test_that("probabilistic: na_anon != 'ignore' in rl_control() warns", {
  x <- data.frame(a = 1:5, b = factor(letters[1:5]))
  expect_warning(
    recordLinkage(x, x, key = c("a", "b"), method = "probabilistic",
                  control = rl_control(na_anon = "match")),
    "na_anon"
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

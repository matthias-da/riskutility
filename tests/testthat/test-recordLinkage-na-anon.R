# na_anon (ignore / match / mismatch) honored across methods.

library(testthat)

test_that("probabilistic: na_anon changes the evidence on a missing key", {
  set.seed(2)
  n <- 10
  x <- data.frame(a = round(rnorm(n, 10, 3)),
                  b = factor(sample(c("x", "y"), n, TRUE)))
  x_anon <- x
  x_anon$b[1] <- NA

  r_m  <- recordLinkage(x, x_anon, key = c("a", "b"), method = "probabilistic",
                        truth = "row", na_anon = "match")
  r_mm <- recordLinkage(x, x_anon, key = c("a", "b"), method = "probabilistic",
                        truth = "row", na_anon = "mismatch")
  r_ig <- recordLinkage(x, x_anon, key = c("a", "b"), method = "probabilistic",
                        truth = "row", na_anon = "ignore")

  # The true-match log-LR for record 1 differs across the three modes.
  expect_false(isTRUE(all.equal(r_m$per_record$d_true[1],
                                r_mm$per_record$d_true[1])))
  expect_false(isTRUE(all.equal(r_ig$per_record$d_true[1],
                                r_mm$per_record$d_true[1])))
})

test_that("na-free data is unaffected by the na_anon default change", {
  set.seed(4)
  n <- 20
  x <- data.frame(age = round(rnorm(n, 45, 10)),
                  sex = factor(sample(c("f", "m"), n, TRUE)))
  x_anon <- x
  x_anon$age <- x_anon$age + sample(c(-1, 0, 1), n, TRUE)
  # Re-seed before each call: .fs_estimate() samples random non-match pairs
  # for u-prob estimation, so both calls must start from the same RNG state
  # to produce identical LRs (and hence identical posteriors).
  set.seed(42)
  r1 <- recordLinkage(x, x_anon, key = c("age", "sex"),
                      method = "probabilistic", truth = "row")
  set.seed(42)
  r2 <- recordLinkage(x, x_anon, key = c("age", "sex"),
                      method = "probabilistic", truth = "row",
                      na_anon = "ignore")
  expect_equal(r1$per_record$risk, r2$per_record$risk)
})

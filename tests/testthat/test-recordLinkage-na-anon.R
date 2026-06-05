# na_anon (ignore / match / mismatch) is now honored by all methods, not only
# the deterministic Gower path. Key regression: PRAM 'ignore' must NOT collapse
# a record with a missing key to 0 risk (the former -Inf -> prob 0 bug).

library(testthat)

test_that("PRAM: na_anon='ignore' keeps nonzero risk on a missing key", {
  lv <- c("a", "b")
  tm <- matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2, dimnames = list(lv, lv))
  pm <- list(k1 = tm, k2 = tm)
  set.seed(1)
  n <- 6
  x <- data.frame(
    k1 = factor(sample(lv, n, TRUE), levels = lv),
    k2 = factor(sample(lv, n, TRUE), levels = lv)
  )
  x_anon <- x
  x_anon$k1[1] <- NA      # candidate 1 (the true match of record 1) misses k1

  res_ig <- recordLinkage(x, x_anon, key = c("k1", "k2"), method = "pram",
                          pram_matrix = pm, truth = "row", na_anon = "ignore")
  res_mm <- recordLinkage(x, x_anon, key = c("k1", "k2"), method = "pram",
                          pram_matrix = pm, truth = "row", na_anon = "mismatch")

  # ignore: variable k1 drops out, risk comes from k2 -> nonzero
  expect_gt(res_ig$per_record$risk[1], 0)
  # mismatch: the true candidate is excluded -> risk 0
  expect_equal(res_mm$per_record$risk[1], 0)
})

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

test_that("mahalanobis: na_anon honored on a missing nominal key", {
  set.seed(3)
  n <- 12
  x <- data.frame(age = round(rnorm(n, 45, 10)),
                  sex = factor(sample(c("f", "m"), n, TRUE)))
  x_anon <- x
  x_anon$sex[1] <- NA

  r_match <- recordLinkage(x, x_anon, key = c("age", "sex"),
                           method = "mahalanobis", truth = "row",
                           na_anon = "match")
  r_mis   <- recordLinkage(x, x_anon, key = c("age", "sex"),
                           method = "mahalanobis", truth = "row",
                           na_anon = "mismatch")
  # d_true (distance of the true match) for record 1 differs: match treats the
  # NA as agreement (smaller distance) than mismatch.
  expect_false(isTRUE(all.equal(r_match$per_record$d_true[1],
                                r_mis$per_record$d_true[1])))
})

test_that("na-free data is unaffected by the na_anon default change", {
  set.seed(4)
  n <- 20
  x <- data.frame(age = round(rnorm(n, 45, 10)),
                  sex = factor(sample(c("f", "m"), n, TRUE)))
  x_anon <- x
  x_anon$age <- x_anon$age + sample(c(-1, 0, 1), n, TRUE)
  r1 <- recordLinkage(x, x_anon, key = c("age", "sex"),
                      method = "probabilistic", truth = "row")
  r2 <- recordLinkage(x, x_anon, key = c("age", "sex"),
                      method = "probabilistic", truth = "row",
                      na_anon = "ignore")
  expect_equal(r1$per_record$risk, r2$per_record$risk)
})

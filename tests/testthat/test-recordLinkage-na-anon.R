# na_anon (ignore / match / mismatch) honored across methods.

library(testthat)

test_that("probabilistic: na_anon in rl_control() has no effect (NAs always ignored)", {
  set.seed(2)
  n <- 10
  x <- data.frame(a = round(rnorm(n, 10, 3)),
                  b = factor(sample(c("x", "y"), n, TRUE)))
  x_anon <- x
  x_anon$b[1] <- NA

  set.seed(42)
  r_ig <- suppressWarnings(
    recordLinkage(x, x_anon, key = c("a", "b"), method = "probabilistic",
                  truth = "row", control = rl_control(na_anon = "ignore"))
  )
  set.seed(42)
  r_mm <- suppressWarnings(
    recordLinkage(x, x_anon, key = c("a", "b"), method = "probabilistic",
                  truth = "row", control = rl_control(na_anon = "mismatch"))
  )

  # Probabilistic always drops NA pairs regardless of na_anon setting.
  expect_equal(r_ig$per_record$lr_true, r_mm$per_record$lr_true)
})

test_that("na-free data is unaffected by the na_anon setting", {
  set.seed(4)
  n <- 20
  x <- data.frame(age = round(rnorm(n, 45, 10)),
                  sex = factor(sample(c("f", "m"), n, TRUE)))
  x_anon <- x
  x_anon$age <- x_anon$age + sample(c(-1, 0, 1), n, TRUE)
  set.seed(42)
  r1 <- recordLinkage(x, x_anon, key = c("age", "sex"),
                      method = "probabilistic", truth = "row")
  set.seed(42)
  r2 <- recordLinkage(x, x_anon, key = c("age", "sex"),
                      method = "probabilistic", truth = "row",
                      control = rl_control(na_anon = "ignore"))
  expect_equal(r1$per_record$risk, r2$per_record$risk)
})

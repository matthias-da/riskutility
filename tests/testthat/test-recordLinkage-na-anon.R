# na_anon (ignore / match / mismatch) honored across methods.

library(testthat)

test_that("probabilistic: na_anon modes give different log-LR for NA pairs", {
  set.seed(2)
  n <- 10
  x <- data.frame(a = round(rnorm(n, 10, 3)),
                  b = factor(sample(c("x", "y"), n, TRUE)))
  x_anon <- x
  x_anon$b[1] <- NA

  set.seed(42)
  r_ig <- recordLinkage(x, x_anon, key = c("a", "b"), method = "probabilistic",
                        truth = "row", control = rl_control(na_anon = "ignore"))
  set.seed(42)
  r_ma <- recordLinkage(x, x_anon, key = c("a", "b"), method = "probabilistic",
                        truth = "row", control = rl_control(na_anon = "match"))
  set.seed(42)
  r_mm <- recordLinkage(x, x_anon, key = c("a", "b"), method = "probabilistic",
                        truth = "row", control = rl_control(na_anon = "mismatch"))

  m_probs <- r_ig$fs_params$m_probs["b"]
  u_probs <- r_ig$fs_params$u_probs["b"]

  # Row 1 has an NA in 'b'; the difference between modes for that row's
  # log-LR against itself (truth = "row") is exactly the assumed contribution
  # of variable 'b': log(m/u) for "match", log((1-m)/(1-u)) for "mismatch",
  # 0 for "ignore".
  expect_equal(r_ma$per_record$lr_true[1] - r_ig$per_record$lr_true[1],
               unname(log(m_probs / u_probs)))
  expect_equal(r_mm$per_record$lr_true[1] - r_ig$per_record$lr_true[1],
               unname(log((1 - m_probs) / (1 - u_probs))))

  # Rows without any NA are unaffected by na_anon.
  expect_equal(r_ig$per_record$lr_true[-1], r_ma$per_record$lr_true[-1])
  expect_equal(r_ig$per_record$lr_true[-1], r_mm$per_record$lr_true[-1])
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

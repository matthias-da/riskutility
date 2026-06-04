# Tests for backward-compatible / consistency argument aliases introduced in the
# Tier-4 API polish:
#   * compare_distributions_cont(): 'vars' is canonical; 'variables' still accepted.
#   * mia_classifier():            'Y'/'holdout'    alias 'synt_data'/'hout_data'.
#   * recordLinkage():             'Y'/'key_vars'   alias 'x_anon'/'key'.
#   * rumap():                     'Y'              aliases 'synthetic'.

library(testthat)

test_that("compare_distributions_cont: vars, variables (legacy) and partial match agree", {
  set.seed(1)
  X <- data.frame(age = rnorm(80, 40, 10), income = rnorm(80, 5e4, 1e4))
  Y <- data.frame(age = rnorm(80, 40, 10), income = rnorm(80, 5e4, 1e4))

  r_vars <- suppressWarnings(compare_distributions_cont(X, Y, vars = "income", n_approx = 100))
  r_old  <- suppressWarnings(compare_distributions_cont(X, Y, variables = "income", n_approx = 100))
  r_pos  <- suppressWarnings(compare_distributions_cont(X, Y, "income", n_approx = 100))
  r_part <- suppressWarnings(compare_distributions_cont(X, Y, var = "income", n_approx = 100))

  expect_s3_class(r_vars, "compare_distributions_cont")
  expect_equal(r_vars$ecdf, r_old$ecdf)
  expect_equal(r_vars$ecdf, r_pos$ecdf)
  expect_equal(r_vars$ecdf, r_part$ecdf)
})

test_that("mia_classifier: Y/holdout aliases match synt_data/hout_data", {
  skip_if_not_installed("ranger")
  set.seed(42)
  real <- data.frame(x1 = rnorm(100), x2 = factor(sample(c("A", "B"), 100, TRUE)))
  synt <- data.frame(x1 = rnorm(100), x2 = factor(sample(c("A", "B"), 100, TRUE)))
  hout <- data.frame(x1 = rnorm(100), x2 = factor(sample(c("A", "B"), 100, TRUE)))

  r_old <- mia_classifier(real, synt_data = synt, hout_data = hout,
                          num_eval_iter = 2, seed = 7)
  r_new <- mia_classifier(real, Y = synt, holdout = hout,
                          num_eval_iter = 2, seed = 7)

  expect_s3_class(r_new, "mia")
  expect_equal(r_new$recall, r_old$recall)
  expect_equal(r_new$precision, r_old$precision)
  expect_equal(r_new$macro_f1, r_old$macro_f1)
})

test_that("mia_classifier errors helpfully when synthetic/holdout missing", {
  real <- data.frame(x1 = rnorm(10))
  expect_error(mia_classifier(real), "Y.*synt_data|synt_data.*Y")
})

test_that("recordLinkage: Y/key_vars aliases match x_anon/key", {
  set.seed(1)
  x <- data.frame(
    age    = sample(18:80, 100, TRUE),
    sex    = factor(sample(c("f", "m"), 100, TRUE)),
    region = factor(sample(paste0("R", 1:5), 100, TRUE))
  )
  x_anon <- x
  for (r in levels(x$region)) {
    idx <- which(x$region == r)
    if (length(idx) > 1L) x_anon$age[idx] <- sample(x$age[idx])
  }
  keys <- c("age", "sex", "region")

  r_old <- recordLinkage(x, x_anon = x_anon, key = keys)
  r_new <- recordLinkage(x, Y = x_anon, key_vars = keys)

  expect_s3_class(r_new, "recordLinkageRisk")
  expect_equal(r_new$per_record, r_old$per_record)
})

test_that("recordLinkage errors helpfully when comparison data / key missing", {
  x <- data.frame(age = 1:5, sex = factor(rep("f", 5)))
  expect_error(recordLinkage(x), "Y.*x_anon|x_anon.*Y")
})

test_that("rumap: Y alias matches synthetic", {
  set.seed(1)
  orig  <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
  synth <- data.frame(x1 = rnorm(100), x2 = rnorm(100))

  r_old <- rumap(orig, synthetic = list(s1 = synth),
                 risk_measures = "ims",
                 utility_measures = c("hellinger", "energy"),
                 seed = 5)
  r_new <- rumap(orig, Y = list(s1 = synth),
                 risk_measures = "ims",
                 utility_measures = c("hellinger", "energy"),
                 seed = 5)

  expect_s3_class(r_new, "rumap")
  expect_equal(r_new$risk, r_old$risk)
  expect_equal(r_new$utility, r_old$utility)
})

# Parity with Oscar Thees' RAPID package (RAPID_Dev @ f436359).
# Reference values: fixtures/rapid_ot_fixture.rds, built by
# fixtures/make-rapid-ot-fixture.R.

fx <- rapid_fixture()

expect_band_parity <- function(m) {
  y <- fx$num$o$y
  rule <- .rapid_rule("band", 0.3, "symmetric", 0.1, 0.01)
  sc <- .rule_score(rule, list(point = fx$point[[m]], dist = .pdist_grid(fx$D[[m]])), y, y)
  ref <- fx$band[[m]]
  expect_equal(sc$band_lo, ref$band_lo, tolerance = 1e-12, info = m)
  expect_equal(sc$band_hi, ref$band_hi, tolerance = 1e-12, info = m)
  expect_equal(sc$g, ref$g_hat, tolerance = 1e-12, info = m)
  expect_equal(sc$b, ref$b_hat, tolerance = 1e-12, info = m)
  expect_equal(sc$r, ref$r_hat, tolerance = 1e-12, info = m)
  expect_equal(sc$predicted_value, fx$point[[m]], tolerance = 1e-12, info = m)
}

test_that("scoring reproduces Oscar's band_risk for every model (lm, cart)", {
  for (m in c("lm", "cart")) expect_band_parity(m)
})

test_that("scoring reproduces Oscar's band_risk for every model (rf, gbm)", {
  skip_on_cran()
  for (m in c("rf", "gbm")) expect_band_parity(m)
})

test_that("bounds reproduce Oscar's band_bounds", {
  for (metric in unique(fx$bounds$metric)) {
    ref <- fx$bounds[fx$bounds$metric == metric, ]
    b <- .rapid_band_bounds(ref$y, metric, ref$epsilon[1], ref$delta[1])
    expect_equal(b$lo, ref$lo, tolerance = 1e-12, info = metric)
    expect_equal(b$hi, ref$hi, tolerance = 1e-12, info = metric)
  }
})

test_that("RCS_marginal reproduces Oscar's normalized gain", {
  skip_on_cran()  # the stored class probabilities come from an rf model
  y <- fx$cat$o$y
  rule <- .rapid_rule("RCS_marginal", 0.3)
  sc <- .rule_score(rule, fx$cat$P, y, y)
  expect_equal(sc$r, fx$cat$gain, tolerance = 1e-12)
  expect_equal(mean(.rule_flag(rule, sc$r)), fx$cat$rapid, tolerance = 1e-12)
})

# ---- rapid() as a whole ------------------------------------------------------

qi <- c("x1", "x2")
# rpart's own defaults, as used by Oscar's package (minbucket = 7 alone would
# make rpart set minsplit = 21)
rpart_defaults <- list(prune = FALSE, cp = 0.01, minsplit = 20, minbucket = 7)

test_that("rapid() reproduces Oscar's RAPID for lm and cart", {
  skip_if_not_installed("rpart")
  lm_ <- rapid(fx$num$o, fx$num$s, qi, "y", model_type = "lm", num_dist = "grid",
               num_scale = "identity")
  expect_equal(lm_$rapid, fx$rapid$lm, tolerance = 1e-8)
  expect_equal(lm_$records$r, fx$band$lm$r_hat, tolerance = 1e-8)
  expect_equal(lm_$records$g, fx$band$lm$g_hat, tolerance = 1e-8)

  cart <- do.call(rapid, c(list(fx$num$o, fx$num$s, qi, "y", model_type = "cart",
                                num_dist = "grid"), rpart_defaults))
  expect_equal(cart$rapid, fx$rapid$cart, tolerance = 1e-8)
  expect_equal(cart$records$r, fx$band$cart$r_hat, tolerance = 1e-8)
})

test_that("rapid() is close to Oscar's RAPID for rf and gbm", {
  skip_on_cran()
  rf <- rapid(fx$num$o, fx$num$s, qi, "y", model_type = "rf", num_dist = "grid", seed = 1)
  expect_lt(abs(rf$rapid - fx$rapid$rf), 0.05)
  skip_if_not_installed("xgboost")
  gbm <- rapid(fx$num$o, fx$num$s, qi, "y", model_type = "gbm", num_dist = "grid", seed = 1)
  expect_lt(abs(gbm$rapid - fx$rapid$gbm), 0.05)
})

test_that("point rule reproduces Oscar's flags", {
  x <- rapid(fx$num$o, fx$num$s, qi, "y", model_type = "lm", num_eval_method = "point",
             num_scale = "identity")
  expect_identical(x$at_risk, fx$point_rule$at_risk)
  expect_equal(x$rapid, fx$point_rule$rapid, tolerance = 1e-12)
  expect_s3_class(x$rule, "rapid_rule_point")
  expect_identical(x$threshold, 0.1)
})

test_that("point-rule baseline equals Oscar's marginal_share", {
  x <- rapid(fx$num$o, fx$num$s, qi, "y", model_type = "lm", num_eval_method = "point",
             num_scale = "identity")
  expect_equal(x$baseline, fx$point_rule$marginal_share, tolerance = 1e-12)
})

test_that("categorical rapid() is close to Oscar's RCS_marginal result", {
  skip_on_cran()
  x <- rapid(fx$cat$o, fx$cat$s, qi, "y", model_type = "rf", seed = 1)
  expect_lt(abs(x$rapid - fx$cat$rapid), 0.05)
})

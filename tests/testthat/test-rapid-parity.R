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

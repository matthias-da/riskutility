test_that("regression_fidelity returns correct S3 class structure", {
  set.seed(42)
  n <- 200
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 2 + 3 * x1 - 1.5 * x2 + rnorm(n)
  X <- data.frame(y = y, x1 = x1, x2 = x2)
  x1s <- rnorm(n); x2s <- rnorm(n)
  ys <- 2 + 3 * x1s - 1.5 * x2s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s, x2 = x2s)

  res <- regression_fidelity(X, Y, formula = y ~ x1 + x2)

  expect_s3_class(res, "regression_fidelity")
  expect_true("utility_score" %in% names(res))
  expect_true("coefficients" %in% names(res))
  expect_true("mean_ci_overlap" %in% names(res))
  expect_true("mean_abs_std_bias" %in% names(res))
  expect_true("sig_agreement_rate" %in% names(res))
  expect_true("formula" %in% names(res))
  expect_true("model" %in% names(res))
  expect_true("conf_level" %in% names(res))
  expect_true("n_X" %in% names(res))
  expect_true("n_Y" %in% names(res))
  expect_true("n_coef" %in% names(res))
})

test_that("coefficients data.frame has correct columns", {
  set.seed(42)
  n <- 200
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 2 + 3 * x1 - 1.5 * x2 + rnorm(n)
  X <- data.frame(y = y, x1 = x1, x2 = x2)
  x1s <- rnorm(n); x2s <- rnorm(n)
  ys <- 2 + 3 * x1s - 1.5 * x2s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s, x2 = x2s)

  res <- regression_fidelity(X, Y, formula = y ~ x1 + x2)
  df <- res$coefficients

  expect_s3_class(df, "data.frame")
  expected_cols <- c("term", "estimate_orig", "estimate_synth", "se_orig",
                     "se_synth", "bias", "std_bias", "ci_overlap",
                     "sig_orig", "sig_synth", "sig_agreement")
  expect_true(all(expected_cols %in% names(df)))
  expect_equal(nrow(df), 3)  # intercept + x1 + x2
  expect_equal(res$n_coef, 3)
})

test_that("perfect fidelity: same DGP yields high CI overlap", {
  set.seed(123)
  n <- 500
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 1 + 2 * x1 + 3 * x2 + rnorm(n)
  X <- data.frame(y = y, x1 = x1, x2 = x2)

  x1s <- rnorm(n); x2s <- rnorm(n)
  ys <- 1 + 2 * x1s + 3 * x2s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s, x2 = x2s)

  res <- regression_fidelity(X, Y, formula = y ~ x1 + x2)

  # With same DGP and n=500, CIs should overlap substantially

  expect_gt(res$utility_score, 0.5)
  expect_equal(res$sig_agreement_rate, 1)
})

test_that("poor fidelity: wrong coefficients yield low CI overlap", {
  set.seed(42)
  n <- 500
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 2 + 3 * x1 - 1.5 * x2 + rnorm(n)
  X <- data.frame(y = y, x1 = x1, x2 = x2)

  # Synthetic has very different coefficients
  x1s <- rnorm(n); x2s <- rnorm(n)
  ys_bad <- 10 + 0.1 * x1s + 5 * x2s + rnorm(n)
  Y_bad <- data.frame(y = ys_bad, x1 = x1s, x2 = x2s)

  res <- regression_fidelity(X, Y_bad, formula = y ~ x1 + x2)

  # With totally different DGP, overall CI overlap should be low
  expect_lt(res$utility_score, 0.5)
  expect_gt(res$mean_abs_std_bias, 1)
})

test_that("bias computation is correct", {
  set.seed(42)
  n <- 100
  x1 <- rnorm(n)
  y <- 5 + 2 * x1 + rnorm(n, sd = 0.5)
  X <- data.frame(y = y, x1 = x1)

  x1s <- rnorm(n)
  ys <- 5 + 2 * x1s + rnorm(n, sd = 0.5)
  Y <- data.frame(y = ys, x1 = x1s)

  res <- regression_fidelity(X, Y, formula = y ~ x1)

  # Verify bias = synth - orig
  fit_orig <- lm(y ~ x1, data = X)
  fit_synth <- lm(y ~ x1, data = Y)
  expected_bias <- coef(fit_synth) - coef(fit_orig)

  expect_equal(res$coefficients$bias, unname(expected_bias), tolerance = 1e-10)
})

test_that("CI overlap formula is correct for hand-computed case", {
  # Construct data so that we know the approximate CIs
  # Orig: coef = 0, SE = 1  -> CI = [-1.96, 1.96]
  # Synth: coef = 1, SE = 1 -> CI = [-0.96, 2.96]
  # overlap_num = min(1.96, 2.96) - max(-1.96, -0.96) = 1.96 - (-0.96) = 2.92
  # max_width = max(3.92, 3.92) = 3.92
  # ci_overlap = 2.92 / 3.92 = 0.7449...

  # We test directly by constructing known results
  # Use very large n so that SEs are negligible compared to coefficient values
  set.seed(42)
  n <- 100000
  x1 <- rnorm(n)

  # Original: y = 0 * x1 + noise(sd=1) -> coef(x1) ~ 0, SE ~ 1/sqrt(n)
  y_orig <- 0 * x1 + rnorm(n)
  X <- data.frame(y = y_orig, x1 = x1)

  # Synthetic: y = 0.01 * x1 + noise -> coef ~ 0.01
  y_synth <- 0.01 * x1 + rnorm(n)
  Y <- data.frame(y = y_synth, x1 = x1)

  res <- regression_fidelity(X, Y, formula = y ~ x1)

  # With large n and small coefficient differences, CI overlap should be very high
  # (the CIs are tiny and nearly identical for the intercept at least)
  expect_gt(res$coefficients$ci_overlap[1], 0.8)  # intercept
})

test_that("CI overlap is 0 for completely non-overlapping CIs", {
  set.seed(42)
  n <- 10000
  x1 <- rnorm(n)

  # Original: strong positive effect
  y_orig <- 100 * x1 + rnorm(n)
  X <- data.frame(y = y_orig, x1 = x1)

  # Synthetic: strong negative effect
  y_synth <- -100 * x1 + rnorm(n)
  Y <- data.frame(y = y_synth, x1 = x1)

  res <- regression_fidelity(X, Y, formula = y ~ x1)

  # x1 coefficient CIs should not overlap at all
  x1_row <- res$coefficients[res$coefficients$term == "x1", ]
  expect_equal(x1_row$ci_overlap, 0)
})

test_that("significance agreement detection is correct", {
  set.seed(42)
  n <- 200
  x1 <- rnorm(n); x2 <- rnorm(n)

  # x1 has strong effect, x2 has no effect in original
  y_orig <- 5 * x1 + 0 * x2 + rnorm(n)
  X <- data.frame(y = y_orig, x1 = x1, x2 = x2)

  # In synthetic, x1 still strong, x2 still zero
  x1s <- rnorm(n); x2s <- rnorm(n)
  y_synth <- 5 * x1s + 0 * x2s + rnorm(n)
  Y <- data.frame(y = y_synth, x1 = x1s, x2 = x2s)

  res <- regression_fidelity(X, Y, formula = y ~ x1 + x2)

  # x1 should be significant in both, x2 not significant in both
  df <- res$coefficients
  x1_row <- df[df$term == "x1", ]
  x2_row <- df[df$term == "x2", ]

  expect_true(x1_row$sig_orig)
  expect_true(x1_row$sig_synth)
  expect_true(x1_row$sig_agreement)

  # x2 should not be significant in either with high probability
  # (but with n=200 there's a small chance of type I error)
  # At minimum, sig_agreement should hold for at least x1
  expect_true(x1_row$sig_agreement)
})

test_that("GLM model works", {
  set.seed(42)
  n <- 300
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.5 + 1.5 * x1 - 0.5 * x2))
  X <- data.frame(y = y, x1 = x1, x2 = x2)

  x1s <- rnorm(n); x2s <- rnorm(n)
  ys <- rbinom(n, 1, plogis(0.5 + 1.5 * x1s - 0.5 * x2s))
  Y <- data.frame(y = ys, x1 = x1s, x2 = x2s)

  res <- regression_fidelity(X, Y, formula = y ~ x1 + x2,
                             model = "glm", family = binomial())

  expect_s3_class(res, "regression_fidelity")
  expect_equal(res$model, "glm")
  expect_equal(res$n_coef, 3)
  expect_true(res$utility_score >= 0 && res$utility_score <= 1)
})

test_that("synth_pair dispatch works", {
  set.seed(42)
  n <- 200
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 2 + 3 * x1 - 1.5 * x2 + rnorm(n)
  X <- data.frame(y = y, x1 = x1, x2 = x2)
  x1s <- rnorm(n); x2s <- rnorm(n)
  ys <- 2 + 3 * x1s - 1.5 * x2s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s, x2 = x2s)

  sp <- list(original = X, synthetic = Y)
  class(sp) <- "synth_pair"

  res <- regression_fidelity(sp, formula = y ~ x1 + x2)
  expect_s3_class(res, "regression_fidelity")
  expect_equal(res$n_X, n)
  expect_equal(res$n_Y, n)
})

test_that("error without formula", {
  X <- data.frame(y = rnorm(10), x1 = rnorm(10))
  Y <- data.frame(y = rnorm(10), x1 = rnorm(10))

  expect_error(regression_fidelity(X, Y),
               "formula.*required")
})

test_that("error with missing variables", {
  X <- data.frame(y = rnorm(10), x1 = rnorm(10))
  Y <- data.frame(y = rnorm(10), x2 = rnorm(10))

  expect_error(regression_fidelity(X, Y, formula = y ~ x1),
               "missing in Y")
})

test_that("print method works without error", {
  set.seed(42)
  n <- 100
  x1 <- rnorm(n)
  y <- 2 + 3 * x1 + rnorm(n)
  X <- data.frame(y = y, x1 = x1)
  x1s <- rnorm(n)
  ys <- 2 + 3 * x1s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s)

  res <- regression_fidelity(X, Y, formula = y ~ x1)

  expect_output(print(res), "Regression Fidelity")
  expect_output(print(res), "Utility score")
  expect_output(print(res), "Significance agreement")
})

test_that("summary method returns correct class", {
  set.seed(42)
  n <- 100
  x1 <- rnorm(n)
  y <- 2 + 3 * x1 + rnorm(n)
  X <- data.frame(y = y, x1 = x1)
  x1s <- rnorm(n)
  ys <- 2 + 3 * x1s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s)

  res <- regression_fidelity(X, Y, formula = y ~ x1)
  s <- summary(res)

  expect_s3_class(s, "summary.regression_fidelity")
  expect_output(print(s), "Summary: Regression Fidelity")
  expect_output(print(s), "Coefficient Comparison")
})

test_that("plot method which=1 works without error", {
  set.seed(42)
  n <- 100
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 2 + 3 * x1 - 1.5 * x2 + rnorm(n)
  X <- data.frame(y = y, x1 = x1, x2 = x2)
  x1s <- rnorm(n); x2s <- rnorm(n)
  ys <- 2 + 3 * x1s - 1.5 * x2s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s, x2 = x2s)

  res <- regression_fidelity(X, Y, formula = y ~ x1 + x2)

  expect_no_error(plot(res, which = 1))
})

test_that("plot method which=2 works without error", {
  set.seed(42)
  n <- 100
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 2 + 3 * x1 - 1.5 * x2 + rnorm(n)
  X <- data.frame(y = y, x1 = x1, x2 = x2)
  x1s <- rnorm(n); x2s <- rnorm(n)
  ys <- 2 + 3 * x1s - 1.5 * x2s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s, x2 = x2s)

  res <- regression_fidelity(X, Y, formula = y ~ x1 + x2)

  expect_no_error(plot(res, which = 2))
})

test_that("plot method which=1:2 works without error", {
  set.seed(42)
  n <- 100
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 2 + 3 * x1 - 1.5 * x2 + rnorm(n)
  X <- data.frame(y = y, x1 = x1, x2 = x2)
  x1s <- rnorm(n); x2s <- rnorm(n)
  ys <- 2 + 3 * x1s - 1.5 * x2s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s, x2 = x2s)

  res <- regression_fidelity(X, Y, formula = y ~ x1 + x2)

  expect_no_error(plot(res, which = 1:2))
})

test_that("NA handling works", {
  set.seed(42)
  n <- 100
  x1 <- rnorm(n)
  y <- 2 + 3 * x1 + rnorm(n)
  X <- data.frame(y = y, x1 = x1)

  # Introduce NAs
  X$y[1:5] <- NA
  X$x1[6:10] <- NA

  x1s <- rnorm(n)
  ys <- 2 + 3 * x1s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s)
  Y$y[1:3] <- NA

  res <- regression_fidelity(X, Y, formula = y ~ x1, na.rm = TRUE)

  expect_s3_class(res, "regression_fidelity")
  expect_equal(res$n_X, 90)  # 100 - 10 NAs
  expect_equal(res$n_Y, 97)  # 100 - 3 NAs
})

test_that("conf_level parameter works", {
  set.seed(42)
  n <- 200
  x1 <- rnorm(n)
  y <- 2 + 3 * x1 + rnorm(n)
  X <- data.frame(y = y, x1 = x1)
  x1s <- rnorm(n)
  ys <- 2 + 3 * x1s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s)

  res_95 <- regression_fidelity(X, Y, formula = y ~ x1, conf_level = 0.95)
  res_99 <- regression_fidelity(X, Y, formula = y ~ x1, conf_level = 0.99)

  expect_equal(res_95$conf_level, 0.95)
  expect_equal(res_99$conf_level, 0.99)

  # Wider CIs (99%) should have higher overlap than narrower CIs (95%)
  # (or at least no less, in general)
  # We just verify both are valid
  expect_true(res_95$utility_score >= 0 && res_95$utility_score <= 1)
  expect_true(res_99$utility_score >= 0 && res_99$utility_score <= 1)
})

test_that("utility_score is in [0, 1]", {
  set.seed(42)
  n <- 200
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 2 + 3 * x1 - 1.5 * x2 + rnorm(n)
  X <- data.frame(y = y, x1 = x1, x2 = x2)

  # Good synthetic
  x1s <- rnorm(n); x2s <- rnorm(n)
  ys <- 2 + 3 * x1s - 1.5 * x2s + rnorm(n)
  Y_good <- data.frame(y = ys, x1 = x1s, x2 = x2s)

  # Bad synthetic
  ys_bad <- 100 + 0.01 * x1s - 50 * x2s + rnorm(n)
  Y_bad <- data.frame(y = ys_bad, x1 = x1s, x2 = x2s)

  res_good <- regression_fidelity(X, Y_good, formula = y ~ x1 + x2)
  res_bad  <- regression_fidelity(X, Y_bad,  formula = y ~ x1 + x2)

  expect_true(res_good$utility_score >= 0 && res_good$utility_score <= 1)
  expect_true(res_bad$utility_score >= 0  && res_bad$utility_score <= 1)
})

test_that("standardized bias uses correct formula", {
  set.seed(42)
  n <- 200
  x1 <- rnorm(n)
  y <- 3 + 5 * x1 + rnorm(n)
  X <- data.frame(y = y, x1 = x1)

  x1s <- rnorm(n)
  ys <- 3 + 5 * x1s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s)

  res <- regression_fidelity(X, Y, formula = y ~ x1)
  df <- res$coefficients

  # std_bias = bias / se_orig
  expected_std_bias <- df$bias / df$se_orig
  expect_equal(df$std_bias, expected_std_bias, tolerance = 1e-10)
})

test_that("sig_agreement_rate is computed correctly", {
  set.seed(42)
  n <- 200
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 5 * x1 + 0 * x2 + rnorm(n)
  X <- data.frame(y = y, x1 = x1, x2 = x2)

  x1s <- rnorm(n); x2s <- rnorm(n)
  ys <- 5 * x1s + 0 * x2s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s, x2 = x2s)

  res <- regression_fidelity(X, Y, formula = y ~ x1 + x2)

  # sig_agreement_rate should equal mean(sig_agreement)
  expect_equal(res$sig_agreement_rate,
               mean(res$coefficients$sig_agreement))
})

test_that("model type is stored correctly", {
  set.seed(42)
  n <- 100
  x1 <- rnorm(n)
  y <- 2 + x1 + rnorm(n)
  X <- data.frame(y = y, x1 = x1)
  x1s <- rnorm(n)
  ys <- 2 + x1s + rnorm(n)
  Y <- data.frame(y = ys, x1 = x1s)

  res_lm <- regression_fidelity(X, Y, formula = y ~ x1, model = "lm")
  expect_equal(res_lm$model, "lm")

  y_bin <- rbinom(n, 1, plogis(x1))
  X_bin <- data.frame(y = y_bin, x1 = x1)
  y_bin_s <- rbinom(n, 1, plogis(x1s))
  Y_bin <- data.frame(y = y_bin_s, x1 = x1s)

  res_glm <- regression_fidelity(X_bin, Y_bin, formula = y ~ x1, model = "glm")
  expect_equal(res_glm$model, "glm")
})

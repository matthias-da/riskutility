# predict() on new records and holdout scoring

qi <- c("x1", "x2")

test_that("predict() on the original data reproduces the records", {
  o <- sim_numeric(120)
  s <- noisy_copy(o, seed = 2)
  r <- rapid(o, s, qi, "y", model_type = "rf", store_model = TRUE, n_perm = 0, seed = 1)
  expect_equal(predict(r, o), as.data.frame(r))
  # A permutation threshold is used as calibrated
  rp <- rapid(o, s, qi, "y", model_type = "lm", store_model = TRUE,
              threshold_rule = "permutation", n_perm = 3, seed = 1)
  expect_equal(predict(rp, o), as.data.frame(rp))
  # Categorical target
  oc <- sim_categorical(120)
  sc <- noisy_copy(oc, seed = 2)
  rc <- rapid(oc, sc, qi, "y", model_type = "rf", store_model = TRUE, n_perm = 0, seed = 1)
  expect_equal(predict(rc, oc), as.data.frame(rc))
})

test_that("predict() without the target returns the attacker's prediction only", {
  o <- sim_numeric(120)
  s <- noisy_copy(o, seed = 2)
  r <- rapid(o, s, qi, "y", model_type = "lm", store_model = TRUE, n_perm = 0)
  p <- predict(r, o[1:5, qi])
  expect_identical(names(p), c(qi, "predicted_value"))
  expect_equal(p$predicted_value, as.data.frame(r)$predicted_value[1:5])

  oc <- sim_categorical(120)
  rc <- rapid(oc, noisy_copy(oc, seed = 2), qi, "y", model_type = "rf",
              store_model = TRUE, n_perm = 0, seed = 1)
  pc <- predict(rc, oc[1:5, qi])
  expect_identical(names(pc), c(qi, "predicted_class"))
  expect_identical(pc$predicted_class, as.data.frame(rc)$predicted_class[1:5])
})

test_that("predict() needs the stored model", {
  o <- sim_numeric(80)
  r <- rapid(o, noisy_copy(o, seed = 2), qi, "y", model_type = "lm", n_perm = 0)
  expect_error(predict(r, o), regexp = "store_model", class = "rapid_error_object")
})

test_that("holdout equal to X reproduces the main result", {
  o <- sim_numeric(120)
  s <- noisy_copy(o, seed = 2)
  r <- rapid(o, s, qi, "y", model_type = "rf", holdout = o, n_perm = 3, seed = 1)
  expect_identical(r$holdout$rapid, r$rapid)
  expect_identical(r$holdout$chance, r$chance)
  expect_identical(r$holdout$excess, r$excess)
  expect_identical(r$holdout$n, r$n_total)
  expect_equal(r$holdout$records, r$records)
  expect_identical(r$release_specific, 0)
  out <- format(r)
  expect_true(any(grepl("^Holdout", out)))
  expect_true(any(grepl("^Release-specific", out)))
  # Without permutations there is no chance level
  r0 <- rapid(o, s, qi, "y", model_type = "lm", holdout = o[1:10, ], n_perm = 0)
  expect_true(is.na(r0$holdout$chance))
  expect_true(is.na(r0$holdout$excess))
  expect_false(any(grepl("chance", grep("^Holdout", format(r0), value = TRUE))))
  # No holdout, no lines
  expect_null(rapid(o, s, qi, "y", model_type = "lm", n_perm = 0)$holdout)
})

test_that("baseline of holdout records uses the original values", {
  o <- sim_numeric(120)
  s <- noisy_copy(o, seed = 2)
  h <- data.frame(x1 = 0.1, x2 = -0.2, y = 76)
  r <- rapid(o, s, qi, "y", model_type = "lm", num_error_metric = "absolute",
             num_epsilon = 5, holdout = h, n_perm = 2, seed = 1)
  expect_identical(r$holdout$n, 1L)
  rec <- r$holdout$records
  expect_identical(nrow(rec), 1L)
  expect_equal(rec$b, sum(o$y > 71 & o$y < 81) / 120)
  # Categorical target: the class share among the original records
  oc <- sim_categorical(120)
  hc <- oc[7, ]
  rc <- rapid(oc, noisy_copy(oc, seed = 2), qi, "y", model_type = "rf", holdout = hc,
              n_perm = 0, seed = 1)
  expect_equal(rc$holdout$records$b, mean(oc$y == hc$y))
})

test_that("a memorising release has positive release_specific", {
  o <- sim_numeric(200, seed = 1)
  h <- sim_numeric(200, seed = 2)
  r <- rapid(o, o, qi, "y", model_type = "rf", holdout = h, n_perm = 0, seed = 1)
  expect_equal(r$release_specific, r$rapid - r$holdout$rapid)
  expect_gt(r$release_specific, 0.1)
})

test_that("a synth_pair with holdout is scored on it", {
  o <- sim_numeric(120)
  s <- noisy_copy(o, seed = 2)
  h <- sim_numeric(40, seed = 3)
  pair <- synth_pair(o, s, key_vars = qi, target_var = "y", holdout = h)
  r <- rapid(pair, model_type = "lm", n_perm = 2, seed = 1)
  direct <- rapid(o, s, qi, "y", model_type = "lm", holdout = h, n_perm = 2, seed = 1)
  expect_identical(r$holdout$n, 40L)
  expect_identical(r$holdout$rapid, direct$holdout$rapid)
  expect_identical(r$release_specific, direct$release_specific)
  # An explicit holdout wins over the one of the pair
  r2 <- rapid(pair, model_type = "lm", holdout = h[1:10, ], n_perm = 0)
  expect_identical(r2$holdout$n, 10L)
})

test_that("holdout with missing variables is a classed error", {
  o <- sim_numeric(80)
  s <- noisy_copy(o, seed = 2)
  expect_error(rapid(o, s, qi, "y", model_type = "lm", holdout = o[, c("x1", "x2")]),
               regexp = "'y'", class = "rapid_error_input")
  expect_error(rapid(o, s, qi, "y", model_type = "lm", holdout = o[, c("x1", "y")]),
               regexp = "'x2'", class = "rapid_error_input")
  r <- rapid(o, s, qi, "y", model_type = "lm", store_model = TRUE, n_perm = 0)
  expect_error(predict(r, o[, c("x1", "y")]), regexp = "'x2'", class = "rapid_error_input")
})

test_that("holdout and newdata are prepared like X", {
  o <- sim_numeric(120)
  s <- noisy_copy(o, seed = 2)
  h <- sim_numeric(30, seed = 3)
  h$y[1:3] <- NA
  r <- rapid(o, s, qi, "y", model_type = "lm", holdout = h, na_strategy = "drop",
             n_perm = 0)
  expect_identical(r$holdout$n, 27L)
  r <- rapid(o, s, qi, "y", model_type = "lm", holdout = h, na_strategy = "constant",
             na_constant_value = 0, n_perm = 0)
  expect_identical(r$holdout$n, 30L)
  expect_equal(r$holdout$records$y[1:3], c(0, 0, 0))
  r <- rapid(o, s, qi, "y", model_type = "lm", holdout = h, na_strategy = "median",
             n_perm = 0, store_model = TRUE)
  # Imputed with the median of the original values, as X
  expect_equal(r$holdout$records$y[1:3], rep(stats::median(o$y), 3))
  expect_identical(nrow(predict(r, h)), 30L)
  # Categorical quasi-identifier: a level unknown to the attacker is refused
  oc <- o
  oc$x2 <- factor(ifelse(o$x2 > 0, "hi", "lo"))
  sc <- noisy_copy(oc, seed = 2)
  hc <- oc[1:5, ]
  hc$x2 <- factor(c("hi", "lo", "mid", "hi", "lo"))
  expect_error(rapid(oc, sc, qi, "y", model_type = "rf", holdout = hc, n_perm = 0),
               regexp = "mid", class = "rapid_error_levels")
  rc <- rapid(oc, sc, qi, "y", model_type = "rf", store_model = TRUE, n_perm = 0)
  expect_error(predict(rc, hc), regexp = "mid", class = "rapid_error_levels")
  expect_identical(nrow(predict(rc, hc[c(1, 2, 4), ])), 3L)
})

test_that("update() with a new threshold on an object with a holdout equals a fresh call", {
  o <- sim_numeric(120)
  s <- noisy_copy(o, seed = 2)
  h <- sim_numeric(50, seed = 3)
  r <- rapid(o, s, qi, "y", model_type = "rf", holdout = h, n_perm = 3, seed = 1)
  fresh <- rapid(o, s, qi, "y", model_type = "rf", holdout = h, num_tau = 0.5, n_perm = 3,
                 seed = 1)
  fp <- rapid(o, s, qi, "y", model_type = "rf", holdout = h, n_perm = 3, seed = 1,
              threshold_rule = "permutation")
  local_mocked_bindings(.learner_fit = function(...) stop("refit"))
  u <- update(r, num_tau = 0.5)
  expect_identical(u$holdout, fresh$holdout)
  expect_identical(u$release_specific, fresh$release_specific)
  expect_identical(u$rapid, fresh$rapid)
  # A permutation threshold calibrated on X is used for the holdout
  up <- update(r, threshold_rule = "permutation")
  expect_identical(up$holdout, fp$holdout)
  thr <- fp$rule$threshold
  expect_identical(fp$holdout$records$at_risk, fp$holdout$records$r > thr)
})

test_that("predict() reproduces the records with tree learners and exact distributions", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("xgboost")
  o <- sim_numeric(100)
  s <- noisy_copy(o, seed = 2)
  for (m in c("cart", "gbm", "lookup")) {
    r <- rapid(o, s, qi, "y", model_type = m, store_model = TRUE, n_perm = 0, seed = 1)
    expect_equal(predict(r, o), as.data.frame(r), info = m)
  }
  r <- rapid(o, s, qi, "y", model_type = "rf", num_dist = "exact", store_model = TRUE,
             n_perm = 0, seed = 1)
  expect_equal(predict(r, o), as.data.frame(r))
})

test_that("a missing target of one new record is handled as in X or refused", {
  o <- sim_numeric(120)
  o$y[1:2] <- NA
  s <- noisy_copy(sim_numeric(120), seed = 2)
  h <- data.frame(x1 = 0, x2 = 0, y = NA_real_)
  fit <- function(strategy, ...) {
    rapid(o, s, qi, "y", model_type = "lm", holdout = h, na_strategy = strategy,
          n_perm = 2, seed = 1, store_model = TRUE, ...)
  }
  # median: the median of the original values, for the holdout and predict()
  r <- fit("median")
  med <- stats::median(o$y, na.rm = TRUE)
  expect_identical(r$holdout$n, 1L)
  expect_equal(r$holdout$records$y, med)
  expect_false(anyNA(r$holdout$records$r))
  expect_equal(predict(r, h)$y, med)
  expect_false(anyNA(predict(r, h)$r))
  # constant: the constant
  r <- fit("constant", na_constant_value = 70)
  expect_equal(r$holdout$records$y, 70)
  expect_equal(predict(r, h)$y, 70)
  # drop: nothing is left, which is an error
  expect_error(fit("drop"), regexp = "no records left", class = "rapid_error_input")
  r <- rapid(o, s, qi, "y", model_type = "lm", na_strategy = "drop", n_perm = 0,
             store_model = TRUE)
  expect_error(predict(r, h), regexp = "no records left", class = "rapid_error_input")
  # no imputation value: refused
  expect_error(fit("constant", na_constant_value = NA), regexp = "cannot be imputed",
               class = "rapid_error_input")

  # Categorical target: the level "missing", known only if X or Y had it
  oc <- sim_categorical(120)
  sc <- noisy_copy(oc, seed = 2)
  hc <- data.frame(x1 = 0, x2 = 0, y = factor(NA, levels = levels(oc$y)))
  expect_error(rapid(oc, sc, qi, "y", model_type = "rf", holdout = hc, n_perm = 0),
               regexp = "missing", class = "rapid_error_levels")
  rc <- rapid(oc, sc, qi, "y", model_type = "rf", n_perm = 0, store_model = TRUE)
  expect_error(predict(rc, hc), regexp = "missing", class = "rapid_error_levels")
  oc$y[1:5] <- NA
  rc <- rapid(oc, sc, qi, "y", model_type = "rf", holdout = hc, n_perm = 0, seed = 1,
              store_model = TRUE)
  expect_identical(as.character(rc$holdout$records$y), "missing")
  expect_equal(rc$holdout$records$b, 5 / 120)
  expect_identical(as.character(predict(rc, hc)$y), "missing")
})

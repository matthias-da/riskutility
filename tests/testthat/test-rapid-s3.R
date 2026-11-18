# The rapid object: constructor, validator and basic methods

qi <- c("x1", "x2")

lm_object <- function(...) {
  o <- sim_numeric(100)
  rapid(o, noisy_copy(o), qi, "y", model_type = "lm", ...)
}

test_that("new_rapid() sets class and version", {
  x <- lm_object()
  expect_s3_class(x, "rapid", exact = TRUE)
  expect_identical(x$version, utils::packageVersion("riskutility"))
  expect_true(all(c("rapid", "n_at_risk", "pct_at_risk", "n_total", "rule", "method",
                    "threshold", "score", "baseline", "at_risk", "records", "attack",
                    "model_type", "attacker", "model_metrics", "formula", "key_vars",
                    "target_var", "is_categorical", "na", "call", "version") %in% names(x)))
  expect_identical(validate_rapid(x), x)
})

test_that("validate_rapid() rejects inconsistent objects", {
  x <- lm_object()
  bad <- x
  bad$score <- bad$score[-1]
  expect_error(validate_rapid(bad), class = "rapid_error_object")
  bad <- x
  bad$n_total <- 99L
  expect_error(validate_rapid(bad), class = "rapid_error_object")
  bad <- x
  bad$at_risk <- NULL
  expect_error(validate_rapid(bad), regexp = "at_risk", class = "rapid_error_object")
  bad <- x
  bad$records <- bad$records[-1, ]
  expect_error(validate_rapid(bad), class = "rapid_error_object")
})

test_that("a 0.2.0 object is refused with advice", {
  old <- structure(
    list(rapid = 0.4, n_at_risk = 4, pct_at_risk = 40, n_total = 10, method = "RCS_conditional",
         threshold = 1, model_type = "rf", records = data.frame(at_risk = rep(TRUE, 4)),
         is_categorical = TRUE),
    class = "rapid"
  )
  calls <- list(print = function() print(old), summary = function() summary(old),
                confint = function() confint(old),
                as.data.frame = function() as.data.frame(old),
                format = function() format(old), nobs = function() nobs(old),
                plot = function() plot(old))
  for (nm in names(calls)) {
    err <- tryCatch(calls[[nm]](), error = function(e) e)
    expect_s3_class(err, "rapid_error_object")
    expect_match(conditionMessage(err), "re-run rapid()", fixed = TRUE, info = nm)
  }
})

test_that("methods that need stored parts say which argument to set", {
  x <- lm_object(store_predictions = FALSE)
  expect_error(.rapid_require(x, "predictions"), regexp = "store_predictions",
               class = "rapid_error_object")
  expect_error(.rapid_require(x, "model"), regexp = "store_model",
               class = "rapid_error_object")
  expect_silent(.rapid_require(x, "version"))
  expect_silent(.rapid_require(lm_object(store_model = TRUE), c("version", "predictions",
                                                                 "model")))
})

test_that("print returns its input invisibly and equals format()", {
  x <- lm_object()
  out <- utils::capture.output(res <- withVisible(print(x)))
  expect_false(res$visible)
  expect_identical(res$value, x)
  expect_identical(out, format(x))
  expect_type(format(x), "character")
  expect_match(format(x)[2], "^Rule: band: symmetric error < 10 %, gain > 0.3$")
  expect_match(format(x)[length(format(x))],
               "^Next: summary\\(\\), plot\\(\\), confint\\(\\), rapid_test\\(\\)$")
})

test_that("print output is stable", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(150)
  s <- noisy_copy(o)
  expect_snapshot(print(rapid(o, s, qi, "y", model_type = "lm", seed = 1)))
  expect_snapshot(print(rapid(o, s, qi, "y", model_type = "lm", num_eval_method = "point", seed = 1)))
  expect_snapshot(print(rapid(to_classes(o), to_classes(s, ref = o), qi, "y",
                              model_type = "cart", seed = 1)))
})

test_that("as.data.frame returns all records, or the flagged ones", {
  x <- lm_object()
  all_rec <- as.data.frame(x)
  expect_s3_class(all_rec, "data.frame")
  expect_identical(all_rec, x$records)
  flagged <- as.data.frame(x, at_risk_only = TRUE)
  expect_equal(nrow(flagged), x$n_at_risk)
  expect_true(all(flagged$at_risk))
  expect_true(x$n_at_risk > 0 && x$n_at_risk < x$n_total)
})

test_that("nobs equals n_total", {
  x <- lm_object()
  expect_identical(stats::nobs(x), x$n_total)
  expect_equal(stats::nobs(x), 100)
})

test_that("the stored call can be re-evaluated", {
  o <- sim_numeric(100)
  s <- noisy_copy(o)
  x <- rapid(o, s, qi, "y", model_type = "rf", seed = 3)
  expect_identical(x$call[[1]], quote(rapid))
  expect_identical(eval(x$call)$score, x$score)

  alias <- rapid(original_data = o, synthetic_data = s, quasi_identifiers = qi,
                 sensitive_attribute = "y", model_type = "rf", seed = 3)
  expect_identical(alias$call[[1]], quote(rapid))
  expect_identical(eval(alias$call)$score, x$score)

  sp <- synth_pair(o, s, key_vars = qi, target_var = "y")
  paired <- rapid(sp, model_type = "rf", seed = 3)
  expect_identical(paired$call[[1]], quote(rapid))
  expect_identical(eval(paired$call)$score, x$score)
})

test_that("the object does not carry the calling environment", {
  x <- lm_object()
  expect_identical(environment(x$formula), globalenv())
})

# ---- update() ----------------------------------------------------------------------

test_that("update() with a new threshold does not refit and equals a fresh call", {
  o <- sim_numeric(150)
  s <- noisy_copy(o, seed = 2)
  r <- rapid(o, s, qi, "y", model_type = "lm", seed = 1)
  fresh <- rapid(o, s, qi, "y", model_type = "lm", num_tau = 0.5, seed = 1)
  local_mocked_bindings(.learner_fit = function(...) stop("refit"))
  u <- update(r, num_tau = 0.5)
  expect_identical(u$rapid, fresh$rapid)
  expect_identical(u$at_risk, fresh$at_risk)
  expect_identical(u$rule$threshold, fresh$rule$threshold)
  expect_identical(u$chance_rates, fresh$chance_rates)
  expect_identical(u$records, fresh$records)
  expect_identical(u$call$num_tau, 0.5)
  expect_false(identical(u$rapid, r$rapid))
})

test_that("update() switches between band and point without refit", {
  o <- sim_numeric(150)
  s <- noisy_copy(o, seed = 2)
  r <- rapid(o, s, qi, "y", model_type = "lm", seed = 1)
  fresh <- rapid(o, s, qi, "y", model_type = "lm", num_eval_method = "point", seed = 1)
  local_mocked_bindings(.learner_fit = function(...) stop("refit"))
  p <- update(r, num_eval_method = "point")
  expect_s3_class(p$rule, "rapid_rule_point")
  expect_identical(p$rapid, fresh$rapid)
  expect_identical(p$at_risk, fresh$at_risk)
  back <- update(p, num_eval_method = "band")
  expect_s3_class(back$rule, "rapid_rule_band")
  expect_identical(back$at_risk, r$at_risk)
  expect_identical(back$rapid, r$rapid)
})

test_that("update() with a model argument refits", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(150)
  s <- noisy_copy(o, seed = 2)
  r <- rapid(o, s, qi, "y", model_type = "lm", seed = 1)
  expect_identical(update(r, model_type = "cart")$model_type, "cart")
  # A rescored object keeps the new scoring arguments in its call
  u <- update(r, num_tau = 0.5)
  refit <- update(u, model_type = "cart")
  expect_identical(refit$model_type, "cart")
  expect_identical(refit$rule$threshold, 0.5)
  # Without stored predictions every update refits
  np <- rapid(o, s, qi, "y", model_type = "lm", seed = 1, store_predictions = FALSE)
  expect_identical(update(np, num_tau = 0.5)$rule$threshold, 0.5)
  expect_identical(update(r, num_tau = 0.5, evaluate = FALSE)$num_tau, 0.5)
})

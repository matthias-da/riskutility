# rapid_suite(): RAPID over a suite of attacker models

qi <- c("x1", "x2")
o_suite <- sim_numeric(100)
s_suite <- noisy_copy(o_suite, seed = 2)
suite_args <- list(o_suite, s_suite, key_vars = qi, target_var = "y", n_perm = 3, seed = 1)

suite_models <- function() {
  c("lm", "lookup", if (requireNamespace("rpart", quietly = TRUE)) "cart")
}
suite_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      fit <<- do.call(rapid_suite, c(suite_args, list(model_types = suite_models())))
    }
    fit
  }
})
single <- function(m) do.call(rapid, c(suite_args, list(model_type = m)))

test_that("max and mean equal the single-model results", {
  s <- suite_fit()
  singles <- lapply(suite_models(), single)
  rates <- vapply(singles, function(f) f$rapid, numeric(1))
  expect_equal(s$results$rapid, rates)
  expect_equal(s$rapid_mean, mean(rates))
  expect_equal(s$rapid_max, max(rates))
  expect_equal(s$results$chance, vapply(singles, function(f) f$chance, numeric(1)))
})

test_that("worst_model names the maximum", {
  s <- suite_fit()
  expect_identical(s$worst_model, s$results$model_type[which.max(s$results$rapid)])
  expect_equal(s$fits[[s$worst_model]]$rapid, s$rapid_max)
})

test_that("excess_max and worst_model_excess follow the excess", {
  s <- suite_fit()
  expect_equal(s$excess_max, max(s$results$excess))
  expect_identical(s$worst_model_excess, s$results$model_type[which.max(s$results$excess)])
  expect_equal(s$results$excess, s$results$rapid - s$results$chance)
  # Without a chance level there is no excess
  s0 <- rapid_suite(o_suite, s_suite, key_vars = qi, target_var = "y", model_types = "lm",
                    n_perm = 0, seed = 1)
  expect_true(is.na(s0$excess_max))
  expect_true(is.na(s0$worst_model_excess))
  expect_false(any(grepl("excess over", format(s0))))
})

test_that("the default suite contains the lookup attacker", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("rpart")
  skip_if_not_installed("xgboost")
  s <- rapid_suite(o_suite, s_suite, key_vars = qi, target_var = "y", n_perm = 2, seed = 1)
  expect_identical(names(s$fits), c("rf", "cart", "gbm", "lookup"))
  expect_true(any(grepl("read the excess over the chance level", format(s))))
})

test_that("lm is skipped for a factor target with a classed warning", {
  oc <- to_classes(o_suite, probs = 0.5)
  sc <- to_classes(s_suite, ref = o_suite, probs = 0.5)
  expect_warning(
    s <- rapid_suite(oc, sc, key_vars = qi, target_var = "y",
                     model_types = c("lm", "logit"), n_perm = 2, seed = 1),
    class = "rapid_warning_skipped_model", regexp = "lm")
  expect_identical(names(s$fits), "logit")
  expect_error(
    rapid_suite(oc, sc, key_vars = qi, target_var = "y", model_types = "lm"),
    class = "rapid_error_model")
})

test_that("a missing backend is skipped with a warning", {
  local_mocked_bindings(.rapid_has_package = function(pkg) FALSE)
  expect_warning(
    s <- rapid_suite(o_suite, s_suite, key_vars = qi, target_var = "y",
                     model_types = c("lm", "rf"), n_perm = 2),
    class = "rapid_warning_skipped_model", regexp = "rf")
  expect_identical(names(s$fits), "lm")
})

test_that("other errors of rapid() pass through", {
  expect_error(rapid_suite(o_suite, s_suite, key_vars = qi, target_var = "y",
                           model_types = "lm", num_dist = "nope"), class = "rapid_error_input")
  expect_error(rapid_suite(o_suite, s_suite, key_vars = qi, target_var = "y",
                           model_types = "nope"), class = "rapid_error_input")
})

test_that("n_flagged_any equals the union of at_risk", {
  s <- suite_fit()
  union <- Reduce(`|`, lapply(s$fits, function(f) f$at_risk))
  expect_identical(s$n_flagged_any, sum(union))
  expect_true(s$n_flagged_any >= max(s$results$n_at_risk))
  expect_equal(summary(s)$n_flagged_all, sum(Reduce(`&`, lapply(s$fits, `[[`, "at_risk"))))
})

test_that("fitted objects combine to the same suite", {
  s <- suite_fit()
  parts <- unname(s$fits)
  s2 <- do.call(rapid_suite, parts)
  expect_identical(s2$results, s$results)
  expect_identical(s2$n_flagged_any, s$n_flagged_any)
  # Names given in the call are used; the same name twice is refused
  s3 <- rapid_suite(a = parts[[1]], b = parts[[1]])
  expect_identical(names(s3$fits), c("a", "b"))
  expect_error(rapid_suite(parts[[1]], parts[[1]]), class = "rapid_error_object")
  expect_error(rapid_suite(parts[[1]], "x"), class = "rapid_error_object")
})

test_that("objects from different targets are refused", {
  other <- rapid(o_suite, s_suite, "x1", "x2", model_type = "lm", n_perm = 3, seed = 1)
  expect_error(rapid_suite(a = single("lm"), b = other), class = "rapid_error_object",
               regexp = "target_var")
})

test_that("objects with a different rule or size are refused", {
  fit <- single("lm")
  expect_error(rapid_suite(a = fit, b = update(fit, num_tau = 0.5)), class = "rapid_error_object",
               regexp = "rule")
  small <- rapid(o_suite[1:50, ], s_suite[1:50, ], qi, "y", model_type = "lm", n_perm = 3)
  expect_error(rapid_suite(a = fit, b = small), class = "rapid_error_object", regexp = "n_total")
})

test_that("suite[['rf']] is a rapid object", {
  s <- suite_fit()
  expect_s3_class(s[["lm"]], "rapid")
  expect_identical(s[["lm"]], s$fits[["lm"]])
  expect_equal(s[[1]]$rapid, s$results$rapid[1])
})

test_that("the suite accepts a synth_pair and the old argument names", {
  sp <- synth_pair(o_suite, s_suite, key_vars = qi, target_var = "y")
  a <- rapid_suite(sp, model_types = "lm", n_perm = 3, seed = 1)
  b <- rapid_suite(original_data = o_suite, synthetic_data = s_suite,
                   quasi_identifiers = qi, sensitive_attribute = "y",
                   model_types = "lm", n_perm = 3, seed = 1)
  expect_equal(a$results, b$results)
})

test_that("as.data.frame returns the results", {
  s <- suite_fit()
  expect_identical(as.data.frame(s), s$results)
  expect_named(s$results, c("model_type", "rapid", "chance", "excess", "n_at_risk"))
})

test_that("a suite keeps the holdout of a synth_pair", {
  h <- sim_numeric(40, seed = 3)
  pair <- synth_pair(o_suite, s_suite, key_vars = qi, target_var = "y", holdout = h)
  su <- rapid_suite(pair, model_types = "lm", n_perm = 2, seed = 1)
  one <- rapid(pair, model_type = "lm", n_perm = 2, seed = 1)
  expect_false(is.null(su$fits[["lm"]]$holdout))
  expect_identical(su$fits[["lm"]]$holdout, one$holdout)
  expect_identical(su$fits[["lm"]]$release_specific, one$release_specific)
})

test_that("the skip check of a suite uses num_dist given in ...", {
  skip_if_not_installed("xgboost")
  w <- expect_warning(su <- rapid_suite(o_suite, s_suite, key_vars = qi, target_var = "y",
                                        model_types = c("lm", "gbm"), num_dist = "exact",
                                        n_perm = 0),
                      class = "rapid_warning_skipped_model")
  expect_match(conditionMessage(w), "gbm", fixed = TRUE)
  expect_identical(names(su$fits), "lm")
})

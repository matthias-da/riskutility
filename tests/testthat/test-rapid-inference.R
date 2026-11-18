# Tests for RAPID inferential framework:
#   confint.rapid(), rapid_test(), rapid_threshold_select(), rapid_synthesizer_cv()

# ---- Shared test data ----

.make_test_data <- function(n = 200, seed = 42) {
  set.seed(seed)
  original <- data.frame(
    age    = sample(20:60, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    income = round(rnorm(n, 50000, 10000))
  )
  synthetic <- data.frame(
    age    = original$age + sample(-3:3, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    income = round(original$income * runif(n, 0.9, 1.1))
  )
  list(original = original, synthetic = synthetic)
}

.make_test_data_cat <- function(n = 200, seed = 42) {
  set.seed(seed)
  original <- data.frame(
    age    = sample(20:60, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    health = factor(sample(c("Good", "Fair", "Poor"), n,
                           replace = TRUE, prob = c(0.6, 0.3, 0.1)))
  )
  synthetic <- data.frame(
    age    = original$age + sample(-3:3, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    health = factor(sample(c("Good", "Fair", "Poor"), n,
                           replace = TRUE, prob = c(0.55, 0.35, 0.1)))
  )
  list(original = original, synthetic = synthetic)
}

# ==========================================================================
# confint.rapid
# ==========================================================================

test_that("confint.rapid() returns correct structure for all methods", {
  d <- .make_test_data()
  r <- rapid(X = d$original, Y = d$synthetic,
             key_vars = c("age", "gender"), target_var = "income",
             model_type = "lm", num_epsilon = 0.1, verbose = FALSE)

  for (method in c("wald", "wilson", "bootstrap")) {
    ci <- confint(r, method = method, n_bootstrap = 200, seed = 1)
    expect_true(is.matrix(ci))
    expect_equal(nrow(ci), 1)
    expect_equal(ncol(ci), 2)
    expect_true(ci[1, 1] <= ci[1, 2])
    # Interval should contain the point estimate (or be close)
    expect_true(ci[1, 1] <= r$rapid + 0.01)
    expect_true(ci[1, 2] >= r$rapid - 0.01)
  }
})

test_that("bootstrap interval is a proper interval", {
  o <- sim_numeric(150)
  r <- rapid(o, noisy_copy(o), key_vars = c("x1", "x2"), target_var = "y",
             model_type = "lm")
  expect_true(r$rapid > 0 && r$rapid < 1)
  ci <- confint(r, method = "bootstrap", n_bootstrap = 2000, seed = 1)
  expect_true(ci[1, 1] > 0 && ci[1, 2] < 1)
  expect_true(ci[1, 1] < r$rapid && r$rapid < ci[1, 2])
  # Close to the Wilson interval
  expect_equal(unname(ci[1, ]), unname(confint(r)[1, ]), tolerance = 0.1)
})

test_that("confint.rapid() Wilson interval stays in [0, 1]", {
  d <- .make_test_data()
  r <- rapid(X = d$original, Y = d$synthetic,
             key_vars = c("age", "gender"), target_var = "income",
             model_type = "lm", num_epsilon = 0.1, verbose = FALSE)

  ci <- confint(r, method = "wilson")
  expect_true(ci[1, 1] >= 0)
  expect_true(ci[1, 2] <= 1)
})

test_that("confint.rapid() respects level parameter", {
  d <- .make_test_data()
  r <- rapid(X = d$original, Y = d$synthetic,
             key_vars = c("age", "gender"), target_var = "income",
             model_type = "lm", num_epsilon = 0.1, verbose = FALSE)

  ci90 <- confint(r, method = "wilson", level = 0.90)
  ci99 <- confint(r, method = "wilson", level = 0.99)
  # Wider confidence level -> wider interval

  width90 <- ci90[1, 2] - ci90[1, 1]
  width99 <- ci99[1, 2] - ci99[1, 1]
  expect_true(width99 > width90)
})

test_that("confint.rapid() works with categorical rapid object", {
  skip_if_not_installed("ranger")
  d <- .make_test_data_cat()
  r <- rapid(X = d$original, Y = d$synthetic,
             key_vars = c("age", "gender"), target_var = "health",
             model_type = "rf", cat_tau = 0.3,
             cat_eval_method = "RCS_marginal", verbose = FALSE)

  ci <- confint(r, method = "wilson")
  expect_true(is.matrix(ci))
  expect_true(ci[1, 1] >= 0 && ci[1, 2] <= 1)
})

# ==========================================================================
# rapid_test: a method on one fit
# ==========================================================================

qi <- c("x1", "x2")

# Linear attacker on data where the quasi-identifiers carry information about y
lm_fit <- function(n = 150, ..., shuffle_target = FALSE) {
  o <- sim_numeric(n)
  rapid(o, noisy_copy(o, seed = 2, shuffle_target = shuffle_target), qi, "y",
        model_type = "lm", ...)
}

test_that("rapid_test() on an object does not refit", {
  x <- lm_fit()
  local_mocked_bindings(.learner_fit = function(...) stop("refit"))
  t <- rapid_test(x, n_permutations = 19, seed = 1)
  expect_s3_class(t, "rapid_test")
  expect_length(t$null_distribution, 19)
  expect_identical(t$observed_rapid, x)
})

test_that("statistic equals the object's RAPID", {
  x <- lm_fit()
  t <- rapid_test(x, n_permutations = 19, seed = 1)
  expect_identical(t$statistic, x$rapid)
  expect_identical(t$rapid, x$rapid)
  expect_identical(t$chance, t$null_mean)
  expect_identical(t$chance, mean(t$null_distribution))
  expect_identical(t$excess, t$rapid - t$chance)
})

test_that("data interface and object interface agree", {
  o <- sim_numeric(150)
  s <- noisy_copy(o, seed = 2)
  by_data <- rapid_test(o, s, qi, "y", model_type = "rf", n_permutations = 19, seed = 5)
  by_object <- rapid_test(rapid(o, s, qi, "y", model_type = "rf", seed = 5, n_perm = 0),
                          n_permutations = 19, seed = 5)
  expect_identical(by_data$null_distribution, by_object$null_distribution)
  expect_identical(by_data$statistic, by_object$statistic)
  expect_identical(by_data$p_value, by_object$p_value)

  by_alias <- rapid_test(original_data = o, synthetic_data = s, quasi_identifiers = qi,
                         sensitive_attribute = "y", model_type = "rf",
                         n_permutations = 19, seed = 5)
  expect_identical(by_alias$null_distribution, by_data$null_distribution)
  by_pair <- rapid_test(synth_pair(o, s, key_vars = qi, target_var = "y"),
                        model_type = "rf", n_permutations = 19, seed = 5)
  expect_identical(by_pair$null_distribution, by_data$null_distribution)
  # The tested object keeps a call that can be evaluated again
  expect_identical(eval(by_data$observed_rapid$call)$score, by_data$observed_rapid$score)
})

test_that("p-value is (1 + #{null >= obs}) / (1 + B)", {
  t <- rapid_test(lm_fit(), n_permutations = 49, seed = 3)
  expect_identical(t$p_value,
                   (1 + sum(t$null_distribution >= t$statistic)) / (1 + 49))
  expect_identical(t$significant, t$p_value <= t$alpha)
})

test_that("no information gives a non-significant test, strong dependence a significant one", {
  strong <- rapid_test(lm_fit(300), n_permutations = 99, seed = 1)
  none <- rapid_test(lm_fit(300, shuffle_target = TRUE), n_permutations = 99, seed = 1)
  expect_true(strong$significant)
  expect_false(none$significant)
  expect_gt(strong$excess, none$excess)
})

test_that("an object without stored predictions is a classed error", {
  x <- lm_fit(store_predictions = FALSE)
  expect_error(rapid_test(x, n_permutations = 9), regexp = "store_predictions",
               class = "rapid_error_object")
  expect_error(rapid_threshold_select(x, n_permutations = 9), regexp = "store_predictions",
               class = "rapid_error_object")
})

test_that("scoring arguments passed to a method on an object are refused", {
  x <- lm_fit()
  expect_error(rapid_test(x, cat_tau = 0.5), regexp = "update()", fixed = TRUE,
               class = "rapid_error_input")
  expect_error(rapid_threshold_select(x, num_epsilon = 0.2), regexp = "update()",
               fixed = TRUE, class = "rapid_error_input")
  expect_error(rapid_test(x, n_permutations = 0), class = "rapid_error_input")
  expect_error(rapid_test(x, alpha = 1.5), class = "rapid_error_input")
})

test_that("NA strategy is the same in rapid() and rapid_test()", {
  o <- sim_numeric(150)
  s <- noisy_copy(o, seed = 2)
  o$y[c(3, 10, 20, 44)] <- NA
  for (st in c("constant", "drop", "median")) {
    expect_warning(r <- rapid(o, s, qi, "y", model_type = "lm", na_strategy = st, seed = 1),
                   class = "rapid_warning_missing_target")
    expect_warning(t <- rapid_test(o, s, qi, "y", model_type = "lm", na_strategy = st,
                                   n_permutations = 9, seed = 1),
                   class = "rapid_warning_missing_target")
    expect_identical(t$statistic, r$rapid, info = st)
    expect_identical(t$observed_rapid$n_total, r$n_total, info = st)
  }
  expect_warning(r <- rapid(o, s, qi, "y", model_type = "lm", na_strategy = "drop"),
                 class = "rapid_warning_missing_target")
  expect_identical(nobs(r), 146L)
})

test_that("the data interface of rapid_test() keeps a given n_perm", {
  o <- sim_numeric(100)
  s <- noisy_copy(o, seed = 2)
  t <- rapid_test(o, s, qi, "y", model_type = "lm", n_perm = 5, n_permutations = 9, seed = 1)
  expect_identical(t$observed_rapid$n_perm, 5)
  expect_false(is.na(t$observed_rapid$chance))
  # A permutation threshold works through the data interface
  tp <- rapid_test(o, s, qi, "y", model_type = "lm", threshold_rule = "permutation",
                   n_permutations = 9, seed = 1)
  expect_identical(tp$observed_rapid$threshold_rule, "permutation")
  expect_identical(tp$statistic,
                   rapid(o, s, qi, "y", model_type = "lm", threshold_rule = "permutation",
                         seed = 1)$rapid)
  # Without n_perm the fit skips the permutations, as documented
  t0 <- rapid_test(o, s, qi, "y", model_type = "lm", n_permutations = 9, seed = 1)
  expect_identical(t0$observed_rapid$n_perm, 0)
})

test_that("format, print and summary of a rapid_test", {
  t <- rapid_test(lm_fit(), n_permutations = 19, seed = 1)
  out <- utils::capture.output(res <- withVisible(print(t)))
  expect_false(res$visible)
  expect_identical(res$value, t)
  expect_identical(out, format(t))
  expect_match(format(t)[1], "RAPID permutation test")
  s <- summary(t)
  expect_s3_class(s, "summary.rapid_test")
  expect_identical(utils::capture.output(print(s)), format(s))
  expect_identical(s$chance, t$chance)
})

test_that("as.data.frame of a rapid_test has one row per permutation", {
  t <- rapid_test(lm_fit(), n_permutations = 19, seed = 1)
  df <- as.data.frame(t)
  expect_identical(names(df), c("permutation", "rapid"))
  expect_identical(df$permutation, 1:19)
  expect_identical(df$rapid, t$null_distribution)
})

# ==========================================================================
# rapid_threshold_select
# ==========================================================================

test_that("threshold selection sweeps tau for a numeric band target", {
  x <- lm_fit()
  taus <- c(0.1, 0.3, 0.5, 0.7)
  sel <- rapid_threshold_select(x, tau_range = taus, n_permutations = 19, seed = 1)
  expect_s3_class(sel, "rapid_threshold")
  expect_identical(sel$results$threshold, taus)
  expect_true(all(c("threshold", "rapid_obs", "null_quantile", "significant")
                  %in% names(sel$results)))
  # Observed RAPID per threshold equals the rescored object
  expect_equal(sel$results$rapid_obs[2], x$rapid)
  expect_equal(sel$threshold_star, max(taus[sel$results$significant]))
  expect_identical(sel$observed_rapid$rule$threshold, sel$threshold_star)
  expect_identical(as.data.frame(sel), sel$results)
})

test_that("permuted scores are computed once", {
  x <- lm_fit()
  calls <- 0
  real <- .rapid_null
  local_mocked_bindings(.rapid_null = function(...) {
    calls <<- calls + 1
    real(...)
  })
  local_mocked_bindings(.learner_fit = function(...) stop("refit"))
  sel <- rapid_threshold_select(x, n_permutations = 19, seed = 1)
  expect_identical(calls, 1)
})

test_that("point rule selects the smallest significant epsilon", {
  x <- lm_fit(num_eval_method = "point")
  eps <- c(0.001, 0.005, 0.02, 0.05, 0.1, 0.2)
  sel <- rapid_threshold_select(x, epsilon_range = eps, n_permutations = 49, seed = 1)
  expect_identical(sel$results$threshold, eps)
  sig <- sel$results$significant
  expect_true(any(sig) && !all(sig))
  expect_identical(sel$threshold_star, min(eps[sig]))
  expect_identical(sel$observed_rapid$rule$threshold, sel$threshold_star)
  expect_identical(sel$observed_rapid$rule$epsilon, sel$threshold_star)
  expect_identical(sel$rapid, sel$observed_rapid$rapid)
})

test_that("threshold selection warns about the grid the rule does not use", {
  w <- expect_warning(rapid_threshold_select(lm_fit(), epsilon_range = c(0.05, 0.1),
                                             n_permutations = 9, seed = 1),
                      class = "rapid_warning_ignored_argument")
  expect_match(conditionMessage(w), "'epsilon_range'", fixed = TRUE)
  expect_match(conditionMessage(w), "tau_range", fixed = TRUE)
  w <- expect_warning(rapid_threshold_select(lm_fit(num_eval_method = "point"),
                                             tau_range = c(0.2, 0.4), n_permutations = 9,
                                             seed = 1),
                      class = "rapid_warning_ignored_argument")
  expect_match(conditionMessage(w), "'tau_range'", fixed = TRUE)
  # The data interface passes on only the grid that was given
  o <- sim_numeric(100)
  s <- noisy_copy(o, seed = 2)
  expect_no_warning(rapid_threshold_select(o, s, qi, "y", model_type = "lm",
                                           tau_range = c(0.2, 0.4), n_permutations = 9,
                                           seed = 1))
  expect_warning(rapid_threshold_select(o, s, qi, "y", model_type = "lm",
                                        epsilon_range = 0.1, n_permutations = 9, seed = 1),
                 regexp = "'epsilon_range'", class = "rapid_warning_ignored_argument")
})

test_that("threshold_star is NA when nothing is significant", {
  sel <- rapid_threshold_select(lm_fit(shuffle_target = TRUE), n_permutations = 19,
                                seed = 1)
  expect_false(any(sel$results$significant))
  expect_identical(sel$threshold_star, NA_real_)
  expect_null(sel$observed_rapid)
  expect_identical(sel$rapid, NA_real_)
})

test_that("threshold selection: data interface, categorical target, format", {
  o <- sim_categorical(150)
  s <- to_classes(noisy_copy(sim_numeric(150), seed = 2), ref = sim_numeric(150))
  taus <- c(0.1, 0.3, 0.5)
  by_data <- rapid_threshold_select(o, s, qi, "y", model_type = "rf", tau_range = taus,
                                    n_permutations = 19, seed = 4)
  by_object <- rapid_threshold_select(rapid(o, s, qi, "y", model_type = "rf", seed = 4,
                                            n_perm = 0),
                                      tau_range = taus, n_permutations = 19, seed = 4)
  expect_identical(by_data$results, by_object$results)
  expect_true(by_data$is_categorical)
  out <- utils::capture.output(res <- withVisible(print(by_data)))
  expect_false(res$visible)
  expect_identical(out, format(by_data))
  sm <- summary(by_data)
  expect_s3_class(sm, "summary.rapid_threshold")
  expect_identical(utils::capture.output(print(sm)), format(sm))
})

test_that("plot methods return a ggplot", {
  t <- rapid_test(lm_fit(), n_permutations = 19, seed = 1)
  p <- withVisible(plot(t))
  expect_false(p$visible)
  expect_s3_class(p$value, "ggplot")
  sel <- rapid_threshold_select(lm_fit(), tau_range = c(0.1, 0.3, 0.5),
                                n_permutations = 19, seed = 1)
  p <- withVisible(plot(sel))
  expect_false(p$visible)
  expect_s3_class(p$value, "ggplot")
})

# ==========================================================================
# rapid_synthesizer_cv
# ==========================================================================

# Adds noise to numeric variables and resamples the others
.cv_synth <- function(data, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  s <- data
  for (v in names(s)) {
    if (is.numeric(s[[v]])) {
      s[[v]] <- s[[v]] + stats::rnorm(nrow(s), 0, stats::sd(s[[v]]) * 0.1)
    } else {
      s[[v]] <- sample(s[[v]], replace = TRUE)
    }
  }
  s
}

.cv_num <- function(..., d = .make_test_data(n = 100)) {
  rapid_synthesizer_cv(d$original, .cv_synth, key_vars = c("age", "gender"),
                       target_var = "income", k = 3, model_type = "lm", ...)
}

test_that("rapid_synthesizer_cv() works with a simple synthesizer", {
  cv <- .cv_num(seed = 1, num_tau = 0.2)
  expect_s3_class(cv, "rapid_cv")
  expect_true(cv$cv_summary$mean >= 0 && cv$cv_summary$mean <= 1)
  expect_true(cv$cv_summary$sd >= 0)
  expect_true(cv$cv_summary$ci_lower <= cv$cv_summary$ci_upper)
  expect_equal(cv$settings$k, 3)
  expect_false(cv$settings$is_categorical)
  expect_equal(cv$rapid, cv$cv_summary$mean)
})

test_that("CV always holds the fold table and the fold objects", {
  cv <- .cv_num(seed = 1)
  expect_identical(names(cv$folds)[1:5], c("fold", "n", "rapid", "chance", "excess"))
  expect_true(all(c("mae", "rmse", "rmae", "rrmse") %in% names(cv$folds)))
  expect_equal(nrow(cv$folds), 3)
  expect_equal(sum(cv$folds$n), 100)
  expect_length(cv$fits, 3)
  expect_true(all(vapply(cv$fits, inherits, logical(1), "rapid")))
  expect_equal(cv$folds$rapid, vapply(cv$fits, function(f) f$rapid, numeric(1)))
  expect_identical(as.data.frame(cv), cv$folds)
})

test_that("CV labels the rule actually used", {
  cv <- .cv_num(seed = 1)
  expect_identical(cv$rule, cv$fits[[1]]$rule)
  expect_identical(cv$method, "band")
  expect_identical(cv$threshold, 0.3)

  d <- .make_test_data_cat(n = 120)
  cat_cv <- function(...) {
    rapid_synthesizer_cv(d$original, .cv_synth, key_vars = c("age", "gender"),
                         target_var = "health", k = 3, model_type = "cart", seed = 1, ...)
  }
  skip_if_not_installed("rpart")
  default <- cat_cv()
  expect_identical(default$method, "RCS_marginal")
  expect_identical(default$threshold, 0.3)
  cond <- cat_cv(cat_eval_method = "RCS_conditional")
  expect_identical(cond$method, "RCS_conditional")
  expect_identical(cond$threshold, 1)
})

test_that("CV reports the mean chance level", {
  cv <- .cv_num(seed = 1, n_perm = 5)
  expect_false(anyNA(cv$folds$chance))
  expect_equal(cv$cv_summary$mean_chance, mean(cv$folds$chance))
  expect_equal(cv$cv_summary$mean_excess, mean(cv$folds$rapid - cv$folds$chance))
  expect_equal(cv$chance, cv$cv_summary$mean_chance)
  expect_equal(cv$excess, cv$cv_summary$mean_excess)
  expect_equal(cv$fits[[1]]$n_perm, 5)
})

test_that("CV has no default seed", {
  expect_null(formals(rapid_synthesizer_cv)$seed)
  expect_false(isTRUE(formals(rapid_synthesizer_cv)$verbose))
  # The data are drawn once: .make_test_data() seeds the generator itself
  d <- .make_test_data(n = 100)
  a <- .cv_num(d = d)
  b <- .cv_num(d = d)
  expect_false(isTRUE(all.equal(a$folds$mae, b$folds$mae)))
  c1 <- .cv_num(seed = 3)
  c2 <- .cv_num(seed = 3)
  expect_equal(c1$folds$rapid, c2$folds$rapid)
  expect_equal(c1$folds$chance, c2$folds$chance)
  expect_equal(c1$folds$mae, c2$folds$mae)
})

test_that("CV accepts the old names and warns about the deprecated arguments", {
  d <- .make_test_data(n = 100)
  cv <- rapid_synthesizer_cv(original_data = d$original, synthesizer = .cv_synth,
                             quasi_identifiers = c("age", "gender"),
                             sensitive_attribute = "income", k = 3, model_type = "lm",
                             seed = 1)
  expect_s3_class(cv, "rapid_cv")
  expect_identical(cv$settings$target_var, "income")
  expect_warning(.cv_num(seed = 1, return_details = TRUE), class = "rapid_warning_deprecated")
  expect_warning(.cv_num(seed = 1, return_all_records = TRUE),
                 class = "rapid_warning_deprecated")
  expect_error(rapid_synthesizer_cv(d$original, .cv_synth, key_vars = "age",
                                    target_var = "income", k = 1),
               class = "rapid_error_input")
})

test_that("CV print, summary and plot methods", {
  cv <- .cv_num(seed = 1)
  out <- utils::capture.output(res <- withVisible(print(cv)))
  expect_false(res$visible)
  expect_identical(out, format(cv))
  expect_match(out[1], "cross-validation", ignore.case = TRUE)
  sm <- summary(cv)
  expect_s3_class(sm, "summary.rapid_cv")
  expect_identical(utils::capture.output(print(sm)), format(sm))
  expect_output(print(sm), "Mean chance level")
  p <- withVisible(plot(cv))
  expect_false(p$visible)
  expect_s3_class(p$value, "ggplot")
})

# ==========================================================================
# Wald test and Clopper-Pearson interval
# ==========================================================================

.wald_fit <- function() {
  d <- .make_test_data(n = 60)
  rapid(X = d$original, Y = d$synthetic, key_vars = c("age", "gender"),
        target_var = "income", model_type = "lm", num_epsilon = 0.5, verbose = FALSE)
}

test_that("Wald z matches the formula", {
  r <- .wald_fit()
  tt <- rapid_test(r, n_permutations = 19, seed = 1)
  pi0 <- tt$null_mean
  n <- tt$observed_rapid$n_total
  z <- (tt$statistic - pi0) / sqrt(pi0 * (1 - pi0) / n)
  expect_equal(tt$wald$pi0, pi0)
  expect_equal(tt$wald$z, z)
  expect_equal(tt$wald$p_value, 1 - stats::pnorm(z))
  out <- format(tt)
  expect_true(any(grepl("Wald", out)))
  expect_true(any(grepl("use n_permutations >= 999 for the Wald test", out, fixed = TRUE)))
  tt2 <- rapid_test(r, n_permutations = 999, seed = 1)
  expect_false(any(grepl("use n_permutations >= 999", format(tt2), fixed = TRUE)))
})

test_that("Wald is NA when the null mean is 0", {
  r <- .wald_fit()
  tt <- rapid_test(r, n_permutations = 9, seed = 1)
  tt$null_mean <- 0
  tt$wald <- .rapid_wald(tt$statistic, 0, tt$observed_rapid$n_total)
  expect_true(is.na(tt$wald$z))
  expect_true(is.na(tt$wald$p_value))
  expect_true(any(grepl("not defined", format(tt))))
})

test_that("Clopper-Pearson equals binom.test", {
  r <- .wald_fit()
  expect_equal(r$n_at_risk > 0 && r$n_at_risk < r$n_total, TRUE)
  for (lev in c(0.95, 0.9)) {
    ci <- confint(r, method = "clopper-pearson", level = lev)
    ref <- as.numeric(stats::binom.test(r$n_at_risk, r$n_total, conf.level = lev)$conf.int)
    expect_equal(as.numeric(ci), ref, tolerance = 1e-12)
  }
})

test_that("Clopper-Pearson at k = 0 and k = n has bounds 0 and 1", {
  r <- .wald_fit()
  r0 <- r
  r0$n_at_risk <- 0L
  r0$rapid <- 0
  ci0 <- confint(r0, method = "clopper-pearson")
  expect_equal(ci0[1, 1], 0)
  expect_gt(ci0[1, 2], 0)
  r1 <- r
  r1$n_at_risk <- r$n_total
  r1$rapid <- 1
  ci1 <- confint(r1, method = "clopper-pearson")
  expect_equal(ci1[1, 2], 1)
  expect_lt(ci1[1, 1], 1)
})

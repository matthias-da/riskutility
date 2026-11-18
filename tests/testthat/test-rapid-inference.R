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
    r <- rapid(o, s, qi, "y", model_type = "lm", na_strategy = st, seed = 1)
    t <- rapid_test(o, s, qi, "y", model_type = "lm", na_strategy = st,
                    n_permutations = 9, seed = 1)
    expect_identical(t$statistic, r$rapid, info = st)
    expect_identical(t$observed_rapid$n_total, r$n_total, info = st)
  }
  expect_identical(nobs(rapid(o, s, qi, "y", model_type = "lm", na_strategy = "drop")),
                   146L)
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

test_that("rapid_synthesizer_cv() works with a simple synthesizer", {
  d <- .make_test_data(n = 100)

  # Simple synthesizer: add noise to numeric, resample categorical
  simple_synth <- function(data, seed = NULL) {
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

  cv <- rapid_synthesizer_cv(
    original_data = d$original,
    synthesizer = simple_synth,
    quasi_identifiers = c("age", "gender"),
    sensitive_attribute = "income",
    k = 3,
    seed = 1,
    verbose = FALSE,
    model_type = "lm",
    num_epsilon = 0.1
  )

  expect_s3_class(cv, "rapid_cv")
  expect_true(cv$cv_summary$mean >= 0 && cv$cv_summary$mean <= 1)
  expect_true(cv$cv_summary$sd >= 0)
  expect_true(cv$cv_summary$ci_lower <= cv$cv_summary$ci_upper)
  expect_equal(cv$settings$k, 3)
  expect_false(cv$settings$is_categorical)
})

test_that("rapid_synthesizer_cv() returns details when requested", {
  d <- .make_test_data(n = 100)

  simple_synth <- function(data, seed = NULL) {
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

  cv <- rapid_synthesizer_cv(
    original_data = d$original,
    synthesizer = simple_synth,
    quasi_identifiers = c("age", "gender"),
    sensitive_attribute = "income",
    k = 3,
    return_details = TRUE,
    seed = 1,
    verbose = FALSE,
    model_type = "lm",
    num_epsilon = 0.1
  )

  expect_true(!is.null(cv$cv_details))
  expect_equal(nrow(cv$cv_details), 3)
  expect_true("rapid" %in% names(cv$cv_details))
})

test_that("print.rapid_cv() runs without error", {
  d <- .make_test_data(n = 100)
  simple_synth <- function(data, seed = NULL) {
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

  cv <- rapid_synthesizer_cv(
    original_data = d$original,
    synthesizer = simple_synth,
    quasi_identifiers = c("age", "gender"),
    sensitive_attribute = "income",
    k = 3,
    seed = 1,
    verbose = FALSE,
    model_type = "lm",
    num_epsilon = 0.1
  )

  expect_output(print(cv), "Cross-Validation Results")
})

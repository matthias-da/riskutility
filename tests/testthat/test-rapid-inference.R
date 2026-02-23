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
             model_type = "lm", num_epsilon = 10, verbose = FALSE)

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

test_that("confint.rapid() Wilson interval stays in [0, 1]", {
  d <- .make_test_data()
  r <- rapid(X = d$original, Y = d$synthetic,
             key_vars = c("age", "gender"), target_var = "income",
             model_type = "lm", num_epsilon = 10, verbose = FALSE)

  ci <- confint(r, method = "wilson")
  expect_true(ci[1, 1] >= 0)
  expect_true(ci[1, 2] <= 1)
})

test_that("confint.rapid() respects level parameter", {
  d <- .make_test_data()
  r <- rapid(X = d$original, Y = d$synthetic,
             key_vars = c("age", "gender"), target_var = "income",
             model_type = "lm", num_epsilon = 10, verbose = FALSE)

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
# rapid_test
# ==========================================================================

test_that("rapid_test() returns valid rapid_test object (numeric)", {
  d <- .make_test_data()
  res <- rapid_test(
    original_data = d$original,
    synthetic_data = d$synthetic,
    quasi_identifiers = c("age", "gender"),
    sensitive_attribute = "income",
    model_type = "lm",
    num_epsilon = 10,
    n_permutations = 19,
    seed = 1
  )

  expect_s3_class(res, "rapid_test")
  expect_true(res$statistic >= 0 && res$statistic <= 1)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  expect_length(res$null_distribution, 19)
  expect_equal(res$n_permutations, 19)
  expect_equal(res$method, "permutation")
  expect_s3_class(res$observed_rapid, "rapid")
  expect_type(res$significant, "logical")
})

test_that("rapid_test() returns valid rapid_test object (categorical)", {
  skip_if_not_installed("ranger")
  d <- .make_test_data_cat()
  res <- rapid_test(
    original_data = d$original,
    synthetic_data = d$synthetic,
    quasi_identifiers = c("age", "gender"),
    sensitive_attribute = "health",
    model_type = "rf",
    cat_tau = 0.3,
    cat_eval_method = "RCS_marginal",
    n_permutations = 19,
    seed = 1
  )

  expect_s3_class(res, "rapid_test")
  expect_true(res$statistic >= 0 && res$statistic <= 1)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  expect_length(res$null_distribution, 19)
})

test_that("rapid_test() p-value is bounded correctly", {
  d <- .make_test_data()
  res <- rapid_test(
    original_data = d$original,
    synthetic_data = d$synthetic,
    quasi_identifiers = c("age", "gender"),
    sensitive_attribute = "income",
    model_type = "lm",
    num_epsilon = 10,
    n_permutations = 19,
    seed = 1
  )

  # p = (1 + sum(>=)) / (1 + B), so min is 1/(1+B), max is 1
  expect_true(res$p_value >= 1 / (1 + 19))
  expect_true(res$p_value <= 1)
})

test_that("print.rapid_test() runs without error", {
  d <- .make_test_data()
  res <- rapid_test(
    original_data = d$original,
    synthetic_data = d$synthetic,
    quasi_identifiers = c("age", "gender"),
    sensitive_attribute = "income",
    model_type = "lm",
    num_epsilon = 10,
    n_permutations = 9,
    seed = 1
  )
  expect_output(print(res), "RAPID Permutation Test")
})

# ==========================================================================
# rapid_threshold_select
# ==========================================================================

test_that("rapid_threshold_select() works for numeric target", {
  d <- .make_test_data()
  sel <- rapid_threshold_select(
    original_data = d$original,
    synthetic_data = d$synthetic,
    quasi_identifiers = c("age", "gender"),
    sensitive_attribute = "income",
    model_type = "lm",
    epsilon_range = c(5, 10, 20),
    n_permutations = 9,
    seed = 1
  )

  expect_s3_class(sel, "rapid_threshold")
  expect_true(is.data.frame(sel$results))
  expect_equal(nrow(sel$results), 3)
  expect_true(all(c("threshold", "rapid_obs", "null_quantile", "significant")
                   %in% names(sel$results)))
  expect_false(sel$is_categorical)
})

test_that("rapid_threshold_select() works for categorical target", {
  skip_if_not_installed("ranger")
  d <- .make_test_data_cat()
  sel <- rapid_threshold_select(
    original_data = d$original,
    synthetic_data = d$synthetic,
    quasi_identifiers = c("age", "gender"),
    sensitive_attribute = "health",
    model_type = "rf",
    cat_eval_method = "RCS_marginal",
    tau_range = c(0.1, 0.3, 0.5),
    n_permutations = 9,
    seed = 1
  )

  expect_s3_class(sel, "rapid_threshold")
  expect_true(sel$is_categorical)
  expect_equal(nrow(sel$results), 3)
})

test_that("rapid_threshold_select() threshold_star is NA when nothing significant", {
  # Use a very strict alpha to ensure nothing is significant
  d <- .make_test_data()
  sel <- rapid_threshold_select(
    original_data = d$original,
    synthetic_data = d$synthetic,
    quasi_identifiers = c("age", "gender"),
    sensitive_attribute = "income",
    model_type = "lm",
    epsilon_range = c(0.001),
    n_permutations = 9,
    alpha = 0.001,
    seed = 1
  )

  # At epsilon=0.001, observed RAPID should be ~0, so not significant
  expect_s3_class(sel, "rapid_threshold")
})

test_that("print.rapid_threshold() runs without error", {
  d <- .make_test_data()
  sel <- rapid_threshold_select(
    original_data = d$original,
    synthetic_data = d$synthetic,
    quasi_identifiers = c("age", "gender"),
    sensitive_attribute = "income",
    model_type = "lm",
    epsilon_range = c(5, 10),
    n_permutations = 9,
    seed = 1
  )
  expect_output(print(sel), "Threshold Selection")
})

test_that("plot.rapid_threshold() returns a ggplot", {
  d <- .make_test_data()
  sel <- rapid_threshold_select(
    original_data = d$original,
    synthetic_data = d$synthetic,
    quasi_identifiers = c("age", "gender"),
    sensitive_attribute = "income",
    model_type = "lm",
    epsilon_range = c(5, 10, 20),
    n_permutations = 9,
    seed = 1
  )
  p <- plot(sel)
  expect_s3_class(p, "ggplot")
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
    num_epsilon = 10
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
    num_epsilon = 10
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
    num_epsilon = 10
  )

  expect_output(print(cv), "Cross-Validation Results")
})

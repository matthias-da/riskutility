# Tests for RAPID (Risk of Attribute Prediction-Induced Disclosure)

test_that("rapid() works with numeric target and lm model", {
  set.seed(123)
  n <- 100

  original <- data.frame(
    age = sample(20:60, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    income = round(rnorm(n, 50000, 10000))
  )

  synthetic <- data.frame(
    age = original$age + sample(-3:3, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    income = round(original$income * runif(n, 0.9, 1.1))
  )

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = c("age", "gender"),
    target_var = "income",
    model_type = "lm",
    num_epsilon = 0.1,
    verbose = FALSE
  )

  expect_s3_class(result, "rapid")
  expect_true(result$rapid >= 0 && result$rapid <= 1)
  expect_equal(result$model_type, "lm")
  expect_true(!result$is_categorical)
  expect_equal(result$n_total, n)
  expect_true(result$n_at_risk <= result$n_total)
  expect_true(!is.null(result$model_metrics$mae))
  expect_true(!is.null(result$model_metrics$rmse))
})

test_that("rapid() works with categorical target and rf model", {
  skip_if_not_installed("ranger")

  set.seed(123)
  n <- 100

  original <- data.frame(
    age = sample(20:60, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    health = factor(sample(c("Good", "Fair", "Poor"), n,
                            replace = TRUE, prob = c(0.6, 0.3, 0.1)))
  )

  synthetic <- data.frame(
    age = original$age + sample(-3:3, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    health = factor(sample(c("Good", "Fair", "Poor"), n,
                            replace = TRUE, prob = c(0.55, 0.35, 0.1)))
  )

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = c("age", "gender"),
    target_var = "health",
    model_type = "rf",
    cat_tau = 0.3,
    cat_eval_method = "RCS_marginal",
    verbose = FALSE
  )

  expect_s3_class(result, "rapid")
  expect_true(result$rapid >= 0 && result$rapid <= 1)
  expect_equal(result$model_type, "rf")
  expect_true(result$is_categorical)
  expect_true(!is.null(result$model_metrics$accuracy))
})

test_that("rapid() works with RCS_conditional method", {
  skip_if_not_installed("ranger")

  set.seed(123)
  n <- 100

  original <- data.frame(
    age = sample(20:60, n, replace = TRUE),
    status = factor(sample(c("A", "B"), n, replace = TRUE))
  )

  synthetic <- data.frame(
    age = original$age + sample(-2:2, n, replace = TRUE),
    status = factor(sample(c("A", "B"), n, replace = TRUE))
  )

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "age",
    target_var = "status",
    model_type = "rf",
    cat_tau = 1.25,
    cat_eval_method = "RCS_conditional",
    verbose = FALSE
  )

  expect_s3_class(result, "rapid")
  expect_equal(result$method, "RCS_conditional")
})

test_that("rapid() works with NCE method", {
  skip_if_not_installed("ranger")

  set.seed(123)
  n <- 100

  original <- data.frame(
    x = rnorm(n),
    y = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )

  synthetic <- data.frame(
    x = rnorm(n),
    y = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "x",
    target_var = "y",
    model_type = "rf",
    cat_tau = 0.5,
    cat_eval_method = "NCE",
    verbose = FALSE
  )

  expect_s3_class(result, "rapid")
  expect_equal(result$method, "NCE")
})

test_that("rapid() works with cart model", {
  skip_if_not_installed("rpart")

  set.seed(123)
  n <- 100

  original <- data.frame(
    age = sample(20:60, n, replace = TRUE),
    income = round(rnorm(n, 50000, 10000))
  )

  synthetic <- data.frame(
    age = original$age + sample(-3:3, n, replace = TRUE),
    income = round(original$income * runif(n, 0.9, 1.1))
  )

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "age",
    target_var = "income",
    model_type = "cart",
    num_epsilon = 0.15,
    verbose = FALSE
  )

  expect_s3_class(result, "rapid")
  expect_equal(result$model_type, "cart")
})

test_that("rapid() keeps all records with the columns of the rule", {
  set.seed(123)
  n <- 50

  original <- data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )

  synthetic <- data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "x",
    target_var = "y",
    model_type = "lm",
    num_epsilon = 0.2,
    verbose = FALSE
  )

  expect_equal(nrow(result$records), n)
  expect_true(all(c("x", "y", "original_value", "predicted_value", "band_lo", "band_hi",
                    "g", "b", "r", "at_risk") %in% names(result$records)))
})

test_that("rapid() print method works", {
  set.seed(123)
  n <- 50

  original <- data.frame(x = rnorm(n), y = rnorm(n))
  synthetic <- data.frame(x = rnorm(n), y = rnorm(n))

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "x",
    target_var = "y",
    model_type = "lm",
    verbose = FALSE
  )

  expect_output(print(result), "RAPID")
  expect_output(print(result), "attacker lm")
})

test_that("rapid() summary method works", {
  set.seed(123)
  n <- 50

  original <- data.frame(x = rnorm(n), y = rnorm(n))
  synthetic <- data.frame(x = rnorm(n), y = rnorm(n))

  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "x",
    target_var = "y",
    model_type = "lm",
    verbose = FALSE
  )

  summ <- summary(result)
  expect_s3_class(summ, "summary.rapid")
  expect_output(print(summ), "Summary: RAPID")
})

test_that("rare summary reproduces Oscar's marking for the point rule", {
  fx <- rapid_fixture()
  x <- rapid(fx$num$o, fx$num$s, c("x1", "x2"), "y", model_type = "lm",
             num_eval_method = "point", num_scale = "identity")
  sm <- summary(x, common_share = 0.7)
  expect_equal(sm$rare_hit_rate, fx$point_rule$rare_hit_rate, tolerance = 1e-12)
  expect_equal(sm$rare_rate, fx$point_rule$rare_confidence_rate, tolerance = 1e-12)
  expect_identical(sm$n_rare, sum(fx$point_rule$rare))
  expect_equal(sm$rare_share, 0.3)
})

test_that("rare fields are NA for a rule without baseline", {
  o <- sim_numeric(100)
  x <- rapid(o, noisy_copy(o), c("x1", "x2"), "y", model_type = "lm",
             num_eval_method = "point")
  x$baseline <- NULL
  sm <- summary(x)
  expect_true(is.na(sm$n_rare) && is.na(sm$rare_hit_rate) && is.na(sm$rare_rate))
  expect_output(print(sm), "Summary: RAPID")
})

test_that("summary reports the Wilson interval, chance level and excess", {
  o <- sim_numeric(100)
  x <- rapid(o, noisy_copy(o), c("x1", "x2"), "y", model_type = "lm", seed = 1)
  tab <- summary(x, level = 0.9)$table
  ci <- confint(x, level = 0.9)
  expect_identical(names(tab), c("target", "model", "rule", "threshold", "n", "rapid",
                                 "lower", "upper", "chance", "excess"))
  expect_equal(c(tab$lower, tab$upper), unname(ci[1, ]))
  expect_equal(tab$excess, x$rapid - x$chance)
  expect_identical(tab$model, "lm")
  expect_error(summary(x, common_share = 1), class = "rapid_error_input")
})

test_that("summary has no risk_level", {
  o <- sim_numeric(100)
  x <- rapid(o, noisy_copy(o), c("x1", "x2"), "y", model_type = "lm")
  sm <- summary(x)
  expect_null(sm$risk_level)
  out <- utils::capture.output(print(sm))
  expect_false(any(grepl("[Rr]isk level|LOW|MEDIUM|HIGH", out)))
})

test_that("summary tables of two releases stack", {
  o <- sim_numeric(120)
  a <- rapid(o, noisy_copy(o, sd = 2, seed = 1), c("x1", "x2"), "y", model_type = "lm",
             seed = 1)
  b <- rapid(o, noisy_copy(o, sd = 8, seed = 1), c("x1", "x2"), "y", model_type = "lm",
             seed = 1)
  stacked <- rbind(as.data.frame(summary(a)), as.data.frame(summary(b)))
  expect_s3_class(stacked, "data.frame")
  expect_equal(nrow(stacked), 2)
  expect_equal(stacked$rapid, c(a$rapid, b$rapid))
  expect_identical(as.data.frame(summary(a)), summary(a)$table)
})

test_that("rapid() validates inputs correctly", {
  set.seed(123)
  n <- 50

  original <- data.frame(x = rnorm(n), y = rnorm(n))
  synthetic <- data.frame(x = rnorm(n), y = rnorm(n))

  # Missing variable

  expect_error(
    rapid(original, synthetic, key_vars = "nonexistent", target_var = "y"),
    "Variables missing", class = "rapid_error_input"
  )

  # lm with categorical target
  original$cat <- factor(sample(c("A", "B"), n, replace = TRUE))
  synthetic$cat <- factor(sample(c("A", "B"), n, replace = TRUE))

  expect_error(
    rapid(original, synthetic, key_vars = "x", target_var = "cat", model_type = "lm"),
    "does not support a categorical target", class = "rapid_error_model"
  )

  # logit with more than 2 levels
  original$cat3 <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  synthetic$cat3 <- factor(sample(c("A", "B", "C"), n, replace = TRUE))

  expect_error(
    rapid(original, synthetic, key_vars = "x", target_var = "cat3", model_type = "logit"),
    "exactly two classes", class = "rapid_error_model"
  )
})

test_that("rapid() works with synth_pair", {
  set.seed(123)
  n <- 50

  original <- data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )
  synthetic <- data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )

  pair <- synth_pair(original, synthetic,
                     key_vars = "x",
                     target_var = "y")

  result <- rapid(pair, model_type = "lm", verbose = FALSE)
  expect_s3_class(result, "rapid")
})

test_that("rapid() still reads the deprecated absolute epsilon type", {
  set.seed(123)
  n <- 50

  original <- data.frame(
    x = rnorm(n),
    y = rnorm(n, 100, 10)
  )
  synthetic <- data.frame(
    x = rnorm(n),
    y = rnorm(n, 100, 10)
  )

  expect_warning(
    result <- rapid(
      X = original,
      Y = synthetic,
      key_vars = "x",
      target_var = "y",
      model_type = "lm",
      num_epsilon = 5,
      num_epsilon_type = "absolute",
      verbose = FALSE
    ),
    class = "rapid_warning_deprecated"
  )

  expect_s3_class(result, "rapid")
  expect_identical(result$rule$metric, "absolute")
  expect_identical(result$rule$epsilon, 5)
  same <- rapid(original, synthetic, key_vars = "x", target_var = "y", model_type = "lm",
                num_error_metric = "absolute", num_epsilon = 5)
  expect_identical(result$at_risk, same$at_risk)
})

test_that("rapid() handles NA values in sensitive attribute", {
  set.seed(123)
  n <- 50

  original <- data.frame(
    x = rnorm(n),
    y = c(rnorm(n - 5), rep(NA, 5))
  )
  synthetic <- data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )

  # With drop strategy
  result <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "x",
    target_var = "y",
    model_type = "lm",
    na_strategy = "drop",
    verbose = FALSE
  )

  expect_s3_class(result, "rapid")
  expect_equal(result$n_total, n - 5)
  expect_identical(result$na$dropped, 46:50)
  expect_identical(result$na$strategy, "drop")

  # With constant strategy
  result2 <- rapid(
    X = original,
    Y = synthetic,
    key_vars = "x",
    target_var = "y",
    model_type = "lm",
    na_strategy = "constant",
    na_constant_value = 0,
    verbose = FALSE
  )

  expect_s3_class(result2, "rapid")
  expect_equal(result2$n_total, n)
})


# =============================================================================
# 0.3.0: band method, marginal baseline, unified core
# =============================================================================

qi <- c("x1", "x2")

test_that("defaults are band, tau 0.3, epsilon 0.1", {
  o <- sim_numeric(100)
  x <- rapid(o, noisy_copy(o), qi, "y", model_type = "lm")
  expect_s3_class(x$rule, "rapid_rule_band")
  expect_identical(x$rule$threshold, 0.3)
  expect_identical(x$rule$epsilon, 0.1)
  expect_identical(x$rule$metric, "symmetric")
  expect_identical(x$rule$delta, 0.01)
  expect_identical(x$method, "band")
  expect_identical(x$threshold, 0.3)
  expect_false(x$is_categorical)
})

test_that("categorical default is RCS_marginal with tau 0.3", {
  o <- sim_categorical(100)
  x <- rapid(o, o[sample(100), ], qi, "y", model_type = "rf", seed = 1)
  expect_s3_class(x$rule, "rapid_rule_RCS_marginal")
  expect_identical(x$threshold, 0.3)
  expect_identical(x$method, "RCS_marginal")
  expect_true(x$is_categorical)
})

test_that("RCS_conditional defaults to tau 1 and NCE to 0.5", {
  o <- sim_categorical(100)
  rc <- rapid(o, o, qi, "y", model_type = "rf", cat_eval_method = "RCS_conditional")
  expect_identical(rc$threshold, 1)
  nce <- rapid(o, o, qi, "y", model_type = "rf", cat_eval_method = "NCE")
  expect_identical(nce$threshold, 0.5)
  expect_identical(rapid(o, o, qi, "y", model_type = "rf", cat_eval_method = "NCE",
                         cat_tau = 0.7)$threshold, 0.7)
})

test_that("at_risk equals r > tau", {
  o <- sim_numeric(120)
  x <- rapid(o, noisy_copy(o), qi, "y", model_type = "lm", num_tau = 0.2)
  expect_identical(x$at_risk, x$records$r > 0.2)
  expect_identical(x$records$at_risk, x$at_risk)
  expect_identical(x$score, x$records$r)
  expect_identical(x$baseline, x$records$b)
  expect_equal(x$rapid, mean(x$at_risk))
  expect_equal(x$n_at_risk, sum(x$at_risk))
  expect_equal(x$pct_at_risk, 100 * mean(x$at_risk))
  expect_equal(x$n_total, 120)
  expect_equal(x$records$original_value, o$y)

  oc <- sim_categorical(120)
  xc <- rapid(oc, oc, qi, "y", model_type = "rf", seed = 2)
  expect_identical(xc$at_risk, xc$records$r > 0.3)
  expect_named(xc$records, c("y", "x1", "x2", "predicted_class", "g", "b", "r", "at_risk"))
})

test_that("alias-only call works", {
  o <- sim_numeric(100)
  s <- noisy_copy(o)
  a <- rapid(original_data = o, synthetic_data = s, quasi_identifiers = qi,
             sensitive_attribute = "y", model_type = "lm")
  b <- rapid(o, s, qi, "y", model_type = "lm")
  expect_s3_class(a, "rapid")
  expect_identical(a$rapid, b$rapid)
  expect_identical(a$key_vars, qi)
  expect_error(rapid(o, s, qi, "y", original_data = o), class = "rapid_error_input")
})

test_that("rapid() runs for every model with numeric, factor, character and ordered targets", {
  skip_if_not_installed("rpart")
  has_gbm <- requireNamespace("xgboost", quietly = TRUE)
  o <- sim_numeric(120)
  s <- noisy_copy(o)
  for (m in c("rf", "lm", "cart", if (has_gbm) "gbm")) {
    x <- rapid(o, s, qi, "y", model_type = m, seed = 1, nrounds = if (m == "gbm") 20)
    expect_s3_class(x, "rapid")
    expect_identical(x$model_type, m)
    expect_true(x$rapid >= 0 && x$rapid <= 1, info = m)
  }
  variants <- list(
    factor = identity,
    character = function(d) { d$y <- as.character(d$y); d },
    ordered = function(d) { d$y <- factor(d$y, levels = levels(d$y), ordered = TRUE); d }
  )
  for (v in names(variants)) {
    o3 <- variants[[v]](to_classes(o))
    s3 <- variants[[v]](to_classes(s, ref = o))
    o2 <- variants[[v]](to_classes(o, probs = 0.5))
    s2 <- variants[[v]](to_classes(s, ref = o, probs = 0.5))
    for (m in c("rf", "cart", if (has_gbm) "gbm", "logit")) {
      x <- if (m == "logit") {
        rapid(o2, s2, qi, "y", model_type = m)
      } else {
        rapid(o3, s3, qi, "y", model_type = m, seed = 1)
      }
      expect_s3_class(x, "rapid")
      expect_true(x$is_categorical, info = paste(v, m))
      expect_true(x$model_metrics$accuracy > 0.4, info = paste(v, m))
    }
  }
})

test_that("store_predictions = FALSE drops the attack", {
  o <- sim_numeric(80)
  s <- noisy_copy(o)
  x <- rapid(o, s, qi, "y", model_type = "lm")
  expect_s3_class(x$attack, "rapid_attack")
  expect_null(x$model)
  expect_null(x$attack$model)
  expect_s3_class(x$attack$pred$dist, "rapid_pdist")
  lean <- rapid(o, s, qi, "y", model_type = "lm", store_predictions = FALSE)
  expect_null(lean$attack)
  expect_identical(lean$rapid, x$rapid)
  full <- rapid(o, s, qi, "y", model_type = "lm", store_model = TRUE)
  expect_s3_class(full$model, "lm")
})

test_that("cat_tau with a numeric target warns", {
  o <- sim_numeric(80)
  expect_warning(rapid(o, noisy_copy(o), qi, "y", model_type = "lm", cat_tau = 0.5),
                 regexp = "cat_tau", class = "rapid_warning_ignored_argument")
  oc <- sim_categorical(80)
  expect_warning(rapid(oc, oc, qi, "y", model_type = "rf", num_tau = 0.5),
                 regexp = "num_tau", class = "rapid_warning_ignored_argument")
  expect_warning(rapid(o, noisy_copy(o), qi, "y", model_type = "rf", num_scale = "log"),
                 regexp = "num_scale", class = "rapid_warning_ignored_argument")
})

test_that("return_all_records is deprecated and ignored", {
  o <- sim_numeric(80)
  expect_warning(x <- rapid(o, noisy_copy(o), qi, "y", model_type = "lm",
                            return_all_records = FALSE),
                 regexp = "at_risk_only", class = "rapid_warning_deprecated")
  expect_equal(nrow(x$records), 80)
})

test_that("a misspelled argument gets a suggestion, model arguments pass", {
  o <- sim_numeric(80)
  s <- noisy_copy(o)
  expect_error(rapid(o, s, qi, "y", model_type = "lm", num_tua = 0.5),
               regexp = "Did you mean 'num_tau'", class = "rapid_error_input")
  expect_error(rapid(o, s, qi, "y", "lm", "band", 5), class = "rapid_error_input")
  x <- rapid(o, s, qi, "y", model_type = "rf", num.trees = 50, store_model = TRUE)
  expect_equal(x$model$num.trees, 50)
})

test_that("an old-style percent tolerance is read with a warning", {
  o <- sim_numeric(80)
  s <- noisy_copy(o)
  expect_warning(x <- rapid(o, s, qi, "y", model_type = "lm", num_epsilon = 10),
                 class = "rapid_warning_legacy_tolerance")
  expect_identical(x$rule$epsilon, 0.1)
  expect_error(rapid(o, s, qi, "y", model_type = "lm", num_error_metric = "absolute"),
               class = "rapid_error_tolerance")
})

test_that("the attacker line of print() names node size, number of leaves or scale", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(100)
  s <- noisy_copy(o)
  rf <- rapid(o, s, qi, "y", model_type = "rf", seed = 1)
  expect_null(rf$model)
  expect_identical(rf$attacker, "rf: 500 trees, min.node.size 5")
  expect_output(print(rf), "min.node.size 5")
  cart <- rapid(o, s, qi, "y", model_type = "cart", seed = 1)
  expect_match(cart$attacker, "^cart: [0-9]+ lea(f|ves), pruned by cross-validation$")
  expect_output(print(cart), "pruned by cross-validation")
  lm_ <- rapid(o, s, qi, "y", model_type = "lm")
  expect_match(lm_$attacker, "^lm: (identity|log) scale \\(chosen automatically\\)$")
  expect_output(print(lm_), "scale")
})

test_that("num_scale changes the lm result and is stored in the call", {
  set.seed(9)
  n <- 300
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  mk <- function() data.frame(x1 = x1, x2 = x2,
                              y = exp(10 + 0.8 * (x1 + x2) / sqrt(2) + 0.1 * stats::rnorm(n)))
  o <- mk()
  s <- mk()
  log_ <- rapid(o, s, qi, "y", model_type = "lm", num_scale = "log")
  id <- rapid(o, s, qi, "y", model_type = "lm", num_scale = "identity")
  auto <- rapid(o, s, qi, "y", model_type = "lm")
  expect_gt(log_$rapid, id$rapid)
  expect_identical(auto$rapid, log_$rapid)
  expect_identical(log_$call$num_scale, "log")
  expect_identical(log_$attacker, "lm: log scale")
  expect_identical(auto$attacker, "lm: log scale (chosen automatically)")
})

test_that("a quasi-identifier named like a record column keeps its name", {
  o <- sim_numeric(80)
  names(o)[1:2] <- c("g", "b")
  x <- rapid(o, noisy_copy(o), c("g", "b"), "y", model_type = "lm")
  expect_equal(x$records$g, o$g)
  expect_equal(x$records$b, o$b)
  expect_true(all(c("g.rapid", "b.rapid", "r") %in% names(x$records)))
  expect_identical(x$records$b.rapid, x$baseline)
  expect_identical(x$records$r, x$score)
})

test_that("the same seed gives the same result", {
  o <- sim_numeric(100)
  s <- noisy_copy(o)
  a <- rapid(o, s, qi, "y", model_type = "rf", seed = 5)
  b <- rapid(o, s, qi, "y", model_type = "rf", seed = 5)
  expect_identical(a$score, b$score)
})

test_that("plot.rapid() draws every type for numeric and categorical targets", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  o <- sim_numeric(150)
  o$grp <- factor(rep(c("u", "v", "w"), 50))
  s <- noisy_copy(o)
  key <- c("x1", "x2", "grp")
  band <- rapid(o, s, key, "y", model_type = "rf", seed = 1, store_model = TRUE)
  point <- rapid(o, s, key, "y", model_type = "lm", num_eval_method = "point")
  cat3 <- rapid(to_classes(o), to_classes(s, ref = o), key, "y", model_type = "rf",
                seed = 1, store_model = TRUE)
  for (x in list(band, point, cat3)) {
    for (w in 1:3) expect_no_error(plot(x, which = w))
    expect_no_error(suppressWarnings(plot(x, which = 5)))
    expect_no_error(suppressWarnings(plot(x, which = 6)))
    expect_s3_class(plot(x, facet_by = "grp"), "ggplot")
  }
  expect_no_error(plot(band, which = 4))
  expect_no_error(plot(cat3, which = 4))
  expect_error(plot(point, which = 4), regexp = "store_model", class = "rapid_error_object")
})

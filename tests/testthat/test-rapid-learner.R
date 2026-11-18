# Tests for user-defined attackers: rapid_learner() and its learner methods

qi <- c("x1", "x2")

# A custom learner wrapping lm with the normal predictive distribution
lm_learner <- function(quantiles = TRUE) {
  rapid_learner(
    "my_lm", fit = stats::lm,
    predict_point = function(fit, newdata) as.numeric(stats::predict(fit, newdata)),
    predict_quantiles = if (quantiles) {
      function(fit, newdata, probs) {
        p <- stats::predict(fit, newdata, se.fit = TRUE)
        sd <- sqrt(p$se.fit^2 + p$residual.scale^2)
        outer(as.numeric(p$fit), rep(1, length(probs))) + outer(sd, stats::qnorm(probs))
      }
    }
  )
}

# A custom classifier wrapping rpart; `...` takes the arguments of rapid()
cart_learner <- function(predict_prob = NULL) {
  if (is.null(predict_prob)) {
    predict_prob <- function(fit, newdata) stats::predict(fit, newdata, type = "prob")
  }
  rapid_learner(
    "my_cart",
    fit = function(formula, data, ...) rpart::rpart(formula, data, method = "class", ...),
    predict_prob = predict_prob, requires = "rpart"
  )
}

# ---- The constructor ---------------------------------------------------------------

test_that("the constructor derives the targets and checks the functions", {
  l <- lm_learner()
  expect_s3_class(l, c("rapid_learner_custom", "rapid_learner"), exact = TRUE)
  expect_identical(l$name, "my_lm")
  expect_identical(l$targets, "numeric")
  expect_false(l$unseen_levels)
  expect_identical(cart_learner()$targets, "categorical")

  both <- rapid_learner("both", fit = function(formula, data) NULL,
                        predict_prob = function(fit, newdata) NULL,
                        predict_point = function(fit, newdata) NULL)
  expect_identical(both$targets, c("numeric", "categorical"))

  expect_error(rapid_learner("x", fit = stats::lm), class = "rapid_error_input")
  expect_error(rapid_learner(c("a", "b"), fit = stats::lm,
                             predict_point = function(fit, newdata) 1),
               class = "rapid_error_input")
  expect_error(rapid_learner("x", fit = "lm", predict_point = function(fit, newdata) 1),
               class = "rapid_error_input")
  expect_error(rapid_learner("x", fit = function(f, d) NULL,
                             predict_point = function(fit, newdata) 1),
               class = "rapid_error_input", regexp = "formula")
  expect_error(rapid_learner("x", fit = stats::lm, predict_point = function(object) 1),
               class = "rapid_error_input", regexp = "newdata")
  expect_error(rapid_learner("x", fit = stats::lm, predict_point = function(fit, newdata) 1,
                             predict_quantiles = function(fit, newdata) 1),
               class = "rapid_error_input", regexp = "probs")
  expect_error(rapid_learner("x", fit = stats::lm, predict_prob = function(fit, newdata) 1,
                             predict_quantiles = function(fit, newdata, probs) 1),
               class = "rapid_error_input", regexp = "predict_point")
  expect_error(rapid_learner("x", fit = stats::lm, predict_point = function(fit, newdata) 1,
                             requires = 1), class = "rapid_error_input")
  # Functions with `...` take any argument
  expect_s3_class(rapid_learner("x", fit = function(...) NULL,
                                predict_point = function(...) NULL), "rapid_learner")
})

test_that("print shows name and targets", {
  out <- capture.output(res <- print(cart_learner()))
  expect_identical(res$name, "my_cart")
  expect_true(any(grepl("my_cart", out)))
  expect_true(any(grepl("categorical", out)))
  expect_true(any(grepl("rpart", out)))
  expect_identical(format(cart_learner()), out)
})

# ---- Equivalence with the built-in learners -------------------------------------------

test_that("a custom learner wrapping lm reproduces model_type = 'lm' with num_dist = 'grid'", {
  o <- sim_numeric(120)
  s <- noisy_copy(o, seed = 2)
  built_in <- rapid(o, s, qi, "y", model_type = "lm", num_dist = "grid",
                    num_scale = "identity", n_perm = 3, seed = 1)
  custom <- rapid(o, s, qi, "y", model_type = lm_learner(), num_dist = "grid",
                  n_perm = 3, seed = 1)
  expect_identical(custom$rapid, built_in$rapid)
  expect_identical(custom$at_risk, built_in$at_risk)
  expect_equal(custom$score, built_in$score)
  expect_identical(custom$model_type, "my_lm")

  # "auto" is the grid as well; "exact" does not exist for a custom learner
  auto <- rapid(o, s, qi, "y", model_type = lm_learner(), n_perm = 3, seed = 1)
  expect_identical(auto$at_risk, custom$at_risk)
  expect_error(rapid(o, s, qi, "y", model_type = lm_learner(), num_dist = "exact"),
               class = "rapid_error_model")
})

test_that("a custom classifier wrapping rpart reproduces model_type = 'cart'", {
  skip_if_not_installed("rpart")
  o <- sim_categorical(150)
  s <- noisy_copy(sim_numeric(150), seed = 3)
  s <- to_classes(s, ref = sim_numeric(150))
  built_in <- rapid(o, s, qi, "y", model_type = "cart", prune = FALSE, cp = 0.01,
                    minbucket = 7, xval = 0, n_perm = 3, seed = 1)
  custom <- rapid(o, s, qi, "y", model_type = cart_learner(), cp = 0.01, minbucket = 7,
                  xval = 0, n_perm = 3, seed = 1)
  expect_identical(custom$rapid, built_in$rapid)
  expect_identical(custom$at_risk, built_in$at_risk)
})

test_that("a class the synthetic data lack gets probability 0", {
  skip_if_not_installed("rpart")
  o <- sim_categorical(150)
  s <- noisy_copy(o, seed = 3)
  s <- s[s$y != "c", ]
  r <- rapid(o, s, qi, "y", model_type = cart_learner(), n_perm = 2, seed = 1,
             store_model = TRUE)
  expect_identical(colnames(r$attack$pred), levels(o$y))
  expect_true(all(r$attack$pred[, "c"] == 0))
})

# ---- Validation of the output of the user's functions ----------------------------------

test_that("band rule without predict_quantiles is a classed error that names the missing function", {
  o <- sim_numeric(80)
  s <- noisy_copy(o, seed = 2)
  expect_error(rapid(o, s, qi, "y", model_type = lm_learner(quantiles = FALSE), n_perm = 2),
               class = "rapid_error_model", regexp = "predict_quantiles")
  # The point rule needs no distribution
  p <- rapid(o, s, qi, "y", model_type = lm_learner(quantiles = FALSE),
             num_eval_method = "point", n_perm = 2, seed = 1)
  expect_s3_class(p, "rapid")
  expect_error(update(p, num_eval_method = "band"), class = "rapid_error_model",
               regexp = "predict_quantiles")
})

test_that("a probability matrix without class names is a classed error", {
  o <- sim_categorical(80)
  s <- noisy_copy(o, seed = 2)
  constant <- function(formula, data) table(data$y) / nrow(data)
  unnamed <- rapid_learner("unnamed", fit = constant, predict_prob = function(fit, newdata) {
    unname(matrix(as.numeric(fit), nrow(newdata), length(fit), byrow = TRUE))
  })
  expect_error(rapid(o, s, qi, "y", model_type = unnamed, n_perm = 2),
               class = "rapid_error_model", regexp = "class names")
})

test_that("the output of the user's functions is checked", {
  o <- sim_numeric(80)
  s <- noisy_copy(o, seed = 2)
  bad_point <- function(f) rapid_learner("bad", fit = stats::lm, predict_point = f)
  expect_error(rapid(o, s, qi, "y", model_type = bad_point(function(fit, newdata) 1),
                     num_eval_method = "point", n_perm = 0),
               class = "rapid_error_model", regexp = "predict_point")
  expect_error(rapid(o, s, qi, "y",
                     model_type = bad_point(function(fit, newdata) rep(NA_real_, nrow(newdata))),
                     num_eval_method = "point", n_perm = 0),
               class = "rapid_error_model", regexp = "NA")
  bad_q <- rapid_learner("bad", fit = stats::lm,
                         predict_point = function(fit, newdata) rep(1, nrow(newdata)),
                         predict_quantiles = function(fit, newdata, probs) {
                           matrix(1, nrow(newdata), length(probs) - 1)
                         })
  expect_error(rapid(o, s, qi, "y", model_type = bad_q, n_perm = 0),
               class = "rapid_error_model", regexp = "predict_quantiles")

  oc <- sim_categorical(80)
  sc <- noisy_copy(oc, seed = 2)
  constant <- function(formula, data) table(data$y) / nrow(data)
  wrong_names <- rapid_learner("wrong", fit = constant, predict_prob = function(fit, newdata) {
    matrix(1 / 3, nrow(newdata), 3, dimnames = list(NULL, c("a", "b", "z")))
  })
  expect_error(rapid(oc, sc, qi, "y", model_type = wrong_names, n_perm = 0),
               class = "rapid_error_model", regexp = "z")

  # A target type the learner does not support
  expect_error(rapid(oc, sc, qi, "y", model_type = lm_learner(), n_perm = 0),
               class = "rapid_error_model")
  # Errors of the user's own functions are passed on
  failing <- rapid_learner("failing", fit = function(formula, data) stop("my own error"),
                           predict_point = function(fit, newdata) 1)
  expect_error(rapid(o, s, qi, "y", model_type = failing), "my own error")
})

test_that("a quasi-identifier level the synthetic data lack is the classed level error", {
  o <- sim_numeric(80)
  o$g <- factor(rep(c("u", "v"), 40))
  s <- noisy_copy(o, seed = 2)
  s$g <- factor(rep("u", 80), levels = "u")
  expect_error(rapid(o, s, c("x1", "g"), "y", model_type = lm_learner(), n_perm = 0),
               class = "rapid_error_levels")
})

# ---- Entry points ----------------------------------------------------------------------

test_that("an object of a custom learner supports predict(), update() and plots", {
  o <- sim_numeric(100)
  s <- noisy_copy(o, seed = 2)
  r <- rapid(o, s, qi, "y", model_type = lm_learner(), n_perm = 2, seed = 1,
             store_model = TRUE)
  expect_s3_class(r$model, "lm")
  p <- predict(r, newdata = o[1:5, ])
  expect_identical(nrow(as.data.frame(p)), 5L)
  u <- update(r, num_tau = 0.5)
  expect_identical(u$model_type, "my_lm")
  refit <- update(r, n_quantiles = 50)
  expect_identical(refit$model_type, "my_lm")
  expect_true(any(grepl("my_lm", format(r))))
  expect_message(plot(r, which = 4), "not available")

  no_model <- rapid(o, s, qi, "y", model_type = lm_learner(), n_perm = 2, seed = 1)
  expect_null(no_model$model)
})

test_that("a custom learner works inside rapid_suite()", {
  o <- sim_numeric(100)
  s <- noisy_copy(o, seed = 2)
  suite <- rapid_suite(o, s, qi, "y", model_types = list("lm", lm_learner()), n_perm = 2,
                       seed = 1)
  expect_identical(names(suite$fits), c("lm", "my_lm"))
  expect_identical(suite$results$model_type, c("lm", "my_lm"))
  expect_identical(suite[["my_lm"]]$at_risk,
                   rapid(o, s, qi, "y", model_type = lm_learner(), n_perm = 2,
                         seed = 1)$at_risk)

  single <- rapid_suite(o, s, qi, "y", model_types = lm_learner(), n_perm = 2, seed = 1)
  expect_identical(names(single$fits), "my_lm")

  # A learner that needs a missing package is skipped like a built-in one
  missing_pkg <- rapid_learner("needs_pkg", fit = stats::lm,
                               predict_point = function(fit, newdata) 1,
                               requires = "a.package.that.does.not.exist")
  expect_warning(
    rapid_suite(o, s, qi, "y", model_types = list(lm_learner(), missing_pkg), n_perm = 2),
    class = "rapid_warning_skipped_model")

  expect_error(rapid_suite(o, s, qi, "y", model_types = list(lm_learner(), lm_learner())),
               class = "rapid_error_input", regexp = "my_lm")
  expect_error(rapid_suite(o, s, qi, "y", model_types = list("lm", 1)),
               class = "rapid_error_input")
})

test_that("the companions accept a custom learner", {
  o <- sim_numeric(80)
  s <- noisy_copy(o, seed = 2)
  l <- lm_learner()
  t <- rapid_test(o, s, qi, "y", model_type = l, n_permutations = 9, seed = 1)
  expect_s3_class(t, "rapid_test")
  th <- rapid_threshold_select(o, s, qi, "y", model_type = l, tau_range = c(0.2, 0.4),
                               n_permutations = 9, seed = 1)
  expect_s3_class(th, "rapid_threshold")
  qc <- rapid_qi_contribution(o, s, qi, "y", model_type = l, n_permutations = 2, seed = 1)
  expect_s3_class(qc, "rapid_contribution")
  cv <- rapid_synthesizer_cv(o, function(data, seed = NULL) noisy_copy(data, seed = seed),
                             key_vars = qi, target_var = "y", model_type = l, k = 2,
                             n_perm = 2, seed = 1)
  expect_identical(cv$model_type, "my_lm")
})

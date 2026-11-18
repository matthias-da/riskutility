# Input layer of RAPID: conditions, aliases, argument checks, data preparation

qi <- c("x1", "x2")

test_that("conditions carry their classes", {
  e <- tryCatch(.rapid_abort(c("a", "b"), "input"), error = function(e) e)
  expect_s3_class(e, c("rapid_error_input", "rapid_error", "error", "condition"), exact = TRUE)
  expect_identical(conditionMessage(e), "ab")
  w <- tryCatch(.rapid_warn("careful", "deprecated"), warning = function(w) w)
  expect_s3_class(w, c("rapid_warning_deprecated", "rapid_warning", "warning", "condition"),
                  exact = TRUE)
})

test_that("aliases and synth_pair resolve to the same inputs", {
  o <- sim_numeric(30)
  s <- noisy_copy(o)
  direct <- .rapid_inputs(o, s, qi, "y", list())
  alias <- .rapid_inputs(NULL, NULL, NULL, NULL,
                         list(original_data = o, synthetic_data = s,
                              quasi_identifiers = qi, sensitive_attribute = "y"))
  pair <- .rapid_inputs(synth_pair(o, s, key_vars = qi, target_var = "y"),
                        NULL, NULL, NULL, list())
  expect_named(direct, c("X", "Y", "key_vars", "target_var", "holdout", "dots"))
  expect_identical(alias, direct)
  expect_identical(pair, direct)
  expect_null(direct$holdout)
})

test_that("a synth_pair holdout is passed on", {
  o <- sim_numeric(30)
  s <- noisy_copy(o)
  h <- sim_numeric(10, seed = 5)
  sp <- synth_pair(o, s, key_vars = qi, target_var = "y", holdout = h)
  expect_identical(.rapid_inputs(sp, NULL, NULL, NULL, list())$holdout, h)
})

test_that("explicit key_vars and target_var override those of a synth_pair", {
  o <- sim_numeric(30)
  s <- noisy_copy(o)
  sp <- synth_pair(o, s, key_vars = qi, target_var = "y")
  inp <- .rapid_inputs(sp, NULL, "x1", "x2", list())
  expect_identical(inp$key_vars, "x1")
  expect_identical(inp$target_var, "x2")
  expect_error(.rapid_inputs(synth_pair(o, s), NULL, NULL, NULL, list()),
               class = "rapid_error_input")
  expect_error(.rapid_inputs(sp, s, NULL, NULL, list()), class = "rapid_error_input")
})

test_that("a name and its alias together is an error", {
  o <- sim_numeric(30)
  s <- noisy_copy(o)
  expect_error(.rapid_inputs(o, s, qi, "y", list(original_data = o)),
               class = "rapid_error_input")
  expect_error(.rapid_inputs(o, s, qi, "y", list(sensitive_attribute = "y")),
               class = "rapid_error_input")
})

test_that("other arguments in the dots are kept", {
  o <- sim_numeric(30)
  s <- noisy_copy(o)
  inp <- .rapid_inputs(o, NULL, qi, "y", list(synthetic_data = s, num.trees = 50))
  expect_identical(inp$dots, list(num.trees = 50))
})

test_that("missing variables are named in the error", {
  o <- sim_numeric(30)
  s <- noisy_copy(o)
  expect_error(.rapid_inputs(o, s, c("x1", "zz"), "y", list()),
               regexp = "zz", class = "rapid_error_input")
  s2 <- s[, c("x1", "x2")]
  err <- tryCatch(.rapid_inputs(o, s2, qi, "y", list()), error = function(e) e)
  expect_s3_class(err, "rapid_error_input")
  expect_match(conditionMessage(err), "'y'")
  expect_match(conditionMessage(err), "Y")
})

test_that("missing or malformed data arguments are classed errors", {
  o <- sim_numeric(30)
  s <- noisy_copy(o)
  expect_error(.rapid_inputs(NULL, s, qi, "y", list()), class = "rapid_error_input")
  expect_error(.rapid_inputs(o, NULL, qi, "y", list()), class = "rapid_error_input")
  expect_error(.rapid_inputs(as.matrix(o), s, qi, "y", list()), class = "rapid_error_input")
  expect_error(.rapid_inputs(o, s, NULL, "y", list()), class = "rapid_error_input")
  expect_error(.rapid_inputs(o, s, qi, c("y", "x1"), list()), class = "rapid_error_input")
  expect_error(.rapid_inputs(o, s, c("x1", "y"), "y", list()), class = "rapid_error_input")
  expect_error(.rapid_inputs(o[0, ], s, qi, "y", list()), class = "rapid_error_input")
})

test_that("a misspelled argument gets a suggestion", {
  known <- c("X", "Y", "key_vars", "target_var", "model_type", "num_tau", "cat_tau",
             "num_epsilon", "seed", "verbose")
  expect_error(.rapid_check_dots(list(num_tua = 0.5), known),
               regexp = "Did you mean 'num_tau'\\?$", class = "rapid_error_input")
  expect_error(.rapid_check_dots(list(model.type = "lm"), known),
               regexp = "Did you mean 'model_type'", class = "rapid_error_input")
  expect_silent(.rapid_check_dots(list(num.trees = 50), known))
})

test_that("model arguments near a one-letter argument are not refused", {
  known <- c("X", "Y", "key_vars", "target_var", "model_type", "num_tau", "seed")
  expect_silent(.rapid_check_dots(
    list(y = TRUE, cp = 0.01, x = TRUE, mtry = 2, eta = 0.1, max_depth = 2, maxdepth = 2,
         min.node.size = 5, nrounds = 20), known))
})

test_that("unnamed extra arguments are an error", {
  known <- c("X", "Y", "num_tau")
  expect_error(.rapid_check_dots(list(5), known), class = "rapid_error_input")
  expect_error(.rapid_check_dots(list(num.trees = 50, "lm"), known),
               class = "rapid_error_input")
  expect_silent(.rapid_check_dots(list(), known))
})

test_that("arguments passed to a method on an object are refused with the hint", {
  expect_silent(.rapid_check_empty_dots(list(), "use update()"))
  err <- tryCatch(.rapid_check_empty_dots(list(cat_tau = 0.5, 3), "use update()"),
                  error = function(e) e)
  expect_s3_class(err, "rapid_error_input")
  expect_match(conditionMessage(err), "cat_tau")
  expect_match(conditionMessage(err), "use update()", fixed = TRUE)
})

test_that("character and ordered targets are harmonised", {
  o <- sim_categorical(60)
  s <- to_classes(noisy_copy(sim_numeric(60)), ref = sim_numeric(60))
  as_chr <- function(d) { d$y <- as.character(d$y); d }
  as_ord <- function(d) { d$y <- factor(d$y, levels = c("a", "b", "c"), ordered = TRUE); d }
  cases <- list(character = list(as_chr(o), as_chr(s)),
                ordered = list(as_ord(o), as_ord(s)),
                mixed = list(as_ord(o), as_chr(s)))
  for (nm in names(cases)) {
    p <- .rapid_prepare_data(cases[[nm]][[1]], cases[[nm]][[2]], qi, "y", "constant", 0)
    expect_true(p$is_categorical, info = nm)
    expect_true(is.factor(p$X$y) && !is.ordered(p$X$y), info = nm)
    expect_true(is.factor(p$Y$y) && !is.ordered(p$Y$y), info = nm)
    expect_identical(levels(p$X$y), levels(p$Y$y), info = nm)
    expect_identical(as.character(p$X$y), as.character(o$y), info = nm)
    expect_identical(as.character(p$Y$y), as.character(s$y), info = nm)
  }
})

test_that("target levels follow the original, extra synthetic levels come last", {
  o <- data.frame(x1 = 1:6, y = factor(c("b", "a", "b", "a", "b", "a"), levels = c("b", "a")))
  s <- data.frame(x1 = 1:6, y = c("a", "z", "b", "a", "z", "b"), stringsAsFactors = FALSE)
  p <- .rapid_prepare_data(o, s, "x1", "y", "constant", 0)
  expect_identical(levels(p$X$y), c("b", "a", "z"))
  expect_identical(levels(p$Y$y), c("b", "a", "z"))
  expect_false(anyNA(p$Y$y))
})

test_that("integer and double targets are harmonised", {
  o <- data.frame(x1 = c(1L, 2L, 3L, 4L), y = c(10L, 20L, 30L, 40L))
  s <- data.frame(x1 = c(1.5, 2, 3, 4), y = c(11.5, 19, 31, 42))
  p <- .rapid_prepare_data(o, s, "x1", "y", "constant", 0)
  expect_false(p$is_categorical)
  expect_type(p$X$y, "double")
  expect_type(p$Y$y, "double")
  expect_type(p$X$x1, "double")
  expect_type(p$Y$x1, "double")
  expect_named(p$X, c("y", "x1"))
})

test_that("character and factor quasi-identifiers get common levels", {
  o <- data.frame(g = c("a", "b", "c", "a"), y = c(1, 2, 3, 4), stringsAsFactors = FALSE)
  s <- data.frame(g = factor(c("a", "b", "b", "d")), y = c(1, 2, 3, 4))
  p <- .rapid_prepare_data(o, s, "g", "y", "constant", 0)
  expect_identical(levels(p$X$g), c("a", "b", "c", "d"))
  expect_identical(levels(p$Y$g), c("a", "b", "c", "d"))
  expect_identical(p$unseen, list(g = "c"))
})

test_that("incompatible column types are a classed error naming the variable", {
  o <- data.frame(x1 = 1:4, y = c(1, 2, 3, 4))
  s <- data.frame(x1 = 1:4, y = factor(c("a", "b", "a", "b")))
  expect_error(.rapid_prepare_data(o, s, "x1", "y", "constant", 0),
               regexp = "'y'", class = "rapid_error_input")
  s2 <- data.frame(x1 = c("u", "v", "u", "v"), y = c(1, 2, 3, 4), stringsAsFactors = FALSE)
  expect_error(.rapid_prepare_data(o, s2, "x1", "y", "constant", 0),
               regexp = "'x1'", class = "rapid_error_input")
})

test_that("missing target values follow the NA strategy", {
  o <- data.frame(x1 = 1:6, y = c(1, NA, 3, 4, NA, 6))
  s <- data.frame(x1 = 1:6, y = c(2, 3, NA, 5, 6, 7))
  drop <- .rapid_prepare_data(o, s, "x1", "y", "drop", 0)
  expect_identical(drop$dropped, c(2L, 5L))
  expect_equal(drop$X$y, c(1, 3, 4, 6))
  expect_equal(nrow(drop$Y), 5)
  const <- .rapid_prepare_data(o, s, "x1", "y", "constant", -1)
  expect_identical(const$dropped, integer(0))
  expect_equal(const$X$y, c(1, -1, 3, 4, -1, 6))
  med <- .rapid_prepare_data(o, s, "x1", "y", "median", 0)
  expect_equal(med$X$y, c(1, 3.5, 3, 4, 3.5, 6))

  oc <- data.frame(x1 = 1:4, y = c("a", NA, "b", "a"), stringsAsFactors = FALSE)
  sc <- data.frame(x1 = 1:4, y = c("a", "b", "b", "a"), stringsAsFactors = FALSE)
  pc <- .rapid_prepare_data(oc, sc, "x1", "y", "drop", 0)
  expect_identical(as.character(pc$X$y), c("a", "missing", "b", "a"))
  expect_identical(levels(pc$X$y), levels(pc$Y$y))
  expect_identical(pc$dropped, integer(0))
})

test_that("numeric NA in a quasi-identifier is a classed error", {
  o <- data.frame(x1 = c(1, NA, 3, 4), y = c(1, 2, 3, 4))
  s <- data.frame(x1 = c(1, 2, 3, 4), y = c(1, 2, 3, 4))
  expect_error(.rapid_prepare_data(o, s, "x1", "y", "constant", 0),
               regexp = "'x1'", class = "rapid_error_input")
})

test_that("unseen quasi-identifier levels never give a raw backend error", {
  skip_if_not_installed("rpart")
  set.seed(3)
  n <- 120
  o <- data.frame(x1 = stats::rnorm(n), g = factor(sample(c("a", "b", "c"), n, TRUE)))
  o$y <- 50 + 5 * o$x1 + 3 * as.integer(o$g) + stats::rnorm(n)
  s <- o[o$g != "c", ]
  s$g <- droplevels(s$g)
  for (m in c("lm", "cart", "rf")) {
    res <- tryCatch(
      suppressWarnings(rapid(o, s, key_vars = c("x1", "g"), target_var = "y",
                             model_type = m)),
      error = function(e) e
    )
    if (inherits(res, "error")) {
      expect_s3_class(res, "rapid_error_levels")
      expect_match(conditionMessage(res), "'g'", info = m)
      expect_match(conditionMessage(res), "c", info = m)
    } else {
      expect_s3_class(res, "rapid")
      expect_equal(res$n_total, n, info = m)
    }
  }
  # A linear model has no coefficient for a level it has never seen
  expect_error(suppressWarnings(rapid(o, s, key_vars = c("x1", "g"), target_var = "y",
                                      model_type = "lm")),
               class = "rapid_error_levels")
})

test_that("tolerance rules", {
  expect_identical(.rapid_resolve_tolerance("symmetric", 0.1, NULL),
                   list(metric = "symmetric", epsilon = 0.1))
  expect_identical(.rapid_resolve_tolerance("stabilised_relative", NULL, NULL)$epsilon, 0.1)
  expect_warning(res <- .rapid_resolve_tolerance("symmetric", 10, NULL),
                 class = "rapid_warning_legacy_tolerance")
  expect_identical(res, list(metric = "symmetric", epsilon = 0.1))
  expect_warning(res <- .rapid_resolve_tolerance("symmetric", 5, "percentage"),
                 class = "rapid_warning_deprecated")
  expect_identical(res, list(metric = "symmetric", epsilon = 0.05))
  expect_warning(res <- .rapid_resolve_tolerance("symmetric", 500, "absolute"),
                 class = "rapid_warning_deprecated")
  expect_identical(res, list(metric = "absolute", epsilon = 500))
  expect_identical(.rapid_resolve_tolerance("absolute", 2, NULL),
                   list(metric = "absolute", epsilon = 2))
  expect_error(.rapid_resolve_tolerance("absolute", NULL, NULL),
               class = "rapid_error_tolerance")
  expect_error(suppressWarnings(.rapid_resolve_tolerance("absolute", 5, "percentage")),
               class = "rapid_error_tolerance")
  expect_error(.rapid_resolve_tolerance("symmetric", 0, NULL), class = "rapid_error_tolerance")
  expect_error(.rapid_resolve_tolerance("symmetric", -0.1, NULL),
               class = "rapid_error_tolerance")
  expect_error(.rapid_resolve_tolerance("symmetric", "0.1", NULL),
               class = "rapid_error_tolerance")
  expect_error(.rapid_resolve_tolerance("symmetric", c(0.1, 0.2), NULL),
               class = "rapid_error_tolerance")
})

test_that("default thresholds per method", {
  expect_identical(.rapid_default_tau("band"), 0.3)
  expect_identical(.rapid_default_tau("RCS_marginal"), 0.3)
  expect_identical(.rapid_default_tau("RCS_conditional"), 1)
  expect_identical(.rapid_default_tau("NCE"), 0.5)
})

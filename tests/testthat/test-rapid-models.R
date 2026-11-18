# Learner classes (attacker models) and predictive distributions of RAPID

qi <- c("x1", "x2")
f_num <- y ~ x1 + x2

# Fits a learner on the synthetic data and predicts the original records.
attack_fit <- function(model, o, s, args = list(), num_dist = "auto", num_scale = "auto",
                       metric = "symmetric", formula = f_num) {
  prep <- .rapid_prepare_data(o, s, setdiff(all.vars(formula), "y"), "y", "constant", 0)
  learner <- .rapid_learner(model, args)
  fit <- .learner_fit(learner, formula, prep$Y, "y", num_dist, num_scale, metric)
  list(learner = learner, fit = fit, X = prep$X, Y = prep$Y,
       point = .learner_predict(learner, fit, prep$X, "y"))
}

pdist_of <- function(a, num_dist = "auto", n_quantiles = 100) {
  .learner_pdist(a$learner, a$fit, a$X, a$Y, "y", num_dist, n_quantiles)
}

# Positive target with multiplicative noise: log(y) = 10 + 0.8 * s
sim_lognormal <- function(n, rho, seed) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  s <- rho * (x1 + x2) / sqrt(2) + sqrt(1 - rho^2) * stats::rnorm(n)
  data.frame(x1 = x1, x2 = x2, s = s, y = exp(10 + 0.8 * s))
}

test_that("all five learners fit and predict", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(120)
  s <- noisy_copy(o)
  models <- c("rf", "lm", "cart", if (requireNamespace("xgboost", quietly = TRUE)) "gbm")
  for (m in models) {
    a <- attack_fit(m, o, s)
    expect_type(a$point, "double")
    expect_length(a$point, nrow(o))
    expect_false(anyNA(a$point), info = m)
  }

  oc <- sim_categorical(120)
  sc <- to_classes(s, ref = o)
  o2 <- to_classes(o, probs = 0.5)
  s2 <- to_classes(s, ref = o, probs = 0.5)
  cases <- list(rf = list(oc, sc), cart = list(oc, sc), logit = list(o2, s2))
  if (requireNamespace("xgboost", quietly = TRUE)) {
    cases$gbm <- list(oc, sc)
    cases$gbm2 <- list(o2, s2)
  }
  for (nm in names(cases)) {
    a <- attack_fit(sub("2$", "", nm), cases[[nm]][[1]], cases[[nm]][[2]])
    P <- a$point
    expect_true(is.matrix(P) && is.numeric(P), info = nm)
    expect_identical(colnames(P), levels(a$X$y), info = nm)
    expect_equal(nrow(P), nrow(oc), info = nm)
    expect_equal(unname(rowSums(P)), rep(1, nrow(P)), tolerance = 1e-6, info = nm)
    expect_true(mean(colnames(P)[max.col(P, "first")] == a$X$y) > 0.5, info = nm)
  }
})

test_that("gbm class probabilities belong to their class labels", {
  skip_if_not_installed("xgboost")
  set.seed(7)
  n <- 240
  d <- data.frame(x1 = stats::runif(n), x2 = stats::rnorm(n))
  for (K in c(4, 2)) {
    d$y <- cut(d$x1, seq(0, 1, length.out = K + 1), labels = LETTERS[seq_len(K)],
               include.lowest = TRUE)
    a <- attack_fit("gbm", d, d)
    P <- a$point
    expect_equal(unname(rowSums(P)), rep(1, n), tolerance = 1e-6)
    cls <- colnames(P)[max.col(P, "first")]
    expect_gt(mean(cls == d$y), 0.9)

    raw <- stats::predict(a$fit, xgboost::xgb.DMatrix(as.matrix(d[, qi])))
    direct <- if (K == 2) {
      levels(d$y)[1 + (raw > 0.5)]
    } else {
      levels(d$y)[max.col(matrix(raw, nrow = n), "first")]
    }
    expect_identical(cls, direct)
  }
})

test_that("a class absent from the synthetic data has probability 0", {
  skip_if_not_installed("rpart")
  o <- sim_categorical(120)
  s <- to_classes(noisy_copy(sim_numeric(120)), ref = sim_numeric(120))
  s <- s[s$y != "c", ]
  for (m in c("rf", "cart", if (requireNamespace("xgboost", quietly = TRUE)) "gbm")) {
    a <- attack_fit(m, o, s)
    expect_identical(colnames(a$point), c("a", "b", "c"), info = m)
    expect_lt(max(a$point[, "c"]), 0.05)
  }
})

test_that("every learner declares targets, dists and requires", {
  for (m in c("rf", "lm", "cart", "gbm", "logit")) {
    l <- .rapid_learner(m)
    expect_s3_class(l, c(paste0("rapid_learner_", m), "rapid_learner"), exact = TRUE)
    expect_identical(l$name, m)
    expect_type(l$args, "list")
    expect_type(l$requires, "character")
    expect_true(length(l$targets) >= 1 && all(l$targets %in% c("numeric", "categorical")))
    expect_true(all(l$dists %in% c("auto", "grid", "exact")), info = m)
    expect_identical("numeric" %in% l$targets, length(l$dists) > 0, info = m)
    expect_identical(.rapid_learner(l), l)
  }
  expect_error(.rapid_learner("svm"), class = "rapid_error_model")
})

test_that("unsupported combinations are classed errors", {
  err <- tryCatch(.learner_check(.rapid_learner("lm"), TRUE, "auto"), error = function(e) e)
  expect_s3_class(err, "rapid_error_model")
  for (m in c("rf", "cart", "gbm", "logit")) expect_match(conditionMessage(err), m)

  err <- tryCatch(.learner_check(.rapid_learner("logit"), FALSE, "auto"), error = function(e) e)
  expect_s3_class(err, "rapid_error_model")
  for (m in c("rf", "lm", "cart", "gbm")) expect_match(conditionMessage(err), m)

  o <- sim_categorical(90)
  err <- tryCatch(attack_fit("logit", o, o), error = function(e) e)
  expect_s3_class(err, "rapid_error_model")
  for (m in c("rf", "cart", "gbm")) expect_match(conditionMessage(err), m)

  expect_error(.learner_check(.rapid_learner("rf"), FALSE, "exact"),
               class = "rapid_error_model")
  expect_silent(.learner_check(.rapid_learner("rf"), FALSE, "grid"))
  expect_silent(.learner_check(.rapid_learner("rf"), TRUE, "auto"))
})

test_that("grid distributions are sorted n x K matrices without NA", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(100)
  s <- noisy_copy(o)
  for (m in c("rf", "lm", "cart", if (requireNamespace("xgboost", quietly = TRUE)) "gbm")) {
    a <- attack_fit(m, o, s, num_dist = "grid", args = if (m == "gbm") list(nrounds = 20))
    pd <- pdist_of(a, "grid", n_quantiles = 20)
    expect_s3_class(pd, c("rapid_pdist_grid", "rapid_pdist"), exact = TRUE)
    expect_identical(dim(pd$D), c(100L, 20L), info = m)
    expect_false(anyNA(pd$D), info = m)
    expect_true(all(pd$D[, -1] >= pd$D[, -20]), info = m)
  }
  a <- attack_fit("lm", o, s)
  expect_error(pdist_of(a, "grid", n_quantiles = 1), class = "rapid_error_input")
})

test_that("grid distributions for lm and cart reproduce Oscar's", {
  skip_if_not_installed("rpart")
  fx <- rapid_fixture()
  lm_a <- attack_fit("lm", fx$num$o, fx$num$s, num_dist = "grid", num_scale = "identity")
  expect_equal(pdist_of(lm_a, "grid")$D, fx$D$lm, tolerance = 1e-8)
  expect_equal(unname(lm_a$point), fx$point$lm, tolerance = 1e-8)

  cart_a <- attack_fit("cart", fx$num$o, fx$num$s, num_dist = "grid",
                       args = list(prune = FALSE, cp = 0.01, minsplit = 20,
                                   minbucket = 7))
  expect_equal(pdist_of(cart_a, "grid")$D, fx$D$cart, tolerance = 1e-8)
  expect_equal(unname(cart_a$point), fx$point$cart, tolerance = 1e-8)
})

test_that("gbm point predictions and grid are close to Oscar's", {
  skip_on_cran()
  skip_if_not_installed("xgboost")
  fx <- rapid_fixture()
  a <- attack_fit("gbm", fx$num$o, fx$num$s, num_dist = "grid")
  expect_equal(unname(a$point), fx$point$gbm, tolerance = 1e-4)
  expect_equal(pdist_of(a, "grid")$D, fx$D$gbm, tolerance = 1e-4)
})

test_that("exact lm mass equals the normal CDF difference", {
  o <- sim_numeric(100)
  s <- noisy_copy(o)
  a <- attack_fit("lm", o, s, num_scale = "identity")
  pd <- pdist_of(a)
  expect_s3_class(pd, c("rapid_pdist_normal", "rapid_pdist"), exact = TRUE)
  p <- stats::predict(a$fit, newdata = a$X, se.fit = TRUE)
  sd <- sqrt(p$se.fit^2 + summary(a$fit)$sigma^2)
  lo <- o$y - 4
  hi <- o$y + 5
  expect_equal(.pdist_mass(pd, lo, hi),
               unname(stats::pnorm(hi, p$fit, sd) - stats::pnorm(lo, p$fit, sd)))
  expect_equal(unname(a$point), unname(p$fit))
})

test_that("exact log-scale lm mass equals the plnorm difference", {
  o <- sim_lognormal(150, 0.9, 1)
  s <- sim_lognormal(150, 0.9, 2)
  a <- attack_fit("lm", o, s, num_scale = "log")
  pd <- pdist_of(a)
  expect_s3_class(pd, c("rapid_pdist_lognormal", "rapid_pdist"), exact = TRUE)
  ref <- stats::lm(log(y) ~ x1 + x2, data = s)
  p <- stats::predict(ref, newdata = o, se.fit = TRUE)
  sdlog <- sqrt(p$se.fit^2 + summary(ref)$sigma^2)
  lo <- o$y * 0.9
  hi <- o$y * 1.1
  expect_equal(.pdist_mass(pd, lo, hi),
               unname(stats::plnorm(hi, p$fit, sdlog) - stats::plnorm(lo, p$fit, sdlog)))
  # The point prediction is the predictive median
  expect_equal(unname(a$point), unname(exp(p$fit)))
  # Bands that reach zero or lie below it
  expect_equal(.pdist_mass(pd, rep(-5, 150), hi), unname(stats::plnorm(hi, p$fit, sdlog)))
  expect_equal(.pdist_mass(pd, rep(-5, 150), rep(0, 150)), rep(0, 150))
})

test_that("the log-scale grid is exp() of the normal grid", {
  o <- sim_lognormal(80, 0.9, 1)
  s <- sim_lognormal(80, 0.9, 2)
  a <- attack_fit("lm", o, s, num_scale = "log", num_dist = "grid")
  ref <- stats::lm(log(y) ~ x1 + x2, data = s)
  p <- stats::predict(ref, newdata = o, se.fit = TRUE)
  sdlog <- sqrt(p$se.fit^2 + summary(ref)$sigma^2)
  probs <- (seq_len(50) - 0.5) / 50
  expect_equal(pdist_of(a, "grid", 50)$D,
               unname(exp(outer(as.numeric(p$fit), rep(1, 50)) +
                            outer(sdlog, stats::qnorm(probs)))))
})

test_that("auto takes the log scale for a lognormal target and the identity for a normal one", {
  ln <- sim_lognormal(500, 0.9, 3)
  nm <- ln
  nm$y <- 50 + 10 * ln$s
  scale_of <- function(d, metric, num_scale = "auto") {
    fit <- attack_fit("lm", d, d, num_scale = num_scale, metric = metric)$fit
    c(attr(fit, "rapid_scale"), attr(fit, "rapid_scale_auto"))
  }
  expect_identical(scale_of(ln, "symmetric"), c("log", "TRUE"))
  expect_identical(scale_of(ln, "stabilised_relative"), c("log", "TRUE"))
  expect_identical(scale_of(nm, "symmetric"), c("identity", "TRUE"))
  expect_identical(scale_of(ln, "absolute"), c("identity", "TRUE"))
  expect_identical(scale_of(nm, "absolute"), c("identity", "TRUE"))
  expect_identical(scale_of(nm, "symmetric", "log"), c("log", "FALSE"))
  expect_identical(scale_of(ln, "symmetric", "identity"), c("identity", "FALSE"))
})

test_that("zeros and negative values: auto falls back, log is refused", {
  ln <- sim_lognormal(200, 0.9, 4)
  ln$y[1:3] <- c(0, -2, 0)
  fit <- attack_fit("lm", ln, ln, num_scale = "auto")$fit
  expect_identical(attr(fit, "rapid_scale"), "identity")
  err <- tryCatch(attack_fit("lm", ln, ln, num_scale = "log"), error = function(e) e)
  expect_s3_class(err, "rapid_error_model")
  expect_match(conditionMessage(err), "3")
})

test_that("the log scale recovers what the raw scale misses", {
  o <- sim_lognormal(500, 0.99, 5)
  s <- sim_lognormal(500, 0.99, 6)
  rule <- .rapid_rule("band", 0.3, "symmetric", 0.1, 0.01)
  flagged <- function(num_scale) {
    a <- attack_fit("lm", o, s, num_scale = num_scale)
    sc <- .rule_score(rule, list(point = a$point, dist = pdist_of(a)), o$y, o$y)
    mean(.rule_flag(rule, sc$r))
  }
  expect_gt(flagged("log"), 0.5)
  expect_lt(flagged("identity"), 0.15)
  expect_equal(flagged("auto"), flagged("log"))
})

test_that("cart is pruned by cross-validation by default", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(200)
  s <- noisy_copy(o)
  cases <- list(numeric = list(o = o, s = s, method = "anova"),
                factor = list(o = to_classes(o), s = to_classes(s, ref = o), method = "class"))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    prep <- .rapid_prepare_data(cs$o, cs$s, qi, "y", "constant", 0)

    set.seed(11)
    fit <- attack_fit("cart", cs$o, cs$s)$fit
    set.seed(11)
    grown <- rpart::rpart(f_num, data = prep$Y, method = cs$method, cp = 0, minbucket = 5,
                          xval = 10)
    best <- grown$cptable[which.min(grown$cptable[, "xerror"]), "CP"]
    ref <- rpart::prune(grown, cp = best)
    expect_equal(fit$frame, ref$frame, info = nm)
    expect_identical(fit$where, ref$where, info = nm)
    expect_true(attr(fit, "rapid_pruned"), info = nm)

    # rpart's own defaults; minbucket = 7 alone would give minsplit = 21
    plain <- attack_fit("cart", cs$o, cs$s,
                        args = list(prune = FALSE, cp = 0.01, minsplit = 20,
                                    minbucket = 7))$fit
    own <- rpart::rpart(f_num, data = prep$Y, method = cs$method)
    expect_equal(plain$frame, own$frame, info = nm)
    expect_false(attr(plain, "rapid_pruned"), info = nm)
  }
})

test_that("user arguments override the cart defaults one by one", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(150)
  fit <- attack_fit("cart", o, noisy_copy(o), args = list(minbucket = 20, prune = FALSE))$fit
  expect_equal(fit$control$minbucket, 20)
  expect_equal(fit$control$cp, 0)
  expect_equal(fit$control$xval, 10)
  ctrl <- rpart::rpart.control(cp = 0.02, xval = 5)
  fit2 <- attack_fit("cart", o, noisy_copy(o), args = list(control = ctrl))$fit
  expect_equal(fit2$control$cp, 0.02)
  expect_equal(fit2$control$minbucket, ctrl$minbucket)
})

test_that("a target without signal gives a one-leaf tree and RAPID 0", {
  skip_if_not_installed("rpart")
  rule <- .rapid_rule("band", 0.3, "symmetric", 0.1, 0.01)
  rapid_of <- function(o, s) {
    a <- attack_fit("cart", o, s)
    sc <- .rule_score(rule, list(point = a$point, dist = pdist_of(a)), o$y, o$y)
    list(leaves = sum(a$fit$frame$var == "<leaf>"), sc = sc,
         rapid = mean(.rule_flag(rule, sc$r)))
  }
  # Same target values, link to the quasi-identifiers broken
  o <- sim_numeric(200, seed = 21)
  s <- o
  set.seed(22)
  s$y <- sample(s$y)
  set.seed(23)
  res <- rapid_of(o, s)
  expect_equal(res$leaves, 1)
  expect_equal(res$sc$g, res$sc$b)
  expect_equal(res$rapid, 0)

  # No split possible at all: the pruning table has a single row
  o2 <- o
  o2$x1 <- 1
  o2$x2 <- 2
  res2 <- rapid_of(o2, o2)
  expect_equal(res2$leaves, 1)
  expect_equal(res2$rapid, 0)
})

test_that("prune = TRUE with xval = 0 is a classed error", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(100)
  expect_error(attack_fit("cart", o, o, args = list(xval = 0)), class = "rapid_error_model")
  expect_s3_class(attack_fit("cart", o, o, args = list(xval = 0, prune = FALSE))$fit, "rpart")
  expect_error(.rapid_learner("cart", list(prune = "yes")), class = "rapid_error_input")
})

test_that("prune and num_scale with another learner warn", {
  expect_warning(l <- .rapid_learner("rf", list(prune = FALSE, num.trees = 20)),
                 class = "rapid_warning_ignored_argument")
  expect_identical(l$args, list(num.trees = 20))
  expect_warning(.learner_check(.rapid_learner("rf"), FALSE, "auto", num_scale = "log"),
                 class = "rapid_warning_ignored_argument")
  expect_silent(.learner_check(.rapid_learner("rf"), FALSE, "auto", num_scale = "auto"))
  expect_silent(.learner_check(.rapid_learner("lm"), FALSE, "auto", num_scale = "log"))
})

test_that("every learner describes itself in one line", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(120)
  s <- noisy_copy(o)
  o2 <- to_classes(o, probs = 0.5)
  s2 <- to_classes(s, ref = o, probs = 0.5)
  rf <- attack_fit("rf", o, s, args = list(num.trees = 40))
  expect_identical(.learner_describe(rf$learner, rf$fit), "rf: 40 trees, min.node.size 5")
  cart <- attack_fit("cart", o, s)
  leaves <- sum(cart$fit$frame$var == "<leaf>")
  expect_match(.learner_describe(cart$learner, cart$fit),
               paste0("^cart: ", leaves, " lea(f|ves), pruned by cross-validation$"))
  cart0 <- attack_fit("cart", o, s, args = list(prune = FALSE))
  expect_match(.learner_describe(cart0$learner, cart0$fit), "not pruned$")
  lm_auto <- attack_fit("lm", o, s)
  expect_identical(.learner_describe(lm_auto$learner, lm_auto$fit),
                   "lm: identity scale (chosen automatically)")
  lm_log <- attack_fit("lm", o, s, num_scale = "log")
  expect_identical(.learner_describe(lm_log$learner, lm_log$fit), "lm: log scale")
  logit <- attack_fit("logit", o2, s2)
  expect_identical(.learner_describe(logit$learner, logit$fit), "logit")
  skip_if_not_installed("xgboost")
  gbm <- attack_fit("gbm", o, s, args = list(nrounds = 15))
  expect_identical(.learner_describe(gbm$learner, gbm$fit), "gbm: 15 rounds")
})

test_that("exact cart mass equals the share of leaf values in the band", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(150)
  s <- noisy_copy(o)
  a <- attack_fit("cart", o, s, args = list(prune = FALSE, cp = 0.01))
  pd <- pdist_of(a)
  expect_s3_class(pd, c("rapid_pdist_leaf", "rapid_pdist"), exact = TRUE)
  lo <- o$y - 3
  hi <- o$y + 3
  leaf_fit <- a$fit
  leaf_fit$frame$yval <- seq_len(nrow(leaf_fit$frame))
  leaf <- stats::predict(leaf_fit, newdata = a$X)
  brute <- vapply(seq_len(nrow(o)), function(i) {
    v <- a$fit$y[a$fit$where == leaf[i]]
    mean(v > lo[i] & v < hi[i])
  }, numeric(1))
  expect_equal(.pdist_mass(pd, lo, hi), brute)
  expect_true(length(unique(leaf)) > 1)
})

test_that("grid mass converges to the exact mass", {
  o <- sim_numeric(100)
  s <- noisy_copy(o)
  a <- attack_fit("lm", o, s, num_scale = "identity")
  band <- .rapid_band_bounds(o$y, "symmetric", 0.1, 0.01)
  exact <- .pdist_mass(pdist_of(a), band$lo, band$hi)
  grid <- .pdist_mass(pdist_of(a, "grid", n_quantiles = 5000), band$lo, band$hi)
  expect_lt(max(abs(grid - exact)), 1e-3)
})

test_that("user arguments reach every learner", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(100)
  s <- noisy_copy(o)
  expect_equal(attack_fit("rf", o, s, args = list(num.trees = 50))$fit$num.trees, 50)
  expect_equal(attack_fit("cart", o, s, args = list(maxdepth = 2))$fit$control$maxdepth, 2)
  wt <- attack_fit("lm", o, s, args = list(weights = rep(c(1, 2), 50)))$fit
  expect_equal(unname(wt$weights), rep(c(1, 2), 50))
  skip_if_not_installed("xgboost")
  gbm <- attack_fit("gbm", o, s, args = list(max_depth = 2, nrounds = 10))$fit
  cfg <- xgboost::xgb.config(gbm)
  expect_identical(cfg$learner$gradient_booster$tree_train_param$max_depth, "2")
  expect_equal(attr(gbm, "rapid_nrounds"), 10)
})

test_that("variable importance is a named vector for every learner", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(120)
  s <- noisy_copy(o)
  for (m in c("rf", "lm", "cart", if (requireNamespace("xgboost", quietly = TRUE)) "gbm")) {
    a <- attack_fit(m, o, s, args = if (m == "cart") list(prune = FALSE))
    imp <- .learner_importance(a$learner, a$fit)
    expect_true(is.numeric(imp) && all(names(imp) %in% qi) && length(imp) >= 1, info = m)
  }
})

test_that("a missing backend package is a classed error", {
  local_mocked_bindings(.rapid_has_package = function(pkg) FALSE)
  for (m in c("rf", "cart", "gbm")) {
    err <- tryCatch(.learner_check(.rapid_learner(m), FALSE, "auto"), error = function(e) e)
    expect_s3_class(err, "rapid_error_dependency")
    expect_match(conditionMessage(err), "install.packages", fixed = TRUE)
  }
  expect_silent(.learner_check(.rapid_learner("lm"), FALSE, "auto"))
})

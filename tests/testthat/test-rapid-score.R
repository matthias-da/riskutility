# Scoring core of RAPID: band bounds, baseline, gain, rule classes

metrics <- c("symmetric", "stabilised_relative", "absolute")
eps_of <- function(metric) if (metric == "absolute") 2 else 0.1
y_grid <- c(-120, -3, -0.004, 0, 0.004, 3, 50, 50000)

test_that("band bounds agree with the error metric for all metrics and signs", {
  for (metric in metrics) {
    eps <- eps_of(metric)
    b <- .rapid_band_bounds(y_grid, metric, eps, 0.01)
    expect_named(b, c("lo", "hi"))
    h <- 1e-6 * (b$hi - b$lo)
    expect_true(all(brute_error(y_grid, b$lo + h, metric) < eps), info = metric)
    expect_true(all(brute_error(y_grid, b$hi - h, metric) < eps), info = metric)
    expect_true(all(brute_error(y_grid, b$lo - h, metric) >= eps), info = metric)
    expect_true(all(brute_error(y_grid, b$hi + h, metric) >= eps), info = metric)
  }
})

test_that("paper example: 50000 with 10 % symmetric tolerance", {
  b <- .rapid_band_bounds(50000, "symmetric", 0.1, 0)
  expect_equal(round(b$lo), 45238)
  expect_equal(round(b$hi), 55263)
})

test_that("a symmetric tolerance of 2 or more gives an unbounded band", {
  b <- .rapid_band_bounds(c(-3, 0, 7), "symmetric", 2, 0.01)
  expect_true(all(b$lo == -Inf) && all(b$hi == Inf))
})

test_that("categorical toy example of the paper", {
  ref <- factor(rep(c("healthy", "diabetic", "hypertensive"), c(60, 20, 20)))
  y <- factor(rep("healthy", 3), levels = levels(ref))
  g <- c(0.70, 0.85, 0.55)
  P <- cbind(healthy = g, diabetic = (1 - g) / 2, hypertensive = (1 - g) / 2)
  rule <- .rapid_rule("RCS_marginal", 0.3)
  sc <- .rule_score(rule, P, y, ref)
  expect_equal(sc$g, g)
  expect_equal(sc$b, rep(0.6, 3))
  expect_equal(sc$r, c(0.25, 0.625, -0.125))
  expect_identical(.rule_flag(rule, sc$r), c(FALSE, TRUE, FALSE))
})

test_that("continuous toy example of the paper", {
  r <- .rapid_gain(c(0.36, 0.60, 0.10), c(0.30, 0.05, 0.20))
  expect_equal(r, c(0.06 / 0.70, 0.55 / 0.95, -0.125))
  rule <- .rapid_rule("band", 0.3, "symmetric", 0.1, 0.01)
  expect_equal(sum(.rule_flag(rule, r)), 1)
})

test_that("g, b, r equal a brute-force count", {
  set.seed(42)
  n <- 40
  y <- stats::rlnorm(n, 3, 1)
  D <- matrix(y * exp(stats::rnorm(n * 50, 0, 0.15)), nrow = n)
  for (metric in metrics) {
    eps <- eps_of(metric)
    rule <- .rapid_rule("band", 0.3, metric, eps, 0.01)
    sc <- .rule_score(rule, list(point = rowMeans(D), dist = .pdist_grid(D)), y, y)
    g <- vapply(seq_len(n), function(i) mean(brute_error(y[i], D[i, ], metric) < eps), 0)
    b <- vapply(seq_len(n), function(i) mean(brute_error(y[i], y, metric) < eps), 0)
    expect_equal(sc$g, g, info = metric)
    expect_equal(sc$b, b, info = metric)
    expect_equal(sc$r, (g - b) / (1 - b), info = metric)
    expect_named(sc, c("predicted_value", "band_lo", "band_hi", "g", "b", "r"))
  }
})

test_that("the baseline counts reference values strictly inside the band", {
  expect_equal(.rapid_baseline_num(lo = c(1, 0), hi = c(3, 10), ref_y = c(1, 2, 3, 4)),
               c(1 / 4, 1))
})

test_that("constant target gives b = 1 and r = 0", {
  y <- rep(5, 20)
  D <- matrix(seq(1, 9, length.out = 200), nrow = 20)
  for (metric in metrics) {
    rule <- .rapid_rule("band", 0.3, metric, eps_of(metric), 0.01)
    sc <- .rule_score(rule, list(point = rowMeans(D), dist = .pdist_grid(D)), y, y)
    expect_equal(sc$b, rep(1, 20), info = metric)
    expect_equal(sc$r, rep(0, 20), info = metric)
    expect_false(anyNA(sc), info = metric)
    expect_false(any(.rule_flag(rule, sc$r)), info = metric)
  }
})

test_that("a single class gives b = 1 and r = 0", {
  y <- factor(rep("a", 5))
  P <- matrix(1, nrow = 5, dimnames = list(NULL, "a"))
  sc <- .rule_score(.rapid_rule("RCS_marginal", 0.3), P, y, y)
  expect_equal(sc$r, rep(0, 5))
})

test_that("bands crossing zero are ordered", {
  for (metric in metrics) {
    b <- .rapid_band_bounds(c(-0.004, 0, 0.004), metric, eps_of(metric), 0.01)
    expect_true(all(b$lo < b$hi), info = metric)
    expect_true(all(b$lo < c(-0.004, 0, 0.004) & c(-0.004, 0, 0.004) < b$hi), info = metric)
  }
})

test_that("point rule flags error below epsilon", {
  rule <- .rapid_rule("point", 0.1, "symmetric", delta = 0.01)
  expect_identical(rule$epsilon, 0.1)
  y <- rep(100, 4)
  pred <- list(point = c(95, 109, 111, 90), dist = NULL)
  sc <- .rule_score(rule, pred, y, y)
  expect_named(sc, c("predicted_value", "error", "b"))
  expect_equal(sc$error, 2 * abs(y - pred$point) / (abs(y) + abs(pred$point) + 0.02))
  expect_identical(.rule_flag(rule, sc$error), c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(sc$b, rep(1, 4))
})

test_that("RCS_conditional and NCE match hand values", {
  y <- factor(c("a", "a", "b", "b"))
  P <- rbind(c(0.8, 0.2), c(0.6, 0.4), c(0.3, 0.7), c(0.5, 0.5))
  colnames(P) <- c("a", "b")
  g <- c(0.8, 0.6, 0.7, 0.5)

  rc <- .rapid_rule("RCS_conditional", 1)
  sc <- .rule_score(rc, P, y, y)
  expect_named(sc, c("predicted_class", "g", "b", "ratio"))
  expect_equal(sc$g, g)
  expect_equal(sc$b, c(0.7, 0.7, 0.6, 0.6))
  expect_equal(sc$ratio, c(8 / 7, 6 / 7, 7 / 6, 5 / 6))
  expect_identical(.rule_flag(rc, sc$ratio), c(TRUE, FALSE, TRUE, FALSE))
  expect_identical(sc$predicted_class, c("a", "a", "b", "a"))

  rn <- .rapid_rule("NCE", 0.5)
  sn <- .rule_score(rn, P, y, y)
  expect_named(sn, c("predicted_class", "g", "nce"))
  expect_equal(sn$nce, 1 + log(g) / log(2))
  expect_identical(.rule_flag(rn, sn$nce), c(TRUE, FALSE, FALSE, FALSE))
})

test_that("every rule has a score method, a score_var column and a format line", {
  y_num <- c(10, 20, 30, 40)
  pred_num <- list(point = c(11, 19, 35, 20),
                   dist = .pdist_grid(matrix(rep(c(11, 19, 35, 20), 5), nrow = 4)))
  y_cat <- factor(c("a", "a", "b", "b"))
  P <- cbind(a = c(0.8, 0.6, 0.3, 0.5), b = c(0.2, 0.4, 0.7, 0.5))
  rules <- list(
    band = .rapid_rule("band", 0.3, "symmetric", 0.1, 0.01),
    point = .rapid_rule("point", 0.1, "symmetric", delta = 0.01),
    RCS_marginal = .rapid_rule("RCS_marginal", 0.3),
    RCS_conditional = .rapid_rule("RCS_conditional", 1),
    NCE = .rapid_rule("NCE", 0.5)
  )
  for (m in names(rules)) {
    rule <- rules[[m]]
    numeric_rule <- m %in% c("band", "point")
    sc <- .rule_score(rule, if (numeric_rule) pred_num else P,
                      if (numeric_rule) y_num else y_cat,
                      if (numeric_rule) y_num else y_cat)
    expect_s3_class(rule, c(paste0("rapid_rule_", m), "rapid_rule"), exact = TRUE)
    expect_identical(rule$method, m)
    expect_s3_class(sc, "data.frame")
    expect_equal(nrow(sc), 4)
    expect_true(rule$score_var %in% names(sc), info = m)
    expect_type(.rule_flag(rule, sc[[rule$score_var]]), "logical")
    expect_length(format(rule), 1)
    expect_type(format(rule), "character")
  }
  expect_identical(vapply(rules, function(r) r$direction, ""),
                   c(band = "above", point = "below", RCS_marginal = "above",
                     RCS_conditional = "above", NCE = "above"))
  expect_identical(vapply(rules, function(r) r$gain, NA),
                   c(band = TRUE, point = FALSE, RCS_marginal = TRUE,
                     RCS_conditional = FALSE, NCE = FALSE))
  expect_identical(format(rules$band), "band: symmetric error < 10 %, gain > 0.3")
})

test_that("a new threshold replaces the old one, and the tolerance of the point rule", {
  band <- .rule_set_threshold(.rapid_rule("band", 0.3, "symmetric", 0.1, 0.01), 0.5)
  expect_identical(band$threshold, 0.5)
  expect_identical(band$epsilon, 0.1)
  point <- .rule_set_threshold(.rapid_rule("point", 0.1, "symmetric", delta = 0.01), 0.2)
  expect_identical(point$threshold, 0.2)
  expect_identical(point$epsilon, 0.2)
})

test_that("an unknown method is a classed error", {
  expect_error(.rapid_rule("foo", 0.3), class = "rapid_error_input")
  expect_error(.rapid_band_bounds(1, "foo", 0.1, 0.01), class = "rapid_error_input")
})

test_that("predictions that do not match the records are a classed error", {
  rule <- .rapid_rule("band", 0.3, "symmetric", 0.1, 0.01)
  pred <- list(point = 1:3, dist = .pdist_grid(matrix(1, nrow = 3, ncol = 2)))
  expect_error(.rule_score(rule, pred, c(1, 2), c(1, 2)), class = "rapid_error_input")
  P <- cbind(a = c(0.5, 0.5), b = c(0.5, 0.5))
  expect_error(.rule_score(.rapid_rule("RCS_marginal", 0.3), P,
                           factor(c("a", "c")), factor(c("a", "c"))),
               class = "rapid_error_input")
})

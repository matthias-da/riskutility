# Shared helpers for the RAPID tests.
#
# The same simulation design is used by fixtures/make-rapid-ot-fixture.R, which
# produces the reference values from Oscar Thees' RAPID package.

# Original data: two quasi-identifiers and a numeric target whose correlation
# with the linear signal is `signal`.
sim_numeric <- function(n = 150, seed = 1, signal = 0.8) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  y <- 75 + 10 * (signal * (x1 + x2) / sqrt(2) + sqrt(1 - signal^2) * stats::rnorm(n))
  data.frame(x1 = x1, x2 = x2, y = y)
}

# Synthetic data: a row-shuffled copy with normal noise added to a numeric
# target. `shuffle_target = TRUE` breaks the link between target and
# quasi-identifiers (an uninformative release).
noisy_copy <- function(d, sd = 2, target = "y", seed = NULL, shuffle_target = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  s <- d[sample(nrow(d)), , drop = FALSE]
  rownames(s) <- NULL
  if (is.numeric(s[[target]]) && sd > 0) {
    s[[target]] <- s[[target]] + stats::rnorm(nrow(s), 0, sd)
  }
  if (shuffle_target) s[[target]] <- s[[target]][sample(nrow(s))]
  s
}

# Turns the numeric target into the classes a, b, c, cut at the 50 % and 80 %
# quantiles of the reference data, so that both data sets share the cut points.
to_classes <- function(d, ref = d, target = "y", probs = c(0.5, 0.8)) {
  breaks <- c(-Inf, stats::quantile(ref[[target]], probs, names = FALSE), Inf)
  d[[target]] <- cut(d[[target]], breaks, labels = letters[seq_len(length(probs) + 1)])
  d
}

sim_categorical <- function(n = 150, seed = 1, signal = 0.8) {
  to_classes(sim_numeric(n, seed, signal))
}

# Error of the value(s) `v` against the true value `y`, written out directly
# from the definitions of the three metrics.
brute_error <- function(y, v, metric, delta = 0.01) {
  switch(metric,
    symmetric = 2 * abs(y - v) / (abs(y) + abs(v) + 2 * delta),
    stabilised_relative = abs(y - v) / (abs(y) + delta),
    absolute = abs(y - v)
  )
}

rapid_fixture <- function() {
  readRDS(testthat::test_path("fixtures", "rapid_ot_fixture.rds"))
}

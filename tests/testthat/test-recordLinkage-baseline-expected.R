# compute_baseline (no-perturbation reference) and PRAM expected_risk
# (perturbation-aware expected risk over the transition distribution).

library(testthat)

test_that("compute_baseline attaches a baseline with risk_reduction", {
  set.seed(1)
  n <- 40
  x <- data.frame(
    age = round(rnorm(n, 45, 10)),
    sex = factor(sample(c("f", "m"), n, TRUE)),
    region = factor(sample(paste0("R", 1:4), n, TRUE))
  )
  x_anon <- x
  for (rg in levels(x$region)) {            # swap age within region
    idx <- which(x$region == rg)
    if (length(idx) > 1) x_anon$age[idx] <- sample(x$age[idx])
  }
  res <- recordLinkage(x, x_anon, key = c("age", "sex", "region"),
                       truth = "row", compute_baseline = TRUE)

  expect_false(is.null(res$baseline))
  expect_length(res$baseline$per_record_risk, nrow(x))
  # Self-linkage (no perturbation) is at least as identifying as the perturbed.
  expect_gte(res$baseline$mean_risk, res$overall$mean_risk - 1e-9)
  expect_equal(res$baseline$risk_reduction,
               res$baseline$mean_risk - res$overall$mean_risk)
})

test_that("compute_baseline short-circuits when X == x_anon", {
  set.seed(2)
  x <- data.frame(a = round(rnorm(20, 10, 3)),
                  b = factor(sample(c("x", "y"), 20, TRUE)))
  res <- recordLinkage(x, x, key = c("a", "b"), truth = "row",
                       compute_baseline = TRUE)
  expect_equal(res$baseline$mean_risk, res$overall$mean_risk)
  expect_equal(res$baseline$risk_reduction, 0)
})

test_that("PRAM expected_risk is bounded and present", {
  lv <- c("a", "b", "c")
  # diagonal-dominant transition matrix
  tm <- matrix(0.1, 3, 3, dimnames = list(lv, lv)); diag(tm) <- 0.8
  pm <- list(k1 = tm, k2 = tm)
  set.seed(3)
  n <- 12
  x <- data.frame(k1 = factor(sample(lv, n, TRUE), levels = lv),
                  k2 = factor(sample(lv, n, TRUE), levels = lv))
  x_anon <- x
  res <- recordLinkage(x, x_anon, key = c("k1", "k2"), method = "pram",
                       pram_matrix = pm, truth = "row", expected_risk = TRUE)
  expect_false(is.null(res$pram_info$expected_risk))
  expect_length(res$pram_info$expected_risk, n)
  expect_true(all(res$pram_info$expected_risk >= 0 &
                  res$pram_info$expected_risk <= 1))
  expect_true(is.finite(res$pram_info$expected_mean_risk))
})

test_that("PRAM expected_risk ~ realized risk for a near-deterministic matrix", {
  lv <- c("a", "b", "c")
  tm <- matrix(0.005, 3, 3, dimnames = list(lv, lv)); diag(tm) <- 0.99
  pm <- list(k1 = tm, k2 = tm)
  set.seed(4)
  n <- 10
  x <- data.frame(k1 = factor(sample(lv, n, TRUE), levels = lv),
                  k2 = factor(sample(lv, n, TRUE), levels = lv))
  res <- recordLinkage(x, x, key = c("k1", "k2"), method = "pram",
                       pram_matrix = pm, truth = "row", expected_risk = TRUE)
  # Near-deterministic mechanism: expectation close to the realized risk.
  expect_lt(mean(abs(res$pram_info$expected_risk - res$per_record$risk)), 0.15)
})

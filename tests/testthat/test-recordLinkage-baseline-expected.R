# compute_baseline (no-perturbation reference).

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


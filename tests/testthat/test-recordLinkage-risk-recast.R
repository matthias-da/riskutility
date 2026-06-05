# Regression tests for the RF / embedding risk recast: `risk` is now the
# probability of identifying the TRUE match (consistent with the other 6
# methods), and the former nearest-neighbour similarity is kept in
# `nn_similarity`. A cyclic shift makes every record's true (row) match hold a
# DIFFERENT person's data, so the nearest released record is never the true one
# -> new risk ~ 0 while nn_similarity stays high (old risk would have been ~ 1).

library(testthat)

.shift_data <- function(n = 60, seed = 7) {
  set.seed(seed)
  x <- data.frame(
    age = round(rnorm(n, 45, 12)),
    sex = factor(sample(c("f", "m"), n, TRUE)),
    inc = round(rnorm(n, 50, 15))
  )
  list(x = x, x_anon = x[c(2:n, 1), , drop = FALSE])  # derangement (no fixed pt)
}

test_that("RF: risk uses the true match, not the nearest neighbour", {
  skip_if_not_installed("ranger")
  d <- .shift_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "inc"),
                       method = "rf", truth = "row", n_trees = 200)

  expect_true("nn_similarity" %in% names(res$per_record))
  expect_true(all(res$per_record$nn_similarity >= 0 &
                  res$per_record$nn_similarity <= 1))
  expect_true(all(res$per_record$risk >= 0 & res$per_record$risk <= 1))
  # New semantics: true row never holds the record's own data -> low risk,
  # but each record's data exists elsewhere -> high nearest-neighbour similarity.
  expect_lt(mean(res$per_record$risk), 0.25)
  expect_gt(mean(res$per_record$nn_similarity), mean(res$per_record$risk))
})

test_that("RF: risk == 1 when the true match IS the unique nearest", {
  skip_if_not_installed("ranger")
  set.seed(3)
  n <- 50
  x <- data.frame(age = round(rnorm(n, 45, 12)),
                  sex = factor(sample(c("f", "m"), n, TRUE)),
                  inc = round(rnorm(n, 50, 15)))
  res <- recordLinkage(x, x, key = c("age", "sex", "inc"),  # identity match
                       method = "rf", truth = "row", n_trees = 200)
  # With identity data the true row is (almost always) the nearest -> risk high.
  expect_gt(mean(res$per_record$risk), 0.5)
})

test_that("embedding: risk uses the true match, not the nearest neighbour", {
  skip_if_no_torch()
  d <- .shift_data(n = 80)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "inc"),
                       method = "embedding", truth = "row",
                       emb_epochs = 15, emb_latent_dim = 4)
  expect_true("nn_similarity" %in% names(res$per_record))
  expect_true(all(res$per_record$risk >= 0 & res$per_record$risk <= 1))
  expect_lt(mean(res$per_record$risk), 0.3)
  expect_gt(mean(res$per_record$nn_similarity), mean(res$per_record$risk))
})

test_that("non-rf/embedding methods do not carry an nn_similarity column", {
  d <- .shift_data(n = 30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "inc"),
                       method = "deterministic", truth = "row")
  expect_false("nn_similarity" %in% names(res$per_record))
})

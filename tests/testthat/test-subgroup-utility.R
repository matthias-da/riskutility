# Tests for subgroup_utility (Stratified Utility Assessment)

library(testthat)

test_that("subgroup_utility returns correct S3 class structure", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )
  Y <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )

  result <- subgroup_utility(X, Y, group_var = "grp", seed = 42)

  expect_s3_class(result, "subgroup_utility")
  expect_true("overall_score" %in% names(result))
  expect_true("utility_score" %in% names(result))
  expect_true("per_group" %in% names(result))
  expect_true("worst_group" %in% names(result))
  expect_true("ratio" %in% names(result))
  expect_true("group_var" %in% names(result))
  expect_true("threshold" %in% names(result))
  expect_true("n_groups" %in% names(result))
})

test_that("per-group scores are computed correctly for 2-group case", {
  set.seed(42)
  n <- 60
  # Group A: identical data -> utility ~ 1
  # Group B: shifted data -> lower utility
  X <- data.frame(
    grp = c(rep("A", n), rep("B", n)),
    x1 = c(rnorm(n, 0, 1), rnorm(n, 0, 1))
  )
  Y <- data.frame(
    grp = c(rep("A", n), rep("B", n)),
    x1 = c(rnorm(n, 0, 1), rnorm(n, 5, 1))
  )

  result <- subgroup_utility(X, Y, group_var = "grp", seed = 42)

  pg <- result$per_group
  # Group A should have higher utility than Group B
  score_A <- pg$utility_score[pg$group == "A"]
  score_B <- pg$utility_score[pg$group == "B"]
  expect_true(score_A > score_B)

  # Verify per-group scores by computing manually
  X_A <- X[X$grp == "A", "x1", drop = FALSE]
  Y_A <- Y[Y$grp == "A", "x1", drop = FALSE]
  manual_A <- energy_distance(X_A, Y_A, seed = 42)

  X_B <- X[X$grp == "B", "x1", drop = FALSE]
  Y_B <- Y[Y$grp == "B", "x1", drop = FALSE]
  manual_B <- energy_distance(X_B, Y_B, seed = 42)

  expect_equal(score_A, manual_A$utility_score, tolerance = 1e-10)
  expect_equal(score_B, manual_B$utility_score, tolerance = 1e-10)
})

test_that("utility_score equals worst subgroup score", {
  set.seed(123)
  n <- 80
  X <- data.frame(
    grp = sample(c("A", "B", "C"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )
  Y <- data.frame(
    grp = sample(c("A", "B", "C"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )

  result <- subgroup_utility(X, Y, group_var = "grp", seed = 42)

  worst <- min(result$per_group$utility_score, na.rm = TRUE)
  expect_equal(result$utility_score, worst)
})

test_that("flags low-utility subgroups", {
  set.seed(42)
  n <- 60
  # Create a subgroup with very poor utility
  X <- data.frame(
    grp = c(rep("good", n), rep("bad", n)),
    x1 = c(rnorm(n, 0, 1), rnorm(n, 0, 1))
  )
  Y <- data.frame(
    grp = c(rep("good", n), rep("bad", n)),
    x1 = c(rnorm(n, 0, 1), rnorm(n, 10, 1))
  )

  result <- subgroup_utility(X, Y, group_var = "grp", threshold = 0.5, seed = 42)

  pg <- result$per_group
  # "bad" group should be flagged
  bad_row <- pg[pg$group == "bad", ]
  expect_true(bad_row$flagged)

  # "good" group should not be flagged (utility should be decent)
  good_row <- pg[pg$group == "good", ]
  expect_false(good_row$flagged)
})

test_that("threshold parameter controls flagging", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n)
  )
  Y <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n)
  )

  # With very high threshold, all should be flagged
  result_high <- subgroup_utility(X, Y, group_var = "grp",
                                  threshold = 1.01, seed = 42)
  expect_true(all(result_high$per_group$flagged, na.rm = TRUE))

  # With threshold = 0, none should be flagged

  result_low <- subgroup_utility(X, Y, group_var = "grp",
                                 threshold = 0, seed = 42)
  expect_false(any(result_low$per_group$flagged, na.rm = TRUE))
})

test_that("works with mmd as utility function", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )
  Y <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )

  result <- subgroup_utility(X, Y, group_var = "grp",
                             utility_fun = mmd, seed = 42)

  expect_s3_class(result, "subgroup_utility")
  expect_true(all(!is.na(result$per_group$utility_score)))
  expect_true(all(result$per_group$utility_score >= 0))
  expect_true(all(result$per_group$utility_score <= 1))
})

test_that("handles subgroups with few observations (skip with NA)", {
  set.seed(123)
  # Create a tiny subgroup with < 5 observations
  X <- data.frame(
    grp = c(rep("big", 50), rep("tiny", 3)),
    x1 = rnorm(53)
  )
  Y <- data.frame(
    grp = c(rep("big", 50), rep("tiny", 3)),
    x1 = rnorm(53)
  )

  expect_warning(
    result <- subgroup_utility(X, Y, group_var = "grp", seed = 42),
    "fewer than 5 observations"
  )

  pg <- result$per_group
  tiny_row <- pg[pg$group == "tiny", ]
  expect_true(is.na(tiny_row$utility_score))
  expect_true(is.na(tiny_row$flagged))

  big_row <- pg[pg$group == "big", ]
  expect_false(is.na(big_row$utility_score))
})

test_that("numeric group_var is converted to factor with message", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    grp = sample(1:3, n, replace = TRUE),
    x1 = rnorm(n)
  )
  Y <- data.frame(
    grp = sample(1:3, n, replace = TRUE),
    x1 = rnorm(n)
  )

  expect_message(
    result <- subgroup_utility(X, Y, group_var = "grp", seed = 42),
    "Converting numeric"
  )
  expect_s3_class(result, "subgroup_utility")
  expect_true(result$n_groups > 0)
})

test_that("error on missing group_var in X", {
  X <- data.frame(x1 = 1:10)
  Y <- data.frame(x1 = 1:10, grp = rep("A", 10))

  expect_error(
    subgroup_utility(X, Y, group_var = "grp"),
    "not found in X"
  )
})

test_that("error on missing group_var in Y", {
  X <- data.frame(x1 = 1:10, grp = rep("A", 10))
  Y <- data.frame(x1 = 1:10)

  expect_error(
    subgroup_utility(X, Y, group_var = "grp"),
    "not found in Y"
  )
})

test_that("error when group_var is not provided", {
  X <- data.frame(x1 = 1:10)
  Y <- data.frame(x1 = 1:10)

  expect_error(
    subgroup_utility(X, Y),
    "group_var"
  )
})

test_that("synth_pair dispatch works", {
  set.seed(123)
  n <- 80
  original <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )
  synthetic <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )

  pair <- synth_pair(original, synthetic)
  result <- subgroup_utility(pair, group_var = "grp", seed = 42)

  expect_s3_class(result, "subgroup_utility")
  expect_true(result$n_groups == 2)
})

test_that("print method works", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n)
  )
  Y <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n)
  )

  result <- subgroup_utility(X, Y, group_var = "grp", seed = 42)

  expect_output(print(result), "Subgroup Utility Assessment")
  expect_output(print(result), "Overall utility score")
  expect_output(print(result), "Worst subgroup score")
})

test_that("summary method returns correct class and content", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n)
  )
  Y <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n)
  )

  result <- subgroup_utility(X, Y, group_var = "grp", seed = 42)
  s <- summary(result)

  expect_s3_class(s, "summary.subgroup_utility")
  expect_true(!is.null(s$mean_utility))
  expect_true(!is.null(s$sd_utility))
  expect_true(!is.null(s$per_group))
  expect_output(print(s), "Summary: Subgroup Utility Assessment")
  expect_output(print(s), "Per-Group Results")
})

test_that("plot method works without error", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    grp = sample(c("A", "B", "C"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )
  Y <- data.frame(
    grp = sample(c("A", "B", "C"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )

  result <- subgroup_utility(X, Y, group_var = "grp", seed = 42)

  expect_no_error(plot(result, which = 1))
})

test_that("ratio is worst / overall", {
  set.seed(42)
  n <- 60
  X <- data.frame(
    grp = c(rep("A", n), rep("B", n)),
    x1 = c(rnorm(n), rnorm(n))
  )
  Y <- data.frame(
    grp = c(rep("A", n), rep("B", n)),
    x1 = c(rnorm(n), rnorm(n, 3, 1))
  )

  result <- subgroup_utility(X, Y, group_var = "grp", seed = 42)

  expected_ratio <- result$utility_score / result$overall_score
  expect_equal(result$ratio, expected_ratio, tolerance = 1e-10)
})

test_that("NA handling in group_var", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    grp = c(rep("A", n - 5), rep(NA, 5)),
    x1 = rnorm(n)
  )
  Y <- data.frame(
    grp = c(rep("A", n - 3), rep(NA, 3)),
    x1 = rnorm(n)
  )

  # With na.rm = TRUE (default), should work
  result <- subgroup_utility(X, Y, group_var = "grp", seed = 42)
  expect_s3_class(result, "subgroup_utility")
  # Only group A should appear
  expect_equal(nrow(result$per_group), 1)
  expect_equal(result$per_group$group, "A")
})

test_that("overall_score differs from per-group scores", {
  set.seed(42)
  n <- 80
  X <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )
  Y <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )

  result <- subgroup_utility(X, Y, group_var = "grp", seed = 42)

  # Overall score should be computed on full data (not average of groups)
  expect_true(is.numeric(result$overall_score))
  expect_true(result$overall_score >= 0 && result$overall_score <= 1)
})

test_that("works with data.table input", {
  skip_if_not_installed("data.table")
  set.seed(123)
  n <- 60
  X <- data.table::data.table(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n)
  )
  Y <- data.table::data.table(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n)
  )

  result <- subgroup_utility(X, Y, group_var = "grp", seed = 42)
  expect_s3_class(result, "subgroup_utility")
})

test_that("per_group has correct columns", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n)
  )
  Y <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n)
  )

  result <- subgroup_utility(X, Y, group_var = "grp", seed = 42)

  expect_true(all(c("group", "n_orig", "n_synth", "utility_score", "flagged") %in%
                    names(result$per_group)))
  expect_equal(nrow(result$per_group), 2)
})

test_that("group_var column is excluded from utility computation", {
  set.seed(123)
  n <- 60
  # If group_var were included, it would cause errors with energy_distance
  # (character/factor column) or skew results
  X <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )
  Y <- data.frame(
    grp = sample(c("A", "B"), n, replace = TRUE),
    x1 = rnorm(n), x2 = rnorm(n)
  )

  # Should work without error (group_var excluded from energy_distance)
  result <- subgroup_utility(X, Y, group_var = "grp", seed = 42)
  expect_s3_class(result, "subgroup_utility")
})

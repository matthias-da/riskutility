# Tests for rapid_qi_contribution()

.make_contrib_data <- function(n = 300, seed = 1) {
  set.seed(seed)
  X <- data.frame(age = sample(20:70, n, replace = TRUE),
                  sex = factor(sample(c("F", "M"), n, replace = TRUE)),
                  region = factor(sample(letters[1:4], n, replace = TRUE)))
  X$smoker <- factor(ifelse(X$age < 40 & X$sex == "M", "yes", "no"))
  Y <- X
  Y$age <- Y$age + sample(-2:2, n, replace = TRUE)
  list(X = X, Y = Y)
}

test_that("rapid_qi_contribution returns one row per drop plus the reference", {
  d <- .make_contrib_data()
  res <- rapid_qi_contribution(d$X, d$Y, c("age", "sex", "region"), "smoker",
                               groups = list("age+sex" = c("age", "sex")),
                               model_type = "cart", n_permutations = 9, seed = 1)
  expect_s3_class(res, "rapid_contribution")
  r <- res$results
  expect_equal(nrow(r), 1 + 3 + 1)
  expect_equal(r$dropped[1], "(none)")
  expect_setequal(r$dropped[-1], c("age", "sex", "region", "age+sex"))
  expect_equal(r$delta_rapid[1], 0)
  expect_equal(r$n_qi[r$dropped == "age+sex"], 1)
  expect_true(all(r$rapid >= 0 & r$rapid <= 1))
  expect_equal(r$excess, r$rapid - r$null_mean)
  # contributions sorted in decreasing order after the reference row
  expect_false(is.unsorted(rev(r$delta_excess[-1])))
})

test_that("informative QIs contribute more than noise QIs", {
  d <- .make_contrib_data()
  res <- rapid_qi_contribution(d$X, d$Y, c("age", "sex", "region"), "smoker",
                               model_type = "cart", n_permutations = 19, seed = 1)
  r <- res$results
  expect_gt(r$delta_rapid[r$dropped == "age"], r$delta_rapid[r$dropped == "region"])
  expect_gt(r$delta_rapid[r$dropped == "sex"], r$delta_rapid[r$dropped == "region"])
})

test_that("n_permutations = 0 skips the null", {
  d <- .make_contrib_data()
  res <- rapid_qi_contribution(d$X, d$Y, c("age", "sex", "region"), "smoker",
                               model_type = "cart", n_permutations = 0, seed = 1)
  expect_true(all(is.na(res$results$excess)))
  expect_s3_class(plot(res), "ggplot")
})

test_that("input validation", {
  d <- .make_contrib_data()
  expect_error(rapid_qi_contribution(d$X, d$Y, "age", "smoker"))
  expect_error(rapid_qi_contribution(d$X, d$Y, c("age", "sex"), "smoker",
                                     groups = list(c("age"))), "named list")
  expect_error(rapid_qi_contribution(d$X, d$Y, c("age", "sex"), "smoker",
                                     groups = list(g = c("age", "foo"))), "not in")
  expect_error(rapid_qi_contribution(d$X, d$Y, c("age", "sex"), "smoker",
                                     groups = list(g = c("age", "sex"))), "all quasi")
})

test_that("print, summary and plot methods work", {
  d <- .make_contrib_data()
  res <- rapid_qi_contribution(d$X, d$Y, c("age", "sex", "region"), "smoker",
                               model_type = "cart", n_permutations = 9, seed = 1)
  expect_output(print(res), "risk contribution")
  expect_s3_class(summary(res), "data.frame")
  expect_s3_class(plot(res), "ggplot")
  expect_s3_class(plot(res, measure = "rapid"), "ggplot")
})

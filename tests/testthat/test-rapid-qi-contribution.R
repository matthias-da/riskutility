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

test_that("an argument the contribution sets itself is a classed error", {
  d <- .make_contrib_data(120)
  expect_error(rapid_qi_contribution(d$X, d$Y, key_vars = c("age", "sex", "region"),
                                     target_var = "smoker", model_type = "rf", n_perm = 5,
                                     n_permutations = 5),
               regexp = "n_permutations", class = "rapid_error_input")
})

test_that("rapid_qi_contribution returns one row per drop plus the reference", {
  d <- .make_contrib_data()
  res <- rapid_qi_contribution(d$X, d$Y, key_vars = c("age", "sex", "region"), target_var = "smoker",
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
  expect_equal(r$excess, r$rapid - r$chance)
  expect_identical(names(r), c("dropped", "n_qi", "rapid", "chance", "excess",
                               "delta_rapid", "delta_excess"))
  # contributions sorted in decreasing order after the reference row
  expect_false(is.unsorted(rev(r$delta_excess[-1])))
})

test_that("informative QIs contribute more than noise QIs", {
  d <- .make_contrib_data()
  res <- rapid_qi_contribution(d$X, d$Y, key_vars = c("age", "sex", "region"), target_var = "smoker",
                               model_type = "cart", n_permutations = 19, seed = 1)
  r <- res$results
  expect_gt(r$delta_rapid[r$dropped == "age"], r$delta_rapid[r$dropped == "region"])
  expect_gt(r$delta_rapid[r$dropped == "sex"], r$delta_rapid[r$dropped == "region"])
})

test_that("n_permutations = 0 skips the null", {
  d <- .make_contrib_data()
  res <- rapid_qi_contribution(d$X, d$Y, key_vars = c("age", "sex", "region"), target_var = "smoker",
                               model_type = "cart", n_permutations = 0, seed = 1)
  expect_true(all(is.na(res$results$excess)))
  expect_s3_class(plot(res), "ggplot")
})

test_that("input validation", {
  d <- .make_contrib_data()
  expect_error(rapid_qi_contribution(d$X, d$Y, key_vars = "age", target_var = "smoker"))
  expect_error(rapid_qi_contribution(d$X, d$Y, key_vars = c("age", "sex"), target_var = "smoker",
                                     groups = list(c("age"))), "named list")
  expect_error(rapid_qi_contribution(d$X, d$Y, key_vars = c("age", "sex"), target_var = "smoker",
                                     groups = list(g = c("age", "foo"))), "not in")
  expect_error(rapid_qi_contribution(d$X, d$Y, key_vars = c("age", "sex"), target_var = "smoker",
                                     groups = list(g = c("age", "sex"))), "all quasi")
})

test_that("print, summary and plot methods work", {
  d <- .make_contrib_data()
  res <- rapid_qi_contribution(d$X, d$Y, key_vars = c("age", "sex", "region"), target_var = "smoker",
                               model_type = "cart", n_permutations = 9, seed = 1)
  expect_output(print(res), "risk contribution")
  expect_s3_class(plot(res), "ggplot")
  expect_s3_class(plot(res, measure = "rapid"), "ggplot")
  both <- plot(res, measure = c("excess", "rapid"))
  expect_length(both, 2)
  expect_s3_class(both[[2]], "ggplot")
  sm <- summary(res)
  expect_s3_class(sm, "summary.rapid_contribution")
  expect_identical(utils::capture.output(print(sm)), format(sm))
  expect_identical(as.data.frame(res), res$results)
})

test_that("QI contribution rows equal separate rapid() calls", {
  d <- .make_contrib_data()
  qi <- c("age", "sex", "region")
  res <- rapid_qi_contribution(d$X, d$Y, key_vars = qi, target_var = "smoker",
                               groups = list("age+sex" = c("age", "sex")),
                               model_type = "cart", n_permutations = 9, seed = 1)
  r <- res$results
  scenarios <- list("(none)" = qi, age = c("sex", "region"), sex = c("age", "region"),
                    region = c("age", "sex"), "age+sex" = "region")
  for (nm in names(scenarios)) {
    one <- rapid(d$X, d$Y, key_vars = scenarios[[nm]], target_var = "smoker",
                 model_type = "cart", n_perm = 9, seed = 1)
    row <- r[r$dropped == nm, ]
    expect_identical(row$rapid, one$rapid, info = nm)
    expect_identical(row$chance, one$chance, info = nm)
    expect_identical(row$excess, one$excess, info = nm)
  }
  expect_identical(res$rapid, r$rapid[1])
  expect_identical(res$rule, rapid(d$X, d$Y, key_vars = qi, target_var = "smoker",
                                   model_type = "cart", n_perm = 0)$rule)
})

test_that("rapid_qi_contribution accepts a synth_pair and the old names", {
  d <- .make_contrib_data()
  qi <- c("age", "sex", "region")
  ref <- rapid_qi_contribution(d$X, d$Y, key_vars = qi, target_var = "smoker",
                               model_type = "cart", n_permutations = 3, seed = 1)
  pair <- synth_pair(d$X, d$Y, key_vars = qi, target_var = "smoker")
  expect_identical(
    rapid_qi_contribution(pair, model_type = "cart", n_permutations = 3, seed = 1)$results,
    ref$results)
  expect_identical(
    rapid_qi_contribution(original_data = d$X, synthetic_data = d$Y, quasi_identifiers = qi,
                          sensitive_attribute = "smoker", model_type = "cart",
                          n_permutations = 3, seed = 1)$results,
    ref$results)
  expect_error(rapid_qi_contribution(d$X, d$Y, key_vars = qi, target_var = "smoker",
                                     num_taus = 0.2), class = "rapid_error_input")
})

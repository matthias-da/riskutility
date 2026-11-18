# Chance level and permutation threshold of rapid()

qi <- c("x1", "x2")

test_that("without information RAPID is close to the chance level", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(150, seed = 1)
  s <- noisy_copy(o, seed = 2, shuffle_target = TRUE)
  r <- rapid(o, s, qi, "y", model_type = "lm", n_perm = 20, seed = 1)
  expect_lt(abs(r$rapid - r$chance), 0.05)

  oc <- sim_categorical(150, seed = 1)
  sc <- noisy_copy(oc, seed = 2, shuffle_target = TRUE)
  rc <- rapid(oc, sc, qi, "y", model_type = "cart", n_perm = 20, seed = 1)
  expect_lt(abs(rc$rapid - rc$chance), 0.05)
})

test_that("with information RAPID exceeds the chance level", {
  o <- sim_numeric(150, seed = 1, signal = 0.9)
  s <- noisy_copy(o, seed = 2)
  r <- rapid(o, s, qi, "y", model_type = "lm", seed = 1)
  expect_gt(r$excess, 0.1)
  expect_equal(r$excess, r$rapid - r$chance)
})

test_that("n_perm controls the permutations", {
  o <- sim_numeric(100, seed = 1)
  s <- noisy_copy(o, seed = 2)
  r3 <- rapid(o, s, qi, "y", model_type = "lm", n_perm = 3, seed = 1)
  expect_length(r3$chance_rates, 3)
  expect_equal(r3$n_perm, 3)
  expect_equal(r3$chance, mean(r3$chance_rates))
  r0 <- rapid(o, s, qi, "y", model_type = "lm", n_perm = 0)
  expect_true(is.na(r0$chance))
  expect_true(is.na(r0$excess))
  expect_length(r0$chance_rates, 0)
  expect_error(
    rapid(o, s, qi, "y", model_type = "lm", n_perm = 0, threshold_rule = "permutation"),
    class = "rapid_error_input"
  )
  expect_error(rapid(o, s, qi, "y", model_type = "lm", n_perm = -1),
               class = "rapid_error_input")
  expect_error(rapid(o, s, qi, "y", model_type = "lm", perm_level = 1),
               class = "rapid_error_input")
})

test_that("chance level is close to Oscar's", {
  skip_on_cran()
  fx <- rapid_fixture()
  y <- fx$num$o$y
  s <- fx$meta$settings
  rule <- .rapid_rule("band", s$tau, s$metric, s$epsilon, s$delta)
  attack <- list(pred = list(point = fx$point$rf, dist = .pdist_grid(fx$D$rf)), y = y,
                 is_categorical = FALSE)
  set.seed(1)
  null <- .rapid_null(attack, rule, 200)
  expect_equal(dim(null), c(length(y), 200L))
  chance <- mean(colMeans(null > s$tau))
  expect_lt(abs(chance - fx$chance$rf), 0.05)
})

test_that("permutation threshold flags about 1 - perm_level under no information", {
  o <- sim_numeric(150, seed = 1)
  s <- noisy_copy(o, seed = 2, shuffle_target = TRUE)
  r <- rapid(o, s, qi, "y", model_type = "lm", n_perm = 20,
             threshold_rule = "permutation", perm_level = 0.95, seed = 1)
  expect_lt(r$rapid, 0.10)
  expect_equal(r$threshold_rule, "permutation")
  expect_equal(r$perm_level, 0.95)
  expect_equal(r$threshold, r$rule$threshold)
  expect_false(r$threshold == 0.3)
  # the chance level is measured at the threshold that is used
  expect_lt(abs(r$chance - 0.05), 0.03)
  r_fixed <- rapid(o, s, qi, "y", model_type = "lm", seed = 1)
  expect_equal(r_fixed$threshold_rule, "fixed")
  expect_equal(r_fixed$threshold, 0.3)
})

test_that("permutation threshold is refused for point, RCS_conditional and NCE", {
  skip_if_not_installed("rpart")
  o <- sim_numeric(100, seed = 1)
  s <- noisy_copy(o, seed = 2)
  expect_error(rapid(o, s, qi, "y", model_type = "lm", num_eval_method = "point",
                     threshold_rule = "permutation"), class = "rapid_error_input")
  oc <- sim_categorical(100, seed = 1)
  sc <- noisy_copy(oc, seed = 2)
  for (m in c("RCS_conditional", "NCE")) {
    expect_error(rapid(oc, sc, qi, "y", model_type = "cart", cat_eval_method = m,
                       threshold_rule = "permutation"), class = "rapid_error_input")
  }
  expect_no_error(rapid(oc, sc, qi, "y", model_type = "cart", n_perm = 5,
                        threshold_rule = "permutation"))
})

test_that("same seed gives the same chance level", {
  o <- sim_numeric(100, seed = 1)
  s <- noisy_copy(o, seed = 2)
  a <- rapid(o, s, qi, "y", model_type = "lm", seed = 7)
  b <- rapid(o, s, qi, "y", model_type = "lm", seed = 7)
  d <- rapid(o, s, qi, "y", model_type = "lm", seed = 8)
  expect_identical(a$chance_rates, b$chance_rates)
  expect_false(identical(a$chance_rates, d$chance_rates))
})

test_that("format shows chance and excess, and the sentence only when indistinguishable", {
  o <- sim_numeric(150, seed = 1)
  s <- noisy_copy(o, seed = 2, shuffle_target = TRUE)
  r <- rapid(o, s, qi, "y", model_type = "lm", n_perm = 10, seed = 1)
  f <- format(r)
  expect_true(any(grepl("^Chance level", f)))
  expect_true(any(grepl("^Excess", f)))
  expect_true(any(grepl("not distinguishable from the chance level", f)))
  r0 <- rapid(o, s, qi, "y", model_type = "lm", n_perm = 0)
  f0 <- format(r0)
  expect_false(any(grepl("Chance|Excess|chance", f0)))
  o2 <- sim_numeric(150, seed = 1, signal = 0.9)
  r2 <- rapid(o2, noisy_copy(o2, seed = 2), qi, "y", model_type = "lm", seed = 1)
  expect_false(any(grepl("not distinguishable", format(r2))))
  r1 <- rapid(o, s, qi, "y", model_type = "lm", n_perm = 1, seed = 1)
  expect_false(any(grepl("not distinguishable", format(r1))))
})

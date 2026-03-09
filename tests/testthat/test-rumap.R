# Tests for rumap (Multivariate Risk-Utility Map)

# Shared test data
make_test_data <- function(n = 100, seed = 42) {
  set.seed(seed)
  original <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE)),
    income = round(rnorm(n, 50000, 15000))
  )
  synthetic <- data.frame(
    age = sample(18:80, n, replace = TRUE),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE)),
    income = round(rnorm(n, 50000, 15000))
  )
  list(original = original, synthetic = synthetic)
}

test_that("rumap returns correct S3 class", {
  d <- make_test_data()
  result <- rumap(d$original, d$synthetic,
                  risk_measures = c("ims"),
                  utility_measures = c("pmse"),
                  key_vars = c("age", "gender", "region"),
                  target_var = "income",
                  seed = 42)

  expect_s3_class(result, "rumap")
  expect_true("risk" %in% names(result))
  expect_true("utility" %in% names(result))
  expect_true("normalized" %in% names(result))
  expect_true("composites" %in% names(result))
  expect_true("pareto" %in% names(result))
  expect_true("n_sdgs" %in% names(result))
  expect_true("sdg_names" %in% names(result))
  expect_true("risk_measures" %in% names(result))
  expect_true("utility_measures" %in% names(result))
  expect_true("metadata" %in% names(result))
})

test_that("rumap works with single SDG (data.frame input)", {
  d <- make_test_data()
  result <- rumap(d$original, d$synthetic,
                  risk_measures = "ims",
                  utility_measures = "pmse",
                  seed = 42)

  expect_equal(result$n_sdgs, 1)
  expect_equal(nrow(result$risk), 1)
  expect_equal(nrow(result$utility), 1)
})

test_that("rumap works with multiple SDGs", {
  d <- make_test_data()
  d2 <- make_test_data(seed = 99)
  result <- rumap(d$original,
                  list(sdg1 = d$synthetic, sdg2 = d2$synthetic),
                  risk_measures = "ims",
                  utility_measures = "pmse",
                  seed = 42)

  expect_equal(result$n_sdgs, 2)
  expect_equal(nrow(result$risk), 2)
  expect_equal(nrow(result$utility), 2)
  expect_equal(result$sdg_names, c("sdg1", "sdg2"))
})

test_that("rumap validates risk measures", {
  d <- make_test_data()
  expect_warning(
    rumap(d$original, d$synthetic,
          risk_measures = c("ims", "nonexistent"),
          utility_measures = "pmse",
          seed = 42),
    "Unknown risk measures"
  )
})

test_that("rumap validates utility measures", {
  d <- make_test_data()
  expect_warning(
    rumap(d$original, d$synthetic,
          risk_measures = "ims",
          utility_measures = c("pmse", "bogus"),
          seed = 42),
    "Unknown utility measures"
  )
})

test_that("rumap warns when key_vars missing for attribution", {
  d <- make_test_data()
  expect_warning(
    rumap(d$original, d$synthetic,
          risk_measures = c("dcap", "ims"),
          utility_measures = "pmse",
          key_vars = NULL,
          seed = 42),
    "key_vars"
  )
})

test_that("rumap normalization works", {
  d <- make_test_data()
  d2 <- make_test_data(seed = 99)

  result <- rumap(d$original,
                  list(s1 = d$synthetic, s2 = d2$synthetic),
                  risk_measures = "ims",
                  utility_measures = "pmse",
                  normalize = TRUE,
                  seed = 42)

  expect_false(is.null(result$normalized))
  # Normalized values should be in [0, 1] or 0.5 (if all same)
  norm_vals <- result$normalized[, result$risk_measures, drop = FALSE]
  for (col in names(norm_vals)) {
    vals <- norm_vals[[col]]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {
      expect_true(all(vals >= 0 & vals <= 1))
    }
  }
})

test_that("rumap without normalization returns NULL normalized", {
  d <- make_test_data()
  result <- rumap(d$original, d$synthetic,
                  risk_measures = "ims",
                  utility_measures = "pmse",
                  normalize = FALSE,
                  seed = 42)

  expect_null(result$normalized)
})

test_that("rumap computes Pareto frontier", {
  d <- make_test_data()
  result <- rumap(d$original, d$synthetic,
                  risk_measures = "ims",
                  utility_measures = "pmse",
                  seed = 42)

  expect_true(is.logical(result$pareto))
  expect_equal(length(result$pareto), result$n_sdgs)
  # Single SDG should be Pareto-optimal
  expect_true(result$pareto[1])
})

test_that("rumap composites contain risk and utility means", {
  d <- make_test_data()
  result <- rumap(d$original, d$synthetic,
                  risk_measures = "ims",
                  utility_measures = "pmse",
                  seed = 42)

  expect_true("risk_mean" %in% names(result$composites))
  expect_true("utility_mean" %in% names(result$composites))
})

test_that("rumap includes rapid as valid risk measure", {
  d <- make_test_data()
  result <- rumap(d$original, d$synthetic,
                  risk_measures = c("rapid", "ims"),
                  utility_measures = "pmse",
                  key_vars = c("age", "gender", "region"),
                  target_var = "income",
                  seed = 42)

  expect_true("rapid" %in% result$risk_measures)
  expect_false(is.na(result$risk[1, "rapid"]))
})

test_that("rumap print method works", {
  d <- make_test_data()
  result <- rumap(d$original, d$synthetic,
                  risk_measures = "ims",
                  utility_measures = "pmse",
                  seed = 42)

  expect_output(print(result), "Risk-Utility Map")
})

test_that("rumap summary method works", {
  d <- make_test_data()
  result <- rumap(d$original, d$synthetic,
                  risk_measures = "ims",
                  utility_measures = "pmse",
                  seed = 42)

  s <- summary(result)
  expect_s3_class(s, "summary.rumap")
})

test_that("rumap plot method runs without error", {
  d <- make_test_data()
  result <- rumap(d$original,
                  list(s1 = d$synthetic, s2 = make_test_data(seed = 99)$synthetic),
                  risk_measures = "ims",
                  utility_measures = "pmse",
                  seed = 42)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(result, which = 1))
})

test_that("rumap synth_pair dispatch works", {
  d <- make_test_data()
  pair <- synth_pair(d$original, d$synthetic,
                     key_vars = c("age", "gender", "region"),
                     target_var = "income")

  result <- rumap(pair,
                  risk_measures = "ims",
                  utility_measures = "pmse")

  expect_s3_class(result, "rumap")
})

test_that("rumap errors on invalid input", {
  expect_error(rumap("not_a_df", data.frame(a = 1)), "data.frame")
  d <- make_test_data()
  expect_error(rumap(d$original, list()), "non-empty")
})

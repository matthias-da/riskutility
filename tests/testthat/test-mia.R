# Tests for mia_classifier (Membership Inference Attack)

library(testthat)

skip_if_not_installed("ranger")

# --- Setup: shared test data ---

make_mia_data <- function(n = 100, seed = 42) {
  set.seed(seed)
  list(
    real = data.frame(
      x1 = rnorm(n),
      x2 = factor(sample(c("A", "B"), n, replace = TRUE))
    ),
    synt = data.frame(
      x1 = rnorm(n),
      x2 = factor(sample(c("A", "B"), n, replace = TRUE))
    ),
    hout = data.frame(
      x1 = rnorm(n),
      x2 = factor(sample(c("A", "B"), n, replace = TRUE))
    )
  )
}

# --- Class structure ---

test_that("mia_classifier returns correct S3 class", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  expect_s3_class(result, "mia")
})

test_that("mia_classifier result has all expected fields", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)

  expect_true("precision" %in% names(result))
  expect_true("precision_se" %in% names(result))
  expect_true("recall" %in% names(result))
  expect_true("recall_se" %in% names(result))
  expect_true("macro_f1" %in% names(result))
  expect_true("macro_f1_se" %in% names(result))
  expect_true("privacy_pass" %in% names(result))
  expect_true("num_eval_iter" %in% names(result))
  expect_true("method" %in% names(result))
})

test_that("privacy_pass is logical", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  expect_type(result$privacy_pass, "logical")
})

test_that("method is stored correctly", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           method = "rf", num_eval_iter = 2, seed = 123)
  expect_equal(result$method, "rf")
})

test_that("num_eval_iter is stored correctly", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 3, seed = 123)
  expect_equal(result$num_eval_iter, 3)
})

# --- Metric value ranges ---

test_that("precision is between 0 and 1", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  expect_true(result$precision >= 0 && result$precision <= 1)
})

test_that("recall is between 0 and 1", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  expect_true(result$recall >= 0 && result$recall <= 1)
})

test_that("macro_f1 is between 0 and 1", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  expect_true(result$macro_f1 >= 0 && result$macro_f1 <= 1)
})

test_that("standard errors are non-negative", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 3, seed = 123)
  expect_true(result$precision_se >= 0)
  expect_true(result$recall_se >= 0)
  expect_true(result$macro_f1_se >= 0)
})

# --- Privacy pass logic ---

test_that("privacy_pass is TRUE when recall <= 0.55", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  if (result$recall <= 0.55) {
    expect_true(result$privacy_pass)
  } else {
    expect_false(result$privacy_pass)
  }
})

# --- Independent data should have recall near 0.5 ---

test_that("recall is near 0.5 for independent random data", {
  set.seed(99)
  n <- 200
  real <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  synt <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  hout <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))

  result <- mia_classifier(real, synt_data = synt, hout_data = hout,
                           num_eval_iter = 3, seed = 42)
  # For truly independent data, recall should be roughly 0.5
  expect_true(result$recall > 0.2 && result$recall < 0.8)
})

# --- Reproducibility via seed ---

test_that("mia_classifier produces reproducible results with same seed", {
  d <- make_mia_data()
  r1 <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                        num_eval_iter = 2, seed = 42)
  r2 <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                        num_eval_iter = 2, seed = 42)
  expect_equal(r1$precision, r2$precision)
  expect_equal(r1$recall, r2$recall)
  expect_equal(r1$macro_f1, r2$macro_f1)
})

# --- synth_pair dispatch ---

test_that("mia_classifier works via synth_pair", {
  d <- make_mia_data()
  pair <- synth_pair(original = d$real, synthetic = d$synt, holdout = d$hout)
  result <- mia_classifier(pair, num_eval_iter = 2, seed = 123)
  expect_s3_class(result, "mia")
  expect_true(!is.na(result$recall))
})

test_that("mia_classifier.synth_pair errors on sdcMicro source", {
  d <- make_mia_data()
  pair <- synth_pair(original = d$real, synthetic = d$synt, holdout = d$hout,
                     source = "sdcMicro")
  expect_error(mia_classifier(pair),
               "not applicable to.*traditionally anonymized")
})

test_that("mia_classifier.synth_pair errors without holdout", {
  d <- make_mia_data()
  pair <- synth_pair(original = d$real, synthetic = d$synt)
  expect_error(mia_classifier(pair),
               "mia_classifier requires holdout data")
})

# --- Input validation ---

test_that("mia_classifier errors on non-data.frame inputs", {
  d <- make_mia_data()
  expect_error(mia_classifier(1:10, synt_data = d$synt, hout_data = d$hout),
               "data.frames")
  expect_error(mia_classifier(d$real, synt_data = 1:10, hout_data = d$hout),
               "data.frames")
  expect_error(mia_classifier(d$real, synt_data = d$synt, hout_data = "bad"),
               "data.frames")
})

test_that("mia_classifier errors on unknown method", {
  d <- make_mia_data()
  expect_error(mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                               method = "unknown"),
               "Unknown method")
})

test_that("mia_classifier errors on invalid num_eval_iter", {
  d <- make_mia_data()
  expect_error(mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                               num_eval_iter = 0),
               "num_eval_iter must be >= 1")
})

test_that("mia_classifier errors on no common columns", {
  real <- data.frame(a = 1:10)
  synt <- data.frame(b = 1:10)
  hout <- data.frame(c = 1:10)
  expect_error(mia_classifier(real, synt_data = synt, hout_data = hout),
               "No common columns")
})

test_that("mia_classifier errors on holdout too small", {
  d <- make_mia_data()
  tiny_hout <- data.frame(x1 = 1, x2 = factor("A"))
  expect_error(mia_classifier(d$real, synt_data = d$synt, hout_data = tiny_hout),
               "Holdout data too small")
})

test_that("mia_classifier errors on invalid cols argument", {
  d <- make_mia_data()
  expect_error(mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                               cols = c("x1", "nonexistent")),
               "cols not found")
})

# --- Column selection ---

test_that("mia_classifier works with explicit cols", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           cols = "x1", num_eval_iter = 2, seed = 123)
  expect_s3_class(result, "mia")
})

test_that("mia_classifier uses only common columns across all three datasets", {
  set.seed(42)
  real <- data.frame(x1 = rnorm(50), x2 = rnorm(50), extra = rnorm(50))
  synt <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  hout <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  # extra is only in real, so it should be dropped automatically
  result <- mia_classifier(real, synt_data = synt, hout_data = hout,
                           num_eval_iter = 2, seed = 123)
  expect_s3_class(result, "mia")
})

# --- S3 methods: print, summary, plot ---

test_that("print.mia runs without error", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  expect_output(print(result), "Membership Inference Attack")
  expect_output(print(result), "Precision")
  expect_output(print(result), "Recall")
  expect_output(print(result), "Macro F1")
  expect_output(print(result), "Privacy pass")
})

test_that("print.mia returns object invisibly", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  out <- capture.output(ret <- print(result))
  expect_identical(ret, result)
})

test_that("summary.mia returns correct class", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  s <- summary(result)
  expect_s3_class(s, "summary.mia")
  expect_true(!is.null(s$precision))
  expect_true(!is.null(s$recall))
  expect_true(!is.null(s$macro_f1))
  expect_true(!is.null(s$recall_excess))
  expect_true(!is.null(s$num_eval_iter))
  expect_true(!is.null(s$method))
})

test_that("summary.mia recall_excess is recall - 0.5", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  s <- summary(result)
  expect_equal(s$recall_excess, result$recall - 0.5)
})

test_that("print.summary.mia runs without error", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  s <- summary(result)
  expect_output(print(s), "Summary: Membership Inference Attack")
  expect_output(print(s), "Metrics")
  expect_output(print(s), "Recall excess")
})

test_that("print.summary.mia returns object invisibly", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  s <- summary(result)
  out <- capture.output(ret <- print(s))
  expect_identical(ret, s)
})

test_that("plot.mia works for which = 1", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(result, which = 1))
  dev.off()
})

test_that("plot.mia returns object invisibly", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 2, seed = 123)
  pdf(tempfile(fileext = ".pdf"))
  ret <- plot(result, which = 1)
  dev.off()
  expect_identical(ret, result)
})

# --- NA handling ---

test_that("mia_classifier warns on NAs in data", {
  set.seed(42)
  n <- 50
  real <- data.frame(x1 = c(rnorm(n - 2), NA, NA), x2 = rnorm(n))
  synt <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  hout <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  expect_warning(
    mia_classifier(real, synt_data = synt, hout_data = hout,
                   num_eval_iter = 2, seed = 123),
    "NAs detected"
  )
})

# --- Numeric-only data ---

test_that("mia_classifier works with numeric-only data", {
  set.seed(42)
  n <- 80
  real <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  synt <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  hout <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  result <- mia_classifier(real, synt_data = synt, hout_data = hout,
                           num_eval_iter = 2, seed = 123)
  expect_s3_class(result, "mia")
})

# --- Categorical-only data ---

test_that("mia_classifier works with categorical-only data", {
  set.seed(42)
  n <- 80
  real <- data.frame(
    x1 = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    x2 = factor(sample(c("X", "Y"), n, replace = TRUE))
  )
  synt <- data.frame(
    x1 = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    x2 = factor(sample(c("X", "Y"), n, replace = TRUE))
  )
  hout <- data.frame(
    x1 = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    x2 = factor(sample(c("X", "Y"), n, replace = TRUE))
  )
  result <- mia_classifier(real, synt_data = synt, hout_data = hout,
                           num_eval_iter = 2, seed = 123)
  expect_s3_class(result, "mia")
})

# --- Single evaluation iteration ---

test_that("mia_classifier works with num_eval_iter = 1 and SE is NA", {
  d <- make_mia_data()
  result <- mia_classifier(d$real, synt_data = d$synt, hout_data = d$hout,
                           num_eval_iter = 1, seed = 123)
  expect_s3_class(result, "mia")
  # With only 1 iteration, SE should be NA
  expect_true(is.na(result$precision_se))
  expect_true(is.na(result$recall_se))
  expect_true(is.na(result$macro_f1_se))
})

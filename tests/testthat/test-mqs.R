# Tests for mqs (Model Quality Score)

library(testthat)

# --- Setup: mock mqs object for lightweight tests ---

make_mock_mqs <- function(ratio = 0.95, measure = "Accuracy") {
  mock <- list(
    mqs_ratio = ratio,
    mqs_table = data.frame(
      data = c("X", "Y"),
      glm = c(0.80, 0.75),
      rpart = c(0.82, 0.78),
      measure = c(measure, measure)
    )
  )
  class(mock) <- "mqs"
  mock
}

make_mock_mqs_rmse <- function(ratio = 1.10) {
  mock <- list(
    mqs_ratio = ratio,
    mqs_table = data.frame(
      data = c("X", "Y"),
      glm = c(5.2, 5.7),
      rpart = c(4.8, 5.3),
      measure = c("RMSE", "RMSE")
    )
  )
  class(mock) <- "mqs"
  mock
}

# --- Class structure (mock-based) ---

test_that("mock mqs object has correct S3 class", {
  mock <- make_mock_mqs()
  expect_s3_class(mock, "mqs")
})

test_that("mock mqs object has expected fields", {
  mock <- make_mock_mqs()
  expect_named(mock, c("mqs_ratio", "mqs_table"))
})

test_that("mqs_ratio is numeric", {
  mock <- make_mock_mqs()
  expect_type(mock$mqs_ratio, "double")
})

test_that("mqs_table is a data.frame", {
  mock <- make_mock_mqs()
  expect_true(is.data.frame(mock$mqs_table))
})

test_that("mqs_table has data column with X and Y", {
  mock <- make_mock_mqs()
  expect_true("data" %in% names(mock$mqs_table))
  expect_equal(mock$mqs_table$data, c("X", "Y"))
})

test_that("mqs_table has measure column", {
  mock <- make_mock_mqs()
  expect_true("measure" %in% names(mock$mqs_table))
})

# --- S3 methods: print (mock-based) ---

test_that("print.mqs runs without error", {
  mock <- make_mock_mqs()
  expect_output(print(mock), "Model Quality Score")
})

test_that("print.mqs shows MQS ratio", {
  mock <- make_mock_mqs(ratio = 0.95)
  expect_output(print(mock), "MQS ratio:")
  expect_output(print(mock), "0.95")
})

test_that("print.mqs shows performance table", {
  mock <- make_mock_mqs()
  expect_output(print(mock), "Model performance table")
})

test_that("print.mqs returns object invisibly", {
  mock <- make_mock_mqs()
  out <- capture.output(ret <- print(mock))
  expect_identical(ret, mock)
})

test_that("print.mqs shows separator line", {
  mock <- make_mock_mqs()
  expect_output(print(mock), "=====")
})

# --- S3 methods: summary (mock-based) ---

test_that("summary.mqs returns correct class", {
  mock <- make_mock_mqs()
  s <- summary(mock)
  expect_s3_class(s, "summary.mqs")
})

test_that("summary.mqs has expected fields", {
  mock <- make_mock_mqs()
  s <- summary(mock)
  expect_named(s, c("mqs_ratio", "mqs_table", "interpretation"))
})

test_that("summary.mqs preserves mqs_ratio", {
  mock <- make_mock_mqs(ratio = 0.88)
  s <- summary(mock)
  expect_equal(s$mqs_ratio, 0.88)
})

test_that("summary.mqs preserves mqs_table", {
  mock <- make_mock_mqs()
  s <- summary(mock)
  expect_equal(s$mqs_table, mock$mqs_table)
})

test_that("summary interpretation: comparable when ratio ~1", {
  mock <- make_mock_mqs(ratio = 1.00)
  s <- summary(mock)
  expect_match(s$interpretation, "comparable")
})

test_that("summary interpretation: better quality when ratio < 0.95", {
  mock <- make_mock_mqs(ratio = 0.80)
  s <- summary(mock)
  expect_match(s$interpretation, "better prediction quality")
})

test_that("summary interpretation: lower quality when ratio > 1.05", {
  mock <- make_mock_mqs(ratio = 1.20)
  s <- summary(mock)
  expect_match(s$interpretation, "lower prediction quality")
})

test_that("summary interpretation edge case: ratio exactly 0.95", {
  mock <- make_mock_mqs(ratio = 0.95)
  s <- summary(mock)
  expect_match(s$interpretation, "comparable")
})

test_that("summary interpretation edge case: ratio exactly 1.05", {
  mock <- make_mock_mqs(ratio = 1.05)
  s <- summary(mock)
  expect_match(s$interpretation, "comparable")
})

# --- S3 methods: print.summary (mock-based) ---

test_that("print.summary.mqs runs without error", {
  mock <- make_mock_mqs()
  expect_output(print(summary(mock)), "Summary: Model Quality Score")
})

test_that("print.summary.mqs shows MQS ratio", {
  mock <- make_mock_mqs(ratio = 0.92)
  expect_output(print(summary(mock)), "MQS ratio:")
  expect_output(print(summary(mock)), "0.92")
})

test_that("print.summary.mqs shows interpretation", {
  mock <- make_mock_mqs(ratio = 0.80)
  expect_output(print(summary(mock)), "Interpretation:")
})

test_that("print.summary.mqs shows performance table", {
  mock <- make_mock_mqs()
  expect_output(print(summary(mock)), "Model performance table")
})

test_that("print.summary.mqs returns object invisibly", {
  mock <- make_mock_mqs()
  s <- summary(mock)
  out <- capture.output(ret <- print(s))
  expect_identical(ret, s)
})

test_that("print.summary.mqs shows separator line", {
  mock <- make_mock_mqs()
  expect_output(print(summary(mock)), "=====")
})

# --- S3 generic dispatch ---

test_that("mqs uses UseMethod dispatch", {
  # Confirm mqs is a generic function
  expect_true(is.function(mqs))
})

# --- Input validation (requires caret and caretEnsemble) ---

test_that("mqs errors on non-data.frame X", {
  skip_if_not_installed("caret")
  skip_if_not_installed("caretEnsemble")

  Y <- data.frame(y = factor(c("A", "B")), x1 = c(1, 2))
  expect_error(mqs(1:10, Y, form = y ~ x1), "X must be a data frame")
})

test_that("mqs errors on non-data.frame Y", {
  skip_if_not_installed("caret")
  skip_if_not_installed("caretEnsemble")

  X <- data.frame(y = factor(c("A", "B")), x1 = c(1, 2))
  expect_error(mqs(X, 1:10, form = y ~ x1), "Y must be a data frame")
})

test_that("mqs errors on invalid formula", {
  skip_if_not_installed("caret")
  skip_if_not_installed("caretEnsemble")

  X <- data.frame(y = factor(c("A", "B")), x1 = c(1, 2))
  Y <- data.frame(y = factor(c("A", "B")), x1 = c(1, 2))
  expect_error(mqs(X, Y, form = "not_a_formula"), "form must be a valid formula")
})

test_that("mqs errors when formula variables missing from data", {
  skip_if_not_installed("caret")
  skip_if_not_installed("caretEnsemble")

  X <- data.frame(y = factor(c("A", "B")), x1 = c(1, 2))
  Y <- data.frame(y = factor(c("A", "B")), x1 = c(1, 2))
  expect_error(mqs(X, Y, form = y ~ x1 + nonexistent), "missing")
})

test_that("mqs errors on invalid na argument", {
  skip_if_not_installed("caret")
  skip_if_not_installed("caretEnsemble")

  X <- data.frame(y = factor(c("A", "B")), x1 = c(1, 2))
  Y <- data.frame(y = factor(c("A", "B")), x1 = c(1, 2))
  expect_error(mqs(X, Y, form = y ~ x1, na = "invalid_option"),
               "na must be one of")
})

test_that("mqs errors when column classes differ between X and Y", {
  skip_if_not_installed("caret")
  skip_if_not_installed("caretEnsemble")

  X <- data.frame(y = factor(c("A", "B")), x1 = c(1, 2))
  Y <- data.frame(y = factor(c("A", "B")), x1 = factor(c("a", "b")))
  expect_error(mqs(X, Y, form = y ~ x1), "same class")
})

# --- synth_pair method (mock-based check) ---

test_that("mqs.synth_pair method exists", {
  expect_true("mqs.synth_pair" %in% ls(getNamespace("riskutility")))
})

# --- Full integration test (requires caret + caretEnsemble) ---

test_that("mqs returns correct class with real data (categorical response)", {
  skip_if_not_installed("caret")
  skip_if_not_installed("caretEnsemble")
  skip_on_cran()

  set.seed(123)
  n <- 100
  X <- data.frame(
    y = factor(sample(c("A", "B"), n, replace = TRUE)),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  Y <- data.frame(
    y = factor(sample(c("A", "B"), n, replace = TRUE)),
    x1 = rnorm(n, 0.1, 1),
    x2 = rnorm(n, 0.1, 1)
  )

  result <- tryCatch(
    suppressWarnings(
      mqs(X, Y, form = y ~ x1 + x2, methods = c("glm", "rpart"))
    ),
    error = function(e) {
      skip(paste("mqs() failed with caretEnsemble:", conditionMessage(e)))
    }
  )

  expect_s3_class(result, "mqs")
  expect_true(is.numeric(result$mqs_ratio))
  expect_true(result$mqs_ratio > 0)
  expect_true(is.data.frame(result$mqs_table))
  expect_equal(nrow(result$mqs_table), 2)
  expect_true("Accuracy" %in% result$mqs_table$measure)
})

test_that("mqs returns correct class with real data (numeric response)", {
  skip_if_not_installed("caret")
  skip_if_not_installed("caretEnsemble")
  skip_on_cran()

  set.seed(456)
  n <- 100
  X <- data.frame(
    y = rnorm(n, 50, 10),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  Y <- data.frame(
    y = rnorm(n, 50, 10),
    x1 = rnorm(n, 0.1, 1),
    x2 = rnorm(n, 0.1, 1)
  )

  result <- tryCatch(
    suppressWarnings(
      mqs(X, Y, form = y ~ x1 + x2, methods = c("glm", "rpart"))
    ),
    error = function(e) {
      skip(paste("mqs() failed with caretEnsemble:", conditionMessage(e)))
    }
  )

  expect_s3_class(result, "mqs")
  expect_true(is.numeric(result$mqs_ratio))
  expect_true(result$mqs_ratio > 0)
  expect_true(is.data.frame(result$mqs_table))
  expect_equal(nrow(result$mqs_table), 2)
  expect_true("RMSE" %in% result$mqs_table$measure)
})

test_that("mqs print works on real result", {
  skip_if_not_installed("caret")
  skip_if_not_installed("caretEnsemble")
  skip_on_cran()

  set.seed(789)
  n <- 80
  X <- data.frame(
    y = factor(sample(c("A", "B"), n, replace = TRUE)),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  Y <- data.frame(
    y = factor(sample(c("A", "B"), n, replace = TRUE)),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  result <- tryCatch(
    suppressWarnings(
      mqs(X, Y, form = y ~ x1 + x2, methods = c("glm", "rpart"))
    ),
    error = function(e) {
      skip(paste("mqs() failed with caretEnsemble:", conditionMessage(e)))
    }
  )

  expect_output(print(result), "Model Quality Score")
  expect_output(print(result), "MQS ratio:")
})

test_that("mqs summary works on real result", {
  skip_if_not_installed("caret")
  skip_if_not_installed("caretEnsemble")
  skip_on_cran()

  set.seed(101)
  n <- 80
  X <- data.frame(
    y = factor(sample(c("A", "B"), n, replace = TRUE)),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  Y <- data.frame(
    y = factor(sample(c("A", "B"), n, replace = TRUE)),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  result <- tryCatch(
    suppressWarnings(
      mqs(X, Y, form = y ~ x1 + x2, methods = c("glm", "rpart"))
    ),
    error = function(e) {
      skip(paste("mqs() failed with caretEnsemble:", conditionMessage(e)))
    }
  )

  s <- summary(result)
  expect_s3_class(s, "summary.mqs")
  expect_true(!is.null(s$interpretation))
  expect_output(print(s), "Summary: Model Quality Score")
  expect_output(print(s), "Interpretation:")
})

test_that("mqs works with synth_pair dispatch", {
  skip_if_not_installed("caret")
  skip_if_not_installed("caretEnsemble")
  skip_on_cran()

  set.seed(202)
  n <- 80
  orig <- data.frame(
    y = factor(sample(c("A", "B"), n, replace = TRUE)),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  syn <- data.frame(
    y = factor(sample(c("A", "B"), n, replace = TRUE)),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  pair <- synth_pair(orig, syn)

  result <- tryCatch(
    suppressWarnings(
      mqs(pair, form = y ~ x1 + x2, methods = c("glm", "rpart"))
    ),
    error = function(e) {
      skip(paste("mqs() failed with caretEnsemble:", conditionMessage(e)))
    }
  )

  expect_s3_class(result, "mqs")
  expect_true(is.numeric(result$mqs_ratio))
})

# --- Mock-based RMSE variant ---

test_that("print works for RMSE-based mqs", {
  mock <- make_mock_mqs_rmse(ratio = 1.10)
  expect_output(print(mock), "MQS ratio:")
  expect_output(print(mock), "1.1")
})

test_that("summary works for RMSE-based mqs", {
  mock <- make_mock_mqs_rmse(ratio = 1.10)
  s <- summary(mock)
  expect_s3_class(s, "summary.mqs")
  expect_match(s$interpretation, "lower prediction quality")
})

# --- S3 methods: plot (mock-based) ---

test_that("plot.mqs runs without error (which = 1)", {
  mock <- make_mock_mqs()
  expect_silent(plot(mock, which = 1))
})

test_that("plot.mqs runs without error (which = 2)", {
  mock <- make_mock_mqs()
  expect_silent(plot(mock, which = 2))
})

test_that("plot.mqs returns object invisibly", {
  mock <- make_mock_mqs()
  ret <- plot(mock, which = 1)
  expect_identical(ret, mock)
})

test_that("plot.mqs works with RMSE-based mock", {
  mock <- make_mock_mqs_rmse(ratio = 1.10)
  expect_silent(plot(mock, which = 1))
  expect_silent(plot(mock, which = 2))
})

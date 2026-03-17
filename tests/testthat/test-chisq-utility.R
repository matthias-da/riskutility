# Tests for chisq_utility (Chi-Square Utility Measures)

library(testthat)

# --- Setup: shared test data ---

make_categorical_data <- function(n = 500, seed = 123) {
  set.seed(seed)
  data.frame(
    age = sample(c("young", "middle", "old"), n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
  )
}

# --- Class structure ---

test_that("chisq_utility returns correct S3 class", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  expect_s3_class(result, "chisq_utility")
})

test_that("chisq_utility result has expected fields", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))

  expect_true("chi2" %in% names(result))
  expect_true("df" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("VW" %in% names(result))
  expect_true("FT" %in% names(result))
  expect_true("G" %in% names(result))
  expect_true("JSD" %in% names(result))
  expect_true("n_cells" %in% names(result))
  expect_true("n_empty_orig" %in% names(result))
  expect_true("n_empty_synth" %in% names(result))
  expect_true("pct_utility" %in% names(result))
  expect_true("table_orig" %in% names(result))
  expect_true("table_synth" %in% names(result))
  expect_true("n_orig" %in% names(result))
  expect_true("n_synth" %in% names(result))
  expect_true("vars" %in% names(result))
})

# --- Basic computation ---

test_that("chisq_utility computes chi2 as a positive number", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  expect_true(result$chi2 >= 0)
})

test_that("JSD is bounded between 0 and 1", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  expect_true(result$JSD >= 0)
  expect_true(result$JSD <= 1)
})

test_that("p_value is between 0 and 1", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  expect_true(result$p_value >= 0)
  expect_true(result$p_value <= 1)
})

test_that("pct_utility is between 0 and 100", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  expect_true(result$pct_utility >= 0)
  expect_true(result$pct_utility <= 100)
})

test_that("FT and G statistics are non-negative", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  expect_true(result$FT >= 0)
  expect_true(result$G >= 0)
})

# --- Identical data yields near-zero statistics ---

test_that("identical data yields low VW and JSD", {
  set.seed(42)
  orig <- make_categorical_data(n = 300)
  result <- chisq_utility(orig, orig, vars = c("age", "gender"))
  expect_equal(result$VW, (result$chi2 - result$df) / 300)
  # Chi2 should be 0 when data is identical
  expect_equal(result$chi2, 0)
  expect_true(result$JSD < 0.001)
})

# --- Number of cells ---

test_that("n_cells matches product of factor levels", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  # age has 3 levels, gender has 2 levels => 6 cells
  expect_equal(result$n_cells, 6)
})

test_that("n_cells is correct for three variables", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender", "region"))
  # 3 * 2 * 4 = 24
  expect_equal(result$n_cells, 24)
})

# --- Different variable sets ---

test_that("chisq_utility works with single variable", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = "age")
  expect_s3_class(result, "chisq_utility")
  expect_equal(result$n_cells, 3)
  expect_equal(result$vars, "age")
})

test_that("chisq_utility works with character columns", {
  set.seed(123)
  orig <- data.frame(x = sample(letters[1:4], 100, replace = TRUE),
                     stringsAsFactors = FALSE)
  syn  <- data.frame(x = sample(letters[1:4], 100, replace = TRUE),
                     stringsAsFactors = FALSE)
  result <- chisq_utility(orig, syn, vars = "x")
  expect_s3_class(result, "chisq_utility")
  expect_equal(result$n_cells, 4)
})

test_that("chisq_utility handles different factor levels in orig and syn", {
  orig <- data.frame(x = factor(c("A", "A", "B", "B", "C")))
  syn  <- data.frame(x = factor(c("A", "B", "B", "D", "D")))
  result <- chisq_utility(orig, syn, vars = "x")
  # Union of levels: A, B, C, D => 4 cells
  expect_equal(result$n_cells, 4)
})

# --- Weights ---

test_that("chisq_utility works with numeric weight vectors", {
  orig <- make_categorical_data(n = 100)
  syn  <- make_categorical_data(n = 100, seed = 456)
  wt_X <- rep(1, 100)
  wt_Y <- rep(1, 100)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"),
                          weight_X = wt_X, weight_Y = wt_Y)
  expect_s3_class(result, "chisq_utility")
})

# --- synth_pair dispatch ---

test_that("chisq_utility works via synth_pair", {
  orig <- make_categorical_data(n = 200)
  syn  <- make_categorical_data(n = 200, seed = 456)
  pair <- synth_pair(orig, syn)
  result <- chisq_utility(pair, vars = c("age", "gender"))
  expect_s3_class(result, "chisq_utility")
  expect_equal(result$vars, c("age", "gender"))
})

test_that("chisq_utility.synth_pair auto-selects categorical vars when vars is NULL", {
  orig <- make_categorical_data(n = 200)
  syn  <- make_categorical_data(n = 200, seed = 456)
  pair <- synth_pair(orig, syn)
  result <- chisq_utility(pair)
  expect_s3_class(result, "chisq_utility")
  expect_true(length(result$vars) > 0)
})

# --- Error handling ---

test_that("chisq_utility errors on non-data.frame X", {
  syn <- make_categorical_data()
  expect_error(chisq_utility(1:10, syn, vars = "age"),
               "X must be a data frame")
})

test_that("chisq_utility errors on non-data.frame Y", {
  orig <- make_categorical_data()
  expect_error(chisq_utility(orig, 1:10, vars = "age"),
               "Y must be a data frame")
})

test_that("chisq_utility errors on missing variables in X", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  expect_error(chisq_utility(orig, syn, vars = c("age", "nonexistent")),
               "Variables missing in X")
})

test_that("chisq_utility errors on missing variables in Y", {
  orig <- data.frame(age = c("young", "old"), extra = c("a", "b"))
  syn  <- data.frame(age = c("young", "old"))
  expect_error(chisq_utility(orig, syn, vars = c("age", "extra")),
               "Variables missing in Y")
})

test_that("chisq_utility warns when max_cells exceeded", {
  set.seed(42)
  # Many-level factors will produce large cross-tabulation
  n <- 50
  orig <- data.frame(
    v1 = factor(1:n),
    v2 = factor(1:n)
  )
  syn <- data.frame(
    v1 = factor(1:n),
    v2 = factor(1:n)
  )
  expect_warning(
    chisq_utility(orig, syn, vars = c("v1", "v2"), max_cells = 100),
    "Cross-tabulation would have"
  )
})

# --- S3 methods: print, summary, plot ---

test_that("print.chisq_utility runs without error", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  expect_output(print(result), "Chi-Square Utility Assessment")
  expect_output(print(result), "Voas-Williamson")
  expect_output(print(result), "Freeman-Tukey")
  expect_output(print(result), "Jensen-Shannon Divergence")
})

test_that("print.chisq_utility returns object invisibly", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  out <- capture.output(ret <- print(result))
  expect_identical(ret, result)
})

test_that("summary.chisq_utility returns correct class", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  s <- summary(result)
  expect_s3_class(s, "summary.chisq_utility")
  expect_true(!is.null(s$VW))
  expect_true(!is.null(s$JSD))
  expect_true(!is.null(s$diff_summary))
  expect_true(!is.null(s$pct_utility))
})

test_that("print.summary.chisq_utility runs without error", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  s <- summary(result)
  expect_output(print(s), "Summary: Chi-Square Utility Assessment")
  expect_output(print(s), "Cell-Level Differences")
  expect_output(print(s), "Estimated Utility")
})

test_that("print.summary.chisq_utility returns object invisibly", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  s <- summary(result)
  out <- capture.output(ret <- print(s))
  expect_identical(ret, s)
})

test_that("plot.chisq_utility works for which = 1", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(result, which = 1))
  dev.off()
})

test_that("plot.chisq_utility works for which = 2", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(result, which = 2))
  dev.off()
})

test_that("plot.chisq_utility works for which = 1:2", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(result, which = 1:2))
  dev.off()
})

# --- Consistency checks ---

test_that("table_orig sums to n_orig", {
  orig <- make_categorical_data(n = 200)
  syn  <- make_categorical_data(n = 300, seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  expect_equal(sum(result$table_orig), 200)
  expect_equal(sum(result$table_synth), 300)
})

test_that("prop_orig and prop_synth sum to 1", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  expect_equal(sum(result$prop_orig), 1)
  expect_equal(sum(result$prop_synth), 1)
})

test_that("df equals number of non-zero original cells minus 1", {
  orig <- make_categorical_data()
  syn  <- make_categorical_data(seed = 456)
  result <- chisq_utility(orig, syn, vars = c("age", "gender"))
  expected_df <- sum(result$prop_orig > 0) - 1
  expect_equal(result$df, expected_df)
})

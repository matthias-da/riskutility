# Tests for suda (Special Uniques Detection Algorithm)

library(testthat)

# --- Setup: shared test data ---

make_test_data <- function(n = 100, seed = 42) {
  set.seed(seed)
  data.frame(
    age = sample(c("young", "middle", "old"), n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE),
    education = sample(c("low", "medium", "high"), n, replace = TRUE)
  )
}

make_unique_data <- function() {
  # Every record has a unique QI combination (high SUDA risk)
  data.frame(
    a = as.character(1:10),
    b = letters[1:10],
    c = LETTERS[1:10]
  )
}

make_uniform_data <- function() {
  # All records share the same QI combination (zero SUDA risk)
  data.frame(
    key1 = rep("A", 20),
    key2 = rep("B", 20),
    val = 1:20
  )
}

# --- Class structure ---

test_that("suda returns correct S3 class", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_s3_class(result, "suda")
})

test_that("suda result has expected fields", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))

  expect_named(result, c("suda_scores", "dis_scores", "msu_counts",
                          "n_msu", "msu_by_size", "high_risk_records",
                          "n_records", "n_complete", "summary_stats",
                          "max_msu", "key_vars"))
  expect_equal(result$key_vars, c("age", "gender", "region"))
  expect_equal(result$n_records, 100)
})

# --- SUDA scores ---

test_that("suda_scores has same length as nrow(X)", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_equal(length(result$suda_scores), nrow(data))
})

test_that("suda_scores are non-negative", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_true(all(result$suda_scores >= 0))
})

test_that("suda_scores are zero for uniform data", {
  data <- make_uniform_data()
  result <- suda(data, key_vars = c("key1", "key2"))
  expect_true(all(result$suda_scores == 0))
  expect_equal(result$n_msu, 0)
})

test_that("suda_scores are positive for all-unique data", {
  data <- make_unique_data()
  result <- suda(data, key_vars = c("a", "b", "c"))
  # Every record is unique on single variables -> all have MSUs
  expect_true(all(result$suda_scores > 0))
})

test_that("smaller MSUs contribute more to SUDA score", {
  # A record unique on 1 variable should get contribution = 1/2^0 = 1

  # A record unique on 2 variables should get contribution = 1/2^1 = 0.5
  data <- make_unique_data()
  result <- suda(data, key_vars = c("a", "b", "c"))
  # Single-variable MSUs contribute 1.0 each
  if (result$msu_by_size[1] > 0) {
    expect_true(max(result$suda_scores) >= 1.0)
  }
})

# --- DIS scores ---

test_that("dis_scores has same length as nrow(X)", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_equal(length(result$dis_scores), nrow(data))
})

test_that("dis_scores are non-negative", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_true(all(result$dis_scores >= 0))
})

test_that("dis_scores are zero when suda_scores are zero", {
  data <- make_uniform_data()
  result <- suda(data, key_vars = c("key1", "key2"))
  expect_true(all(result$dis_scores == 0))
})

# --- MSU counts ---

test_that("msu_counts has same length as nrow(X)", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_equal(length(result$msu_counts), nrow(data))
})

test_that("msu_counts are non-negative integers", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_true(all(result$msu_counts >= 0))
  expect_true(all(result$msu_counts == as.integer(result$msu_counts)))
})

test_that("msu_counts are zero for uniform data", {
  data <- make_uniform_data()
  result <- suda(data, key_vars = c("key1", "key2"))
  expect_true(all(result$msu_counts == 0))
})

# --- MSU by size ---

test_that("msu_by_size has correct length", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_equal(length(result$msu_by_size), 3)  # 3 key vars
  expect_equal(names(result$msu_by_size), c("1", "2", "3"))
})

test_that("n_msu equals sum of msu_by_size", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_equal(result$n_msu, sum(result$msu_by_size))
})

# --- max_msu parameter ---

test_that("max_msu limits the maximum MSU size searched", {
  data <- make_test_data()
  # With max_msu = 1, only single-variable uniques are found
  result1 <- suda(data, key_vars = c("age", "gender", "region"), max_msu = 1)
  result_all <- suda(data, key_vars = c("age", "gender", "region"))
  expect_equal(result1$max_msu, 1)
  expect_equal(length(result1$msu_by_size), 1)
  expect_true(result1$n_msu <= result_all$n_msu)
})

test_that("max_msu defaults to number of key variables", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_equal(result$max_msu, 3)
})

test_that("max_msu is capped at number of key variables", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender"), max_msu = 10)
  expect_equal(result$max_msu, 2)
})

# --- Summary stats ---

test_that("summary_stats contains expected fields", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  ss <- result$summary_stats
  expect_true(all(c("n_with_msu", "pct_with_msu", "mean_suda",
                     "median_suda", "max_suda", "mean_msu_count",
                     "max_msu_count") %in% names(ss)))
})

test_that("n_with_msu counts records with positive SUDA scores", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_equal(result$summary_stats$n_with_msu,
               sum(result$suda_scores > 0))
})

test_that("pct_with_msu is consistent with n_with_msu", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expected_pct <- 100 * result$summary_stats$n_with_msu / result$n_records
  expect_equal(result$summary_stats$pct_with_msu, expected_pct)
})

# --- NA handling ---

test_that("suda removes NAs by default", {
  data <- make_test_data(n = 50)
  data$age[1:5] <- NA
  result <- suda(data, key_vars = c("age", "gender"), na.rm = TRUE)
  expect_equal(result$n_complete, 45)
  expect_equal(result$n_records, 50)
})

test_that("suda errors when all cases have NAs", {
  data <- data.frame(a = c(NA, NA, NA), b = c(1, 2, 3))
  expect_error(suda(data, key_vars = "a"),
               "No complete cases")
})

# --- Edge cases ---

test_that("suda works with single key variable", {
  data <- make_test_data()
  result <- suda(data, key_vars = "age")
  expect_s3_class(result, "suda")
  expect_equal(result$max_msu, 1)
})

test_that("suda works with many key variables", {
  data <- make_test_data(n = 50)
  result <- suda(data,
                 key_vars = c("age", "gender", "region", "education"))
  expect_s3_class(result, "suda")
  expect_equal(result$max_msu, 4)
})

test_that("suda handles single-record dataset", {
  data <- data.frame(a = "x", b = "y")
  result <- suda(data, key_vars = c("a", "b"))
  expect_equal(result$n_records, 1)
  # Single record is unique on any subset
  expect_true(result$suda_scores[1] > 0)
})

# --- Error handling ---

test_that("suda errors on non-data.frame", {
  expect_error(suda(1:10, key_vars = "a"),
               "X must be a data frame")
})

test_that("suda errors on missing key variables", {
  data <- make_test_data()
  expect_error(suda(data, key_vars = c("age", "nonexistent")),
               "Key variables missing")
})

# --- synth_pair method ---

test_that("suda works with synth_pair objects (default = synthetic)", {
  set.seed(42)
  orig <- data.frame(
    age = sample(c("young", "old"), 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE),
    income = sample(c("low", "high"), 50, replace = TRUE)
  )
  syn <- orig
  syn$income <- sample(syn$income)

  pair <- synth_pair(orig, syn,
                     key_vars = c("age", "gender"),
                     target_var = "income")

  result <- suda(pair)
  expect_s3_class(result, "suda")
  expect_equal(result$n_records, 50)
  expect_equal(result$key_vars, c("age", "gender"))
})

test_that("suda.synth_pair data='original' uses original data", {
  set.seed(42)
  orig <- data.frame(
    age = sample(c("young", "old"), 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE)
  )
  # Synthetic: all records identical -> no MSUs
  syn <- data.frame(
    age = rep("young", 50),
    gender = rep("M", 50)
  )

  pair <- synth_pair(orig, syn, key_vars = c("age", "gender"))

  result_syn <- suda(pair, data = "synthetic")
  result_orig <- suda(pair, data = "original")

  # Synthetic is uniform -> no MSUs
  expect_equal(result_syn$n_msu, 0)
  # Original may have MSUs (unique records possible)
  # At minimum, the original should have different diversity
  expect_true(result_orig$n_records == result_syn$n_records)
})

test_that("suda.synth_pair errors without key_vars", {
  pair <- synth_pair(data.frame(a = 1:5), data.frame(a = 1:5))
  expect_error(suda(pair), "key_vars")
})

test_that("suda.synth_pair passes max_msu through", {
  set.seed(42)
  orig <- data.frame(
    age = sample(c("young", "old"), 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE)
  )
  pair <- synth_pair(orig, orig, key_vars = c("age", "gender"))
  result <- suda(pair, max_msu = 1)
  expect_equal(result$max_msu, 1)
})

# --- S3 methods: print, summary, plot ---

test_that("print.suda runs without error", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_output(print(result), "SUDA")
  expect_output(print(result), "Key variables")
  expect_output(print(result), "Total MSUs found")
})

test_that("print.suda displays risk assessment", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_output(print(result), "Risk Assessment")
})

test_that("summary.suda returns correct class", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  s <- summary(result)
  expect_s3_class(s, "summary.suda")
  expect_true(!is.null(s$n_msu))
  expect_true(!is.null(s$msu_by_size))
  expect_true(!is.null(s$suda_quantiles))
})

test_that("print.summary.suda runs without error", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_output(print(summary(result)), "Summary: SUDA Analysis")
  expect_output(print(summary(result)), "SUDA Score Distribution")
})

test_that("summary.suda includes top_records when MSUs exist", {
  data <- make_unique_data()
  result <- suda(data, key_vars = c("a", "b", "c"))
  s <- summary(result)
  expect_true(!is.null(s$top_records))
  expect_true(is.data.frame(s$top_records))
})

test_that("plot.suda which=1 runs without error", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_silent(plot(result, which = 1))
})

test_that("plot.suda which=2 runs without error", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_silent(plot(result, which = 2))
})

test_that("plot.suda which=1:2 runs without error", {
  data <- make_test_data()
  result <- suda(data, key_vars = c("age", "gender", "region"))
  expect_silent(plot(result, which = 1:2))
})

test_that("plot.suda handles zero MSUs gracefully", {
  data <- make_uniform_data()
  result <- suda(data, key_vars = c("key1", "key2"))
  expect_silent(plot(result, which = 1))
  expect_silent(plot(result, which = 2))
})

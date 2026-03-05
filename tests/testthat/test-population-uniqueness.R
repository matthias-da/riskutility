# Tests for population_uniqueness (Population Uniqueness Risk)

# --- Setup: shared test data ---

make_test_data <- function(n = 200, seed = 123) {
  set.seed(seed)
  data.frame(
    age = sample(c("young", "middle", "old"), n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE),
    income = sample(c("low", "medium", "high"), n, replace = TRUE)
  )
}

make_high_risk_data <- function(n = 100, seed = 42) {
  # Data with many unique combinations (high granularity keys)
  set.seed(seed)
  data.frame(
    id1 = seq_len(n),
    id2 = sample(letters, n, replace = TRUE),
    cat = sample(c("A", "B"), n, replace = TRUE)
  )
}

make_low_risk_data <- function(n = 200, seed = 123) {
  # Data with very few unique combinations (low granularity keys)
  set.seed(seed)
  data.frame(
    binary1 = sample(c("A", "B"), n, replace = TRUE),
    binary2 = sample(c("X", "Y"), n, replace = TRUE),
    val = rnorm(n)
  )
}

# --- Class structure ---

test_that("population_uniqueness returns correct S3 class", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01)
  expect_s3_class(result, "population_uniqueness")
})

test_that("population_uniqueness result has expected fields", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.05,
                                   method = "zayatz")

  expect_named(result, c("risk_per_record", "global_risk",
                           "n_sample_uniques", "pct_sample_uniques",
                           "n_population_uniques_est", "freq_table",
                           "method", "key_vars", "sampling_fraction",
                           "n_records", "privacy_pass", "comparison"))
  expect_equal(result$method, "zayatz")
  expect_equal(result$key_vars, c("age", "gender"))
  expect_equal(result$sampling_fraction, 0.05)
})

test_that("result contains comparison when method='all'", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01,
                                   method = "all")
  expect_true(!is.null(result$comparison))
  expect_named(result$comparison, c("pitman", "zayatz", "snb"))
})

test_that("result has NULL comparison for single method", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01,
                                   method = "pitman")
  expect_null(result$comparison)
})

# --- Core functionality ---

test_that("risk_per_record has correct length", {
  data <- make_test_data(n = 100)
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01,
                                   method = "zayatz")
  expect_equal(length(result$risk_per_record), 100)
})

test_that("risk_per_record is zero for non-unique records", {
  data <- make_test_data(n = 200)
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01,
                                   method = "zayatz")

  # Get record frequencies
  key_sig <- apply(data[, c("age", "gender"), drop = FALSE], 1, paste, collapse = "|")
  freq_tab <- table(key_sig)
  record_freq <- as.integer(freq_tab[key_sig])

  # Non-uniques should have zero risk
  non_unique_idx <- which(record_freq > 1)
  if (length(non_unique_idx) > 0) {
    expect_true(all(result$risk_per_record[non_unique_idx] == 0))
  }
})

test_that("risk_per_record values are between 0 and 1", {
  data <- make_test_data()
  for (m in c("pitman", "zayatz", "snb")) {
    result <- population_uniqueness(data,
                                     key_vars = c("age", "gender", "region"),
                                     sampling_fraction = 0.05,
                                     method = m)
    expect_true(all(result$risk_per_record >= 0), info = m)
    expect_true(all(result$risk_per_record <= 1), info = m)
  }
})

test_that("global_risk is mean of risk_per_record", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01,
                                   method = "pitman")
  expect_equal(result$global_risk, mean(result$risk_per_record))
})

test_that("n_sample_uniques is correctly counted", {
  data <- make_test_data(n = 100, seed = 99)
  kv <- c("age", "gender", "region")
  result <- population_uniqueness(data, key_vars = kv,
                                   sampling_fraction = 0.01)

  key_sig <- apply(data[, kv, drop = FALSE], 1, paste, collapse = "|")
  expected_uniques <- sum(table(key_sig) == 1)
  expect_equal(result$n_sample_uniques, expected_uniques)
})

test_that("pct_sample_uniques is consistent with n_sample_uniques", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01)
  expected_pct <- 100 * result$n_sample_uniques / result$n_records
  expect_equal(result$pct_sample_uniques, expected_pct)
})

test_that("n_population_uniques_est is non-negative", {
  data <- make_test_data()
  for (m in c("pitman", "zayatz", "snb")) {
    result <- population_uniqueness(data,
                                     key_vars = c("age", "gender"),
                                     sampling_fraction = 0.01,
                                     method = m)
    expect_true(result$n_population_uniques_est >= 0, info = m)
  }
})

# --- Zayatz method ---

test_that("Zayatz gives constant risk for all sample uniques", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender", "region", "income"),
                                   sampling_fraction = 0.05,
                                   method = "zayatz")
  # All nonzero values should be the same
  nonzero <- result$risk_per_record[result$risk_per_record > 0]
  expect_true(length(nonzero) > 0)
  expect_equal(length(unique(nonzero)), 1)
})

test_that("Zayatz formula is correct", {
  pi_frac <- 0.05
  expected_p <- (1 - pi_frac)^(1 / pi_frac - 1)

  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender", "region", "income"),
                                   sampling_fraction = pi_frac,
                                   method = "zayatz")

  nonzero <- result$risk_per_record[result$risk_per_record > 0]
  expect_true(length(nonzero) > 0)
  expect_equal(nonzero[1], expected_p)
})

test_that("Zayatz P(pop unique) increases with sampling fraction", {
  # Use 4 key vars to ensure sample uniques exist
  data <- make_test_data()
  kv <- c("age", "gender", "region", "income")

  r1 <- population_uniqueness(data, key_vars = kv,
                                sampling_fraction = 0.01, method = "zayatz")
  r2 <- population_uniqueness(data, key_vars = kv,
                                sampling_fraction = 0.1, method = "zayatz")

  # Zayatz P = (1-pi)^(1/pi - 1) increases with pi:
  # with larger samples, sample uniques are more likely to be truly unique
  # Only compare if there are sample uniques
  if (r1$n_sample_uniques > 0) {
    nonzero_1 <- max(r1$risk_per_record)
    nonzero_2 <- max(r2$risk_per_record)
    expect_true(nonzero_2 > nonzero_1)
  }
})

# --- Pitman method ---

test_that("Pitman method runs without error", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01,
                                   method = "pitman")
  expect_s3_class(result, "population_uniqueness")
  expect_equal(result$method, "pitman")
})

test_that("Pitman gives risk between 0 and 1", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender", "region"),
                                   sampling_fraction = 0.05,
                                   method = "pitman")
  expect_true(all(result$risk_per_record >= 0))
  expect_true(all(result$risk_per_record <= 1))
})

test_that("Pitman handles data with no pairs (f2=0)", {
  # All QI combos are unique
  data <- data.frame(
    k1 = letters[1:20],
    k2 = LETTERS[1:20],
    val = 1:20
  )
  result <- population_uniqueness(data,
                                   key_vars = c("k1", "k2"),
                                   sampling_fraction = 0.01,
                                   method = "pitman")
  expect_s3_class(result, "population_uniqueness")
  expect_true(result$global_risk > 0)
})

# --- SNB method ---

test_that("SNB method runs without error", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01,
                                   method = "snb")
  expect_s3_class(result, "population_uniqueness")
  expect_equal(result$method, "snb")
})

test_that("SNB gives risk between 0 and 1", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender", "region"),
                                   sampling_fraction = 0.05,
                                   method = "snb")
  expect_true(all(result$risk_per_record >= 0))
  expect_true(all(result$risk_per_record <= 1))
})

test_that("SNB handles low-variance data", {
  data <- make_low_risk_data()
  result <- population_uniqueness(data,
                                   key_vars = c("binary1", "binary2"),
                                   sampling_fraction = 0.1,
                                   method = "snb")
  expect_s3_class(result, "population_uniqueness")
})

# --- Method "all" ---

test_that("method='all' returns results for all three methods", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01,
                                   method = "all")

  expect_true(!is.null(result$comparison))
  expect_length(result$comparison, 3)

  for (nm in c("pitman", "zayatz", "snb")) {
    expect_true("risk_per_record" %in% names(result$comparison[[nm]]),
                info = nm)
    expect_true("global_risk" %in% names(result$comparison[[nm]]),
                info = nm)
  }
})

test_that("method='all' primary result uses pitman", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender", "region"),
                                   sampling_fraction = 0.01,
                                   method = "all")
  pitman_risk <- result$comparison$pitman$global_risk
  expect_equal(result$global_risk, pitman_risk)
})

# --- Edge cases ---

test_that("works with a single key variable", {
  data <- make_test_data()
  result <- population_uniqueness(data, key_vars = "age",
                                   sampling_fraction = 0.01)
  expect_s3_class(result, "population_uniqueness")
})

test_that("works with many key variables (many uniques)", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender", "region", "income"),
                                   sampling_fraction = 0.01)
  expect_s3_class(result, "population_uniqueness")
  # With 4 key vars and 200 records, should have more uniques
  expect_true(result$n_sample_uniques > 0)
})

test_that("works with all records unique", {
  data <- make_high_risk_data()
  result <- population_uniqueness(data,
                                   key_vars = c("id1", "id2"),
                                   sampling_fraction = 0.01)
  expect_s3_class(result, "population_uniqueness")
  expect_equal(result$n_sample_uniques, nrow(data))
})

test_that("works with no sample uniques", {
  # Every QI combo appears at least twice
  data <- data.frame(
    key = rep(c("A", "B", "C"), each = 10),
    val = 1:30
  )
  result <- population_uniqueness(data, key_vars = "key",
                                   sampling_fraction = 0.01)
  expect_equal(result$n_sample_uniques, 0)
  expect_equal(result$global_risk, 0)
  expect_true(all(result$risk_per_record == 0))
})

test_that("n_records reflects data after NA removal", {
  data <- make_test_data(n = 100)
  data$age[1:10] <- NA
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01)
  expect_equal(result$n_records, 90)
})

test_that("privacy_pass is TRUE when global_risk <= 0.1", {
  # Low-risk data: few uniques with small sampling fraction
  data <- make_low_risk_data()
  result <- population_uniqueness(data,
                                   key_vars = c("binary1", "binary2"),
                                   sampling_fraction = 0.01)
  # With only 4 QI combos and 200 records, no uniques expected
  if (result$global_risk <= 0.1) {
    expect_true(result$privacy_pass)
  } else {
    expect_false(result$privacy_pass)
  }
})

test_that("privacy_pass is FALSE when global_risk > 0.1", {
  # High risk: all unique + high sampling fraction
  data <- data.frame(k = seq_len(50), v = 1:50)
  result <- population_uniqueness(data, key_vars = "k",
                                   sampling_fraction = 0.5,
                                   method = "zayatz")
  # All records are unique, high sampling fraction -> high risk
  expect_true(result$global_risk > 0.1)
  expect_false(result$privacy_pass)
})

# --- Error handling ---

test_that("errors on non-data.frame input", {
  expect_error(population_uniqueness(1:10, key_vars = "a"),
               "X must be a data frame")
})

test_that("errors on missing key variables", {
  data <- data.frame(a = 1:5, b = letters[1:5])
  expect_error(population_uniqueness(data, key_vars = "c"),
               "Key variables missing")
  expect_error(population_uniqueness(data, key_vars = c("a", "z")),
               "Key variables missing")
})

test_that("errors on invalid sampling_fraction", {
  data <- make_test_data()
  expect_error(population_uniqueness(data, key_vars = "age",
                                      sampling_fraction = 0),
               "sampling_fraction must be")
  expect_error(population_uniqueness(data, key_vars = "age",
                                      sampling_fraction = 1),
               "sampling_fraction must be")
  expect_error(population_uniqueness(data, key_vars = "age",
                                      sampling_fraction = -0.1),
               "sampling_fraction must be")
  expect_error(population_uniqueness(data, key_vars = "age",
                                      sampling_fraction = 1.5),
               "sampling_fraction must be")
  expect_error(population_uniqueness(data, key_vars = "age",
                                      sampling_fraction = "abc"),
               "sampling_fraction must be")
})

test_that("errors on invalid method", {
  data <- make_test_data()
  expect_error(population_uniqueness(data, key_vars = "age",
                                      method = "invalid"))
})

test_that("errors when no complete cases remain", {
  data <- data.frame(a = c(NA, NA, NA), b = c(NA, NA, NA))
  expect_error(population_uniqueness(data, key_vars = "a"),
               "No complete cases")
})

# --- synth_pair method ---

test_that("works with synth_pair objects", {
  set.seed(123)
  orig <- data.frame(
    age = sample(c("young", "old"), 100, replace = TRUE),
    gender = sample(c("M", "F"), 100, replace = TRUE),
    income = sample(c("low", "high"), 100, replace = TRUE)
  )
  syn <- orig
  syn$income <- sample(syn$income)

  pair <- synth_pair(orig, syn,
                     key_vars = c("age", "gender"),
                     target_var = "income")

  result <- population_uniqueness(pair, sampling_fraction = 0.01)
  expect_s3_class(result, "population_uniqueness")
  expect_equal(result$key_vars, c("age", "gender"))
  # Assesses original data
  expect_equal(result$n_records, 100)
})

test_that("synth_pair errors without key_vars", {
  pair <- synth_pair(data.frame(a = 1:5), data.frame(a = 1:5))
  expect_error(population_uniqueness(pair), "key_vars")
})

# --- S3 methods: print, summary, plot ---

test_that("print.population_uniqueness runs without error", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01)
  expect_output(print(result), "Population Uniqueness Risk Assessment")
  expect_output(print(result), "Sampling fraction")
})

test_that("print works for method='all'", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01,
                                   method = "all")
  expect_output(print(result), "Method Comparison")
  expect_output(print(result), "pitman")
  expect_output(print(result), "zayatz")
  expect_output(print(result), "snb")
})

test_that("print works for single method", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01,
                                   method = "zayatz")
  expect_output(print(result), "Risk Estimate")
})

test_that("summary.population_uniqueness returns correct class", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01)
  s <- summary(result)
  expect_s3_class(s, "summary.population_uniqueness")
  expect_true("freq_distribution" %in% names(s))
  expect_true("global_risk" %in% names(s))
})

test_that("print.summary.population_uniqueness runs without error", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01)
  expect_output(print(summary(result)),
                "Summary: Population Uniqueness Risk Assessment")
})

test_that("summary print includes method comparison for method='all'", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender"),
                                   sampling_fraction = 0.01,
                                   method = "all")
  expect_output(print(summary(result)), "Method Comparison")
})

test_that("plot which=1 runs without error", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender", "region"),
                                   sampling_fraction = 0.05)
  expect_silent(plot(result, which = 1))
})

test_that("plot which=2 runs without error for method='all'", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender", "region"),
                                   sampling_fraction = 0.05,
                                   method = "all")
  expect_silent(plot(result, which = 2))
})

test_that("plot which=2 runs without error for single method", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender", "region"),
                                   sampling_fraction = 0.05,
                                   method = "zayatz")
  expect_silent(plot(result, which = 2))
})

test_that("plot which=1:2 runs without error", {
  data <- make_test_data()
  result <- population_uniqueness(data,
                                   key_vars = c("age", "gender", "region"),
                                   sampling_fraction = 0.05,
                                   method = "all")
  expect_silent(plot(result, which = 1:2))
})

test_that("plot which=1 handles no sample uniques", {
  data <- data.frame(
    key = rep(c("A", "B"), each = 20),
    val = 1:40
  )
  result <- population_uniqueness(data, key_vars = "key",
                                   sampling_fraction = 0.01)
  expect_silent(plot(result, which = 1))
})

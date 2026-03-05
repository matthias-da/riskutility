# Tests for attacker_risk (Prosecutor/Journalist/Marketer Models)

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

make_unique_data <- function() {
  # Every record has a unique QI combination
  data.frame(
    id = as.character(1:10),
    val = letters[1:10]
  )
}

make_uniform_data <- function() {
  # All records share the same QI combination
  data.frame(
    key = rep("A", 20),
    val = 1:20
  )
}

# --- Class structure ---

test_that("attacker_risk returns correct S3 class", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  expect_s3_class(result, "attacker_risk")
})

test_that("attacker_risk result has expected fields (model = 'all')", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"),
                          sampling_fraction = 0.05)

  expect_named(result, c("risk_per_record", "global_risk", "n_uniques",
                          "pct_uniques", "freq_table", "model", "key_vars",
                          "sampling_fraction", "n_records", "privacy_pass"))
  expect_equal(result$model, "all")
  expect_equal(result$key_vars, c("age", "gender"))
  expect_equal(result$sampling_fraction, 0.05)
  expect_equal(result$n_records, 200)
})

test_that("risk_per_record is a data.frame with correct columns", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  rpr <- result$risk_per_record
  expect_true(is.data.frame(rpr))
  expect_true("freq" %in% names(rpr))
  expect_true("prosecutor" %in% names(rpr))
  expect_true("journalist" %in% names(rpr))
  expect_equal(nrow(rpr), 200)
})

test_that("global_risk has all three models when model = 'all'", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  expect_true(all(c("prosecutor", "journalist", "marketer") %in%
                     names(result$global_risk)))
})

# --- Prosecutor model ---

test_that("prosecutor per-record risk equals 1/f_k", {
  data <- make_test_data(n = 50, seed = 42)
  result <- attacker_risk(data, key_vars = c("age", "gender"),
                          model = "prosecutor")
  expect_equal(result$risk_per_record$prosecutor,
               1 / result$risk_per_record$freq)
})

test_that("prosecutor global risk equals mean of per-record risks", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"),
                          model = "prosecutor")
  expect_equal(result$global_risk$prosecutor,
               mean(result$risk_per_record$prosecutor))
})

test_that("prosecutor risk is 1 for all-unique data", {
  data <- make_unique_data()
  result <- attacker_risk(data, key_vars = "id", model = "prosecutor")
  expect_equal(result$global_risk$prosecutor, 1)
  expect_true(all(result$risk_per_record$prosecutor == 1))
})

test_that("prosecutor risk is 1/n for uniform data", {
  data <- make_uniform_data()
  result <- attacker_risk(data, key_vars = "key", model = "prosecutor")
  expect_equal(result$global_risk$prosecutor, 1 / 20)
  expect_true(all(result$risk_per_record$prosecutor == 1 / 20))
})

test_that("prosecutor model only computes prosecutor risk", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"),
                          model = "prosecutor")
  expect_true("prosecutor" %in% names(result$global_risk))
  expect_false("journalist" %in% names(result$global_risk))
  expect_false("marketer" %in% names(result$global_risk))
  expect_true("prosecutor" %in% names(result$risk_per_record))
  expect_false("journalist" %in% names(result$risk_per_record))
})

# --- Journalist model ---

test_that("journalist per-record risk equals sampling_fraction/f_k", {
  data <- make_test_data(n = 50, seed = 42)
  sf <- 0.05
  result <- attacker_risk(data, key_vars = c("age", "gender"),
                          model = "journalist", sampling_fraction = sf)
  expected_risk <- sf / result$risk_per_record$freq
  expect_equal(result$risk_per_record$journalist, expected_risk)
})

test_that("journalist global risk equals mean of per-record journalist risks", {
  data <- make_test_data()
  sf <- 0.03
  result <- attacker_risk(data, key_vars = c("age", "gender"),
                          model = "journalist", sampling_fraction = sf)
  expect_equal(result$global_risk$journalist,
               mean(result$risk_per_record$journalist))
})

test_that("journalist risk equals sampling_fraction * prosecutor risk", {
  data <- make_test_data()
  sf <- 0.02
  result <- attacker_risk(data, key_vars = c("age", "gender"),
                          sampling_fraction = sf)
  expect_equal(result$global_risk$journalist,
               sf * result$global_risk$prosecutor,
               tolerance = 1e-12)
})

test_that("journalist risk decreases with smaller sampling fraction", {
  data <- make_test_data()
  r1 <- attacker_risk(data, key_vars = c("age", "gender"),
                       model = "journalist", sampling_fraction = 0.1)
  r2 <- attacker_risk(data, key_vars = c("age", "gender"),
                       model = "journalist", sampling_fraction = 0.01)
  expect_true(r1$global_risk$journalist > r2$global_risk$journalist)
})

test_that("journalist model only computes journalist risk", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"),
                          model = "journalist")
  expect_true("journalist" %in% names(result$global_risk))
  expect_false("prosecutor" %in% names(result$global_risk))
  expect_false("marketer" %in% names(result$global_risk))
})

# --- Marketer model ---

test_that("marketer global risk equals (1/n) * sum(1/f_k)", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"),
                          model = "marketer")
  f_k <- result$risk_per_record$freq
  expected <- (1 / result$n_records) * sum(1 / f_k)
  expect_equal(result$global_risk$marketer, expected)
})

test_that("marketer global risk equals prosecutor global risk", {
  # Both are mean(1/f_k) = (1/n) * sum(1/f_k)
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  expect_equal(result$global_risk$marketer,
               result$global_risk$prosecutor)
})

test_that("marketer model only computes marketer risk", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"),
                          model = "marketer")
  expect_true("marketer" %in% names(result$global_risk))
  expect_false("prosecutor" %in% names(result$global_risk))
  expect_false("journalist" %in% names(result$global_risk))
})

# --- Unique records ---

test_that("n_uniques is correct for all-unique data", {
  data <- make_unique_data()
  result <- attacker_risk(data, key_vars = "id")
  expect_equal(result$n_uniques, 10)
  expect_equal(result$pct_uniques, 1)
})

test_that("n_uniques is zero for uniform data", {
  data <- make_uniform_data()
  result <- attacker_risk(data, key_vars = "key")
  expect_equal(result$n_uniques, 0)
  expect_equal(result$pct_uniques, 0)
})

test_that("pct_uniques is fraction of unique ECs", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  ec_sizes <- as.numeric(result$freq_table)
  expected <- sum(ec_sizes == 1) / length(ec_sizes)
  expect_equal(result$pct_uniques, expected)
})

# --- Privacy pass ---

test_that("privacy_pass is TRUE when prosecutor risk <= 0.1", {
  # Large equivalence classes => low risk
  data <- make_uniform_data()
  result <- attacker_risk(data, key_vars = "key")
  expect_true(result$privacy_pass)
})

test_that("privacy_pass is FALSE when prosecutor risk > 0.1", {
  # All unique records => risk = 1
  data <- make_unique_data()
  result <- attacker_risk(data, key_vars = "id")
  expect_false(result$privacy_pass)
})

test_that("privacy_pass uses journalist when prosecutor not computed", {
  data <- make_unique_data()
  # journalist risk = sampling_fraction / 1 = 0.01 for unique records
  result <- attacker_risk(data, key_vars = "id", model = "journalist",
                          sampling_fraction = 0.01)
  # 0.01 <= 0.1, so should pass
  expect_true(result$privacy_pass)
})

# --- freq_table ---

test_that("freq_table sums to n_records", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  expect_equal(sum(result$freq_table), result$n_records)
})

test_that("freq_table matches per-record freqs", {
  data <- make_test_data(n = 50, seed = 1)
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  # Each record's freq should be the size of its EC
  for (i in 1:50) {
    sig <- paste(data$age[i], data$gender[i], sep = "|")
    expect_equal(unname(result$risk_per_record$freq[i]),
                 unname(as.numeric(result$freq_table[sig])))
  }
})

# --- NA handling ---

test_that("attacker_risk removes NAs by default", {
  data <- make_test_data(n = 50)
  data$age[1:5] <- NA
  result <- attacker_risk(data, key_vars = c("age", "gender"), na.rm = TRUE)
  expect_equal(result$n_records, 45)
})

test_that("attacker_risk errors when all cases have NAs", {
  data <- data.frame(a = c(NA, NA, NA), b = c(1, 2, 3))
  expect_error(attacker_risk(data, key_vars = "a"),
               "No complete cases")
})

# --- Edge cases ---

test_that("attacker_risk works with single key variable", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = "age")
  expect_s3_class(result, "attacker_risk")
  expect_true(result$n_records == 200)
})

test_that("attacker_risk works with many key variables", {
  data <- make_test_data()
  result <- attacker_risk(data,
                          key_vars = c("age", "gender", "region", "income"))
  expect_s3_class(result, "attacker_risk")
})

test_that("attacker_risk works with single-record dataset", {
  data <- data.frame(a = "x", b = "y")
  result <- attacker_risk(data, key_vars = c("a", "b"))
  expect_equal(result$n_records, 1)
  expect_equal(result$global_risk$prosecutor, 1)
})

test_that("attacker_risk handles sampling_fraction = 1", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"),
                          sampling_fraction = 1)
  # journalist risk should equal prosecutor risk when sf = 1
  expect_equal(result$global_risk$journalist,
               result$global_risk$prosecutor)
})

# --- Error handling ---

test_that("attacker_risk errors on non-data.frame", {
  expect_error(attacker_risk(1:10, key_vars = "a"),
               "X must be a data frame")
})

test_that("attacker_risk errors on missing key variables", {
  data <- make_test_data()
  expect_error(attacker_risk(data, key_vars = c("age", "nonexistent")),
               "Key variables missing")
})

test_that("attacker_risk errors on invalid sampling_fraction", {
  data <- make_test_data()
  expect_error(attacker_risk(data, key_vars = "age", sampling_fraction = 0),
               "sampling_fraction must be")
  expect_error(attacker_risk(data, key_vars = "age", sampling_fraction = -0.1),
               "sampling_fraction must be")
  expect_error(attacker_risk(data, key_vars = "age", sampling_fraction = 1.5),
               "sampling_fraction must be")
  expect_error(attacker_risk(data, key_vars = "age", sampling_fraction = "a"),
               "sampling_fraction must be")
})

test_that("attacker_risk errors on invalid model argument", {
  data <- make_test_data()
  expect_error(attacker_risk(data, key_vars = "age", model = "invalid"))
})

# --- synth_pair method ---

test_that("attacker_risk works with synth_pair objects", {
  set.seed(123)
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

  result <- attacker_risk(pair)
  expect_s3_class(result, "attacker_risk")
  # synth_pair method uses synthetic data
  expect_equal(result$n_records, 50)
  expect_equal(result$key_vars, c("age", "gender"))
})

test_that("attacker_risk.synth_pair errors without key_vars", {
  pair <- synth_pair(data.frame(a = 1:5), data.frame(a = 1:5))
  expect_error(attacker_risk(pair), "key_vars")
})

test_that("attacker_risk.synth_pair passes sampling_fraction through", {
  set.seed(123)
  orig <- data.frame(
    age = sample(c("young", "old"), 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE)
  )
  syn <- orig

  pair <- synth_pair(orig, syn, key_vars = c("age", "gender"))
  result <- attacker_risk(pair, sampling_fraction = 0.05)
  expect_equal(result$sampling_fraction, 0.05)
})

# --- S3 methods: print, summary, plot ---

test_that("print.attacker_risk runs without error (model = 'all')", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  expect_output(print(result), "Attacker Risk Models")
  expect_output(print(result), "Prosecutor")
  expect_output(print(result), "Journalist")
  expect_output(print(result), "Marketer")
})

test_that("print.attacker_risk runs for single model", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"),
                          model = "prosecutor")
  expect_output(print(result), "Prosecutor")
})

test_that("summary.attacker_risk returns correct class", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  s <- summary(result)
  expect_s3_class(s, "summary.attacker_risk")
  expect_true(!is.null(s$ec_size_distribution))
  expect_true(!is.null(s$risk_quantiles))
  expect_true(!is.null(s$n_high_risk))
})

test_that("print.summary.attacker_risk runs without error", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  expect_output(print(summary(result)), "Summary: Attacker Risk Models")
  expect_output(print(summary(result)), "Equivalence Class Size Distribution")
})

test_that("plot.attacker_risk which=1 runs without error", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  expect_silent(plot(result, which = 1))
})

test_that("plot.attacker_risk which=2 runs without error", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  expect_silent(plot(result, which = 2))
})

test_that("plot.attacker_risk which=1:2 runs without error", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  expect_silent(plot(result, which = 1:2))
})

test_that("plot.attacker_risk works with single model", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"),
                          model = "journalist")
  expect_silent(plot(result, which = 1))
  expect_silent(plot(result, which = 2))
})

# --- Consistency checks ---

test_that("risk values are bounded in [0, 1]", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  expect_true(all(result$risk_per_record$prosecutor >= 0))
  expect_true(all(result$risk_per_record$prosecutor <= 1))
  expect_true(all(result$risk_per_record$journalist >= 0))
  expect_true(all(result$risk_per_record$journalist <= 1))
  expect_true(result$global_risk$prosecutor >= 0)
  expect_true(result$global_risk$prosecutor <= 1)
  expect_true(result$global_risk$journalist >= 0)
  expect_true(result$global_risk$journalist <= 1)
  expect_true(result$global_risk$marketer >= 0)
  expect_true(result$global_risk$marketer <= 1)
})

test_that("journalist risk is always <= prosecutor risk", {
  data <- make_test_data()
  for (sf in c(0.01, 0.1, 0.5)) {
    result <- attacker_risk(data, key_vars = c("age", "gender"),
                            sampling_fraction = sf)
    expect_true(result$global_risk$journalist <=
                  result$global_risk$prosecutor + 1e-12)
  }
})

test_that("more key variables leads to higher or equal prosecutor risk", {
  data <- make_test_data(seed = 42)
  r1 <- attacker_risk(data, key_vars = "age", model = "prosecutor")
  r2 <- attacker_risk(data, key_vars = c("age", "gender"),
                       model = "prosecutor")
  r3 <- attacker_risk(data, key_vars = c("age", "gender", "region"),
                       model = "prosecutor")
  # More keys => smaller ECs => higher risk (generally)
  expect_true(r1$global_risk$prosecutor <= r2$global_risk$prosecutor + 1e-10)
  expect_true(r2$global_risk$prosecutor <= r3$global_risk$prosecutor + 1e-10)
})

test_that("consistent n_records across components", {
  data <- make_test_data()
  result <- attacker_risk(data, key_vars = c("age", "gender"))
  expect_equal(result$n_records, nrow(result$risk_per_record))
  expect_equal(result$n_records, sum(result$freq_table))
})

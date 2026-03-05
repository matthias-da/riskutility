# Tests for tcloseness (t-Closeness Assessment)

# --- Setup: shared test data ---

make_cat_data <- function(n = 200, seed = 123) {
  set.seed(seed)
  data.frame(
    age = sample(c("young", "middle", "old"), n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S"), n, replace = TRUE),
    disease = sample(c("healthy", "cold", "flu", "covid"), n, replace = TRUE,
                     prob = c(0.5, 0.2, 0.2, 0.1))
  )
}

make_num_data <- function(n = 200, seed = 123) {
  set.seed(seed)
  data.frame(
    age = sample(c("young", "middle", "old"), n, replace = TRUE),
    gender = sample(c("M", "F"), n, replace = TRUE),
    region = sample(c("N", "S"), n, replace = TRUE),
    salary = rnorm(n, mean = 50000, sd = 15000)
  )
}

# --- Class structure ---

test_that("tcloseness returns correct S3 class", {
  data <- make_cat_data()
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_s3_class(result, "tcloseness")
})

test_that("tcloseness result has expected fields", {
  data <- make_cat_data()
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease", t = 0.3)

  expect_named(result, c("t_achieved", "t_threshold", "satisfies_t",
                          "n_violating", "pct_violating", "n_records",
                          "n_ec", "per_ec", "sensitive_type",
                          "key_vars", "sensitive_var"))
  expect_equal(result$t_threshold, 0.3)
  expect_equal(result$key_vars, c("age", "gender"))
  expect_equal(result$sensitive_var, "disease")
})

test_that("per_ec data.frame has correct columns", {
  data <- make_cat_data()
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_true(is.data.frame(result$per_ec))
  expect_named(result$per_ec, c("key", "size", "emd", "violates_t"))
})

# --- Categorical sensitive variable ---

test_that("tcloseness works with categorical sensitive variable", {
  data <- make_cat_data()
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease", t = 0.3)

  expect_equal(result$sensitive_type, "categorical")
  expect_true(result$t_achieved >= 0)
  expect_true(result$t_achieved <= 1)
  expect_true(all(result$per_ec$emd >= 0))
  expect_true(all(result$per_ec$emd <= 1))
})

test_that("categorical EMD is zero when EC has same distribution as global", {
  # Build data where every EC has exactly the same category distribution
  data <- data.frame(
    key = rep(c("A", "B"), each = 100),
    sens = rep(c("x", "x", "y", "y"), 50)
  )
  result <- tcloseness(data, key_vars = "key", sensitive_var = "sens")
  expect_equal(result$t_achieved, 0)
  expect_true(result$satisfies_t)
})

test_that("categorical EMD is nonzero when EC distribution differs", {
  data <- data.frame(
    key = c(rep("A", 10), rep("B", 10)),
    sens = c(rep("x", 10), rep("y", 10))
  )
  result <- tcloseness(data, key_vars = "key", sensitive_var = "sens")
  expect_true(result$t_achieved > 0)
})

# --- Numeric sensitive variable ---

test_that("tcloseness works with numeric sensitive variable", {
  data <- make_num_data()
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "salary", t = 0.3)

  expect_equal(result$sensitive_type, "numeric")
  expect_true(result$t_achieved >= 0)
  expect_true(all(result$per_ec$emd >= 0))
})

test_that("numeric EMD is zero when all values are the same", {
  data <- data.frame(
    key = rep(c("A", "B"), each = 10),
    val = rep(42, 20)
  )
  result <- tcloseness(data, key_vars = "key", sensitive_var = "val")
  expect_equal(result$t_achieved, 0)
  expect_equal(result$sensitive_type, "numeric")
})

test_that("numeric EMD is zero when EC distribution equals global", {
  vals <- c(1, 2, 3, 4)
  data <- data.frame(
    key = rep(c("A", "B"), each = 4),
    val = rep(vals, 2)
  )
  result <- tcloseness(data, key_vars = "key", sensitive_var = "val")
  expect_equal(result$t_achieved, 0, tolerance = 1e-10)
})

# --- Violation detection ---

test_that("violations correctly detected when EMD exceeds threshold", {
  # One EC has very different distribution
  data <- data.frame(
    key = c(rep("A", 50), rep("B", 50)),
    sens = c(rep("x", 50), rep("y", 50))
  )
  result <- tcloseness(data, key_vars = "key", sensitive_var = "sens",
                       t = 0.1)
  expect_false(result$satisfies_t)
  expect_true(result$n_violating > 0)
  expect_equal(result$n_violating, 100)  # all records in violating ECs
})

test_that("no violations when t is very large", {
  data <- make_cat_data()
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease", t = 1.0)
  expect_true(result$satisfies_t)
  expect_equal(result$n_violating, 0)
})

test_that("pct_violating is consistent with n_violating", {
  data <- make_cat_data()
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease", t = 0.1)
  expected_pct <- 100 * result$n_violating / result$n_records
  expect_equal(result$pct_violating, expected_pct)
})

# --- Edge cases ---

test_that("tcloseness works with a single key variable", {
  data <- make_cat_data()
  result <- tcloseness(data, key_vars = "age",
                       sensitive_var = "disease")
  expect_s3_class(result, "tcloseness")
  expect_true(result$n_ec > 0)
})

test_that("tcloseness works with many key variables (small ECs)", {
  data <- make_cat_data()
  result <- tcloseness(data, key_vars = c("age", "gender", "region"),
                       sensitive_var = "disease")
  expect_s3_class(result, "tcloseness")
})

test_that("tcloseness handles single-record ECs", {
  set.seed(42)
  data <- data.frame(
    id = 1:10,
    sens = sample(c("a", "b"), 10, replace = TRUE)
  )
  result <- tcloseness(data, key_vars = "id", sensitive_var = "sens")
  expect_s3_class(result, "tcloseness")
  expect_equal(result$n_ec, 10)
})

test_that("tcloseness handles single EC (all same key)", {
  data <- data.frame(
    key = rep("A", 20),
    sens = sample(c("x", "y"), 20, replace = TRUE)
  )
  result <- tcloseness(data, key_vars = "key", sensitive_var = "sens")
  expect_equal(result$n_ec, 1)
  expect_equal(result$t_achieved, 0)
})

test_that("n_records reflects data after NA removal", {
  data <- make_cat_data(n = 50)
  data$disease[1:5] <- NA
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease", na.rm = TRUE)
  expect_equal(result$n_records, 45)
})

# --- Error handling ---

test_that("tcloseness errors on non-data.frame", {
  expect_error(tcloseness(1:10, key_vars = "a", sensitive_var = "b"),
               "X must be a data frame")
})

test_that("tcloseness errors on missing variables", {
  data <- data.frame(a = 1:5, b = letters[1:5])
  expect_error(tcloseness(data, key_vars = "c", sensitive_var = "b"),
               "Variables missing")
  expect_error(tcloseness(data, key_vars = "a", sensitive_var = "z"),
               "Variables missing")
})

test_that("tcloseness errors on invalid t parameter", {
  data <- make_cat_data()
  expect_error(tcloseness(data, key_vars = "age", sensitive_var = "disease",
                          t = -0.1), "t must be")
  expect_error(tcloseness(data, key_vars = "age", sensitive_var = "disease",
                          t = 1.5), "t must be")
  expect_error(tcloseness(data, key_vars = "age", sensitive_var = "disease",
                          t = "abc"), "t must be")
})

test_that("tcloseness errors when no complete cases remain", {
  data <- data.frame(a = c(NA, NA), b = c(NA, NA))
  expect_error(tcloseness(data, key_vars = "a", sensitive_var = "b"),
               "No complete cases")
})

# --- synth_pair method ---

test_that("tcloseness works with synth_pair objects", {
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

  result <- tcloseness(pair)
  expect_s3_class(result, "tcloseness")
  expect_equal(result$sensitive_var, "income")
})

test_that("tcloseness.synth_pair errors without key_vars", {
  pair <- synth_pair(data.frame(a = 1:5), data.frame(a = 1:5))
  expect_error(tcloseness(pair), "key_vars")
})

test_that("tcloseness.synth_pair errors without target_var", {
  pair <- synth_pair(data.frame(a = 1:5, b = 1:5),
                     data.frame(a = 1:5, b = 1:5),
                     key_vars = "a")
  expect_error(tcloseness(pair), "target_var")
})

# --- S3 methods: print, summary, plot ---

test_that("print.tcloseness runs without error", {
  data <- make_cat_data()
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_output(print(result), "t-Closeness Assessment")
  expect_output(print(result), "Maximum EMD")
})

test_that("summary.tcloseness returns correct class", {
  data <- make_cat_data()
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  s <- summary(result)
  expect_s3_class(s, "summary.tcloseness")
  expect_true(is.data.frame(s$emd_summary))
  expect_true(is.data.frame(s$worst_ec))
})

test_that("print.summary.tcloseness runs without error", {
  data <- make_cat_data()
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_output(print(summary(result)), "Summary: t-Closeness Assessment")
})

test_that("plot.tcloseness which=1 runs without error", {
  data <- make_num_data()
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "salary")
  expect_silent(plot(result, which = 1))
})

test_that("plot.tcloseness which=2 runs without error", {
  data <- make_cat_data()
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_silent(plot(result, which = 2))
})

test_that("plot.tcloseness which=1:2 runs without error", {
  data <- make_cat_data()
  result <- tcloseness(data, key_vars = c("age", "gender"),
                       sensitive_var = "disease")
  expect_silent(plot(result, which = 1:2))
})

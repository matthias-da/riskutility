# Tests for nnaa (Nearest-Neighbor Adversarial Accuracy)

# --- Class structure and fields ---

test_that("nnaa returns correct S3 class structure", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    age = rnorm(n, 40, 10),
    income = rnorm(n, 50000, 15000)
  )
  Y <- data.frame(
    age = rnorm(80, 40, 10),
    income = rnorm(80, 50000, 15000)
  )

  result <- nnaa(X, Y, method = "euclidean", seed = 42)

  expect_s3_class(result, "nnaa")
  expect_true("aa_train" %in% names(result))
  expect_true("aa_holdout" %in% names(result))
  expect_true("privacy_loss" %in% names(result))
  expect_true("privacy_pass" %in% names(result))
  expect_true("aa_train_left" %in% names(result))
  expect_true("aa_train_right" %in% names(result))
  expect_true("aa_holdout_left" %in% names(result))
  expect_true("aa_holdout_right" %in% names(result))
  expect_true("d_TS" %in% names(result))
  expect_true("d_TT" %in% names(result))
  expect_true("n_train" %in% names(result))
  expect_true("n_synthetic" %in% names(result))
  expect_true("n_holdout" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("vars" %in% names(result))
})

test_that("nnaa privacy_pass is logical", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- nnaa(X, Y, method = "euclidean", seed = 42)

  expect_type(result$privacy_pass, "logical")
})

# --- AA values for random/independent data ---

test_that("AA is approximately 0.5 for random independent data", {
  set.seed(42)
  n <- 300
  X <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )
  Y <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )

  result <- nnaa(X, Y, method = "euclidean", seed = 123)

  # AA should be around 0.5 for data from same distribution
  expect_true(result$aa_train > 0.3 && result$aa_train < 0.7)
  expect_true(result$aa_holdout > 0.3 && result$aa_holdout < 0.7)
})

test_that("AA components are between 0 and 1", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  result <- nnaa(X, Y, method = "euclidean", seed = 42)

  expect_true(result$aa_train >= 0 && result$aa_train <= 1)
  expect_true(result$aa_train_left >= 0 && result$aa_train_left <= 1)
  expect_true(result$aa_train_right >= 0 && result$aa_train_right <= 1)
  expect_true(result$aa_holdout >= 0 && result$aa_holdout <= 1)
  expect_true(result$aa_holdout_left >= 0 && result$aa_holdout_left <= 1)
  expect_true(result$aa_holdout_right >= 0 && result$aa_holdout_right <= 1)
})

test_that("AA is mean of left and right components", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  result <- nnaa(X, Y, method = "euclidean", seed = 42)

  expect_equal(result$aa_train,
               0.5 * (result$aa_train_left + result$aa_train_right))
  expect_equal(result$aa_holdout,
               0.5 * (result$aa_holdout_left + result$aa_holdout_right))
})

# --- Memorization detection ---

test_that("AA detects copied/memorized data", {
  set.seed(123)
  n <- 100
  X <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  # Y is an exact copy of X (severe memorization)
  Y <- X[sample(nrow(X), n, replace = TRUE), ]

  result <- nnaa(X, Y, method = "euclidean", seed = 42)

  # With copied data, AA(train) should be low (< 0.5)
  # because synthetic is closer to training than to itself

  expect_true(result$aa_train < 0.5)
})

test_that("privacy loss is positive for memorized data", {
  set.seed(123)
  n <- 100
  train <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  holdout <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  # Synthetic is exact copy of training
  Y <- train[sample(nrow(train), n, replace = TRUE), ]

  result <- nnaa(train, Y, holdout = holdout, method = "euclidean")

  # Privacy loss should be positive (holdout AA > train AA)
  expect_true(result$privacy_loss > 0)
})

test_that("privacy loss is near zero for independent data", {
  set.seed(42)
  n <- 200
  train <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  holdout <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  Y <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  result <- nnaa(train, Y, holdout = holdout, method = "euclidean")

  # Privacy loss should be close to 0 for independent data
  expect_true(abs(result$privacy_loss) < 0.15)
  expect_true(result$privacy_pass)
})

test_that("privacy_loss equals aa_holdout minus aa_train", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  result <- nnaa(X, Y, method = "euclidean", seed = 42)

  expect_equal(result$privacy_loss, result$aa_holdout - result$aa_train)
})

# --- Distance methods ---

test_that("nnaa works with gower method", {
  set.seed(123)
  n <- 60
  X <- data.frame(
    age = rnorm(n, 40, 10),
    gender = sample(c("M", "F"), n, replace = TRUE),
    income = rnorm(n, 50000, 15000)
  )
  Y <- data.frame(
    age = rnorm(40, 40, 10),
    gender = sample(c("M", "F"), 40, replace = TRUE),
    income = rnorm(40, 50000, 15000)
  )

  result <- nnaa(X, Y, method = "gower", seed = 42)

  expect_s3_class(result, "nnaa")
  expect_equal(result$method, "gower")
  expect_true(result$aa_train >= 0 && result$aa_train <= 1)
})

test_that("nnaa works with euclidean method", {
  set.seed(123)
  n <- 60
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- nnaa(X, Y, method = "euclidean", seed = 42)

  expect_s3_class(result, "nnaa")
  expect_equal(result$method, "euclidean")
})

test_that("nnaa euclidean errors with non-numeric data", {
  X <- data.frame(x = c("a", "b", "c", "d", "e"))
  Y <- data.frame(x = c("a", "b", "c"))

  expect_error(nnaa(X, Y, method = "euclidean"), "numeric")
})

# --- Holdout handling ---

test_that("nnaa works with explicit holdout", {
  set.seed(123)
  n <- 60
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))
  holdout <- data.frame(x1 = rnorm(30), x2 = rnorm(30))

  result <- nnaa(X, Y, holdout = holdout, method = "euclidean")

  expect_s3_class(result, "nnaa")
  # With explicit holdout, all of X is used as training
  expect_equal(result$n_train, nrow(X))
  expect_equal(result$n_holdout, nrow(holdout))
})

test_that("nnaa works with holdout_fraction", {
  set.seed(123)
  n <- 100
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(50), x2 = rnorm(50))

  result <- nnaa(X, Y, holdout_fraction = 0.3, method = "euclidean", seed = 42)

  expect_s3_class(result, "nnaa")
  expect_equal(result$n_holdout, floor(n * 0.3))
  expect_equal(result$n_train, n - floor(n * 0.3))
})

test_that("nnaa seed produces reproducible results", {
  set.seed(99)
  X <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  Y <- data.frame(x1 = rnorm(60), x2 = rnorm(60))

  r1 <- nnaa(X, Y, method = "euclidean", seed = 42)
  r2 <- nnaa(X, Y, method = "euclidean", seed = 42)

  expect_equal(r1$aa_train, r2$aa_train)
  expect_equal(r1$aa_holdout, r2$aa_holdout)
  expect_equal(r1$privacy_loss, r2$privacy_loss)
})

# --- synth_pair method ---

test_that("nnaa.synth_pair works", {
  set.seed(123)
  n <- 60
  original <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  synthetic <- data.frame(
    x1 = rnorm(40),
    x2 = rnorm(40)
  )

  pair <- synth_pair(original = original, synthetic = synthetic)
  result <- nnaa(pair, method = "euclidean", seed = 42)

  expect_s3_class(result, "nnaa")
  expect_true(!is.null(result$aa_train))
})

test_that("nnaa.synth_pair uses holdout when available", {
  set.seed(123)
  original <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  synthetic <- data.frame(x1 = rnorm(40), x2 = rnorm(40))
  holdout <- data.frame(x1 = rnorm(30), x2 = rnorm(30))

  pair <- synth_pair(original = original, synthetic = synthetic, holdout = holdout)
  result <- nnaa(pair, method = "euclidean")

  expect_equal(result$n_train, 60)
  expect_equal(result$n_holdout, 30)
})

# --- Input validation ---

test_that("nnaa input validation works", {
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  Y <- data.frame(x1 = rnorm(10), x2 = rnorm(10))

  expect_error(nnaa(1:10, Y), "X must be a data frame")
  expect_error(nnaa(X, 1:10), "Y must be a data frame")
  expect_error(nnaa(X, Y, holdout = "not a df"), "holdout must be a data frame")
})

test_that("nnaa errors when no common variables", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(b = 1:10)

  expect_error(nnaa(X, Y), "No common variables")
})

test_that("nnaa errors with mismatched variable types", {
  X <- data.frame(x = 1:10)
  Y <- data.frame(x = letters[1:10])

  expect_error(nnaa(X, Y), "different class")
})

test_that("nnaa errors with missing variables", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10)

  expect_error(nnaa(X, Y, vars = c("a", "b")), "Variables missing in X")
})

test_that("nnaa handles NA values", {
  set.seed(123)
  X <- data.frame(x1 = c(rnorm(18), NA, NA), x2 = rnorm(20))
  Y <- data.frame(x1 = c(rnorm(9), NA), x2 = rnorm(10))

  result <- nnaa(X, Y, method = "euclidean", seed = 42, na.rm = TRUE)

  expect_s3_class(result, "nnaa")
})

# --- NN distance vectors ---

test_that("d_TS and d_TT have correct length", {
  set.seed(123)
  n <- 60
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- nnaa(X, Y, method = "euclidean", seed = 42)

  expect_equal(length(result$d_TS), result$n_train)
  expect_equal(length(result$d_TT), result$n_train)
  expect_equal(length(result$d_ST), result$n_synthetic)
  expect_equal(length(result$d_SS), result$n_synthetic)
})

test_that("all distances are non-negative", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- nnaa(X, Y, method = "euclidean", seed = 42)

  expect_true(all(result$d_TS >= 0))
  expect_true(all(result$d_TT >= 0))
  expect_true(all(result$d_ST >= 0))
  expect_true(all(result$d_SS >= 0))
})

# --- print, summary, plot ---

test_that("nnaa print method works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- nnaa(X, Y, method = "euclidean", seed = 42)

  expect_output(print(result), "Nearest-Neighbor Adversarial Accuracy")
  expect_output(print(result), "Adversarial Accuracy")
  expect_output(print(result), "Privacy Assessment")
  expect_output(print(result), "Privacy Loss")
})

test_that("nnaa summary returns correct class", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- nnaa(X, Y, method = "euclidean", seed = 42)
  s <- summary(result)

  expect_s3_class(s, "summary.nnaa")
  expect_true(!is.null(s$aa_train))
  expect_true(!is.null(s$aa_holdout))
  expect_true(!is.null(s$privacy_loss))
  expect_true(!is.null(s$quantiles_d_TS))
  expect_true(!is.null(s$quantiles_d_TT))
})

test_that("nnaa print.summary works", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- nnaa(X, Y, method = "euclidean", seed = 42)
  s <- summary(result)

  expect_output(print(s), "Summary: Nearest-Neighbor Adversarial Accuracy")
  expect_output(print(s), "Key Metrics")
  expect_output(print(s), "AA Components")
})

test_that("nnaa plot method works without error", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40))

  result <- nnaa(X, Y, method = "euclidean", seed = 42)

  expect_no_error(plot(result, which = 1))
  expect_no_error(plot(result, which = 2))
  expect_no_error(plot(result, which = 1:2))
})

# --- Variable selection ---

test_that("nnaa uses specified vars", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  Y <- data.frame(x1 = rnorm(40), x2 = rnorm(40), x3 = rnorm(40))

  result <- nnaa(X, Y, vars = c("x1", "x2"), method = "euclidean", seed = 42)

  expect_equal(result$vars, c("x1", "x2"))
})

test_that("nnaa auto-selects common vars", {
  set.seed(123)
  X <- data.frame(a = rnorm(60), b = rnorm(60), c_only = rnorm(60))
  Y <- data.frame(a = rnorm(40), b = rnorm(40), d_only = rnorm(40))

  result <- nnaa(X, Y, method = "euclidean", seed = 42)

  expect_equal(sort(result$vars), c("a", "b"))
})

# --- Dataset sizes stored correctly ---

test_that("nnaa stores correct dataset sizes", {
  set.seed(123)
  X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
  Y <- data.frame(x1 = rnorm(80), x2 = rnorm(80))

  result <- nnaa(X, Y, holdout_fraction = 0.4, method = "euclidean", seed = 42)

  expect_equal(result$n_synthetic, 80)
  expect_equal(result$n_holdout, 40)
  expect_equal(result$n_train, 60)
})

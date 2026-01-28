# Tests for ims (Identical Match Share)

test_that("ims returns correct S3 class structure", {
  set.seed(123)
  X <- data.frame(
    age = sample(20:60, 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), 50, replace = TRUE)
  )
  Y <- data.frame(
    age = sample(20:60, 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE),
    region = sample(c("N", "S", "E", "W"), 50, replace = TRUE)
  )

  result <- ims(X, Y)

  expect_s3_class(result, "ims")
  expect_true("ims" %in% names(result))
  expect_true("ims_pct" %in% names(result))
  expect_true("n_identical" %in% names(result))
  expect_true("privacy_pass" %in% names(result))
})

test_that("ims values are in valid range", {
  set.seed(123)
  X <- data.frame(
    age = sample(20:60, 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE)
  )
  Y <- data.frame(
    age = sample(20:60, 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE)
  )

  result <- ims(X, Y)

  expect_true(result$ims >= 0 && result$ims <= 1)
  expect_true(result$ims_pct >= 0 && result$ims_pct <= 100)
  expect_true(result$n_identical >= 0)
  expect_true(result$n_identical <= result$n_synthetic)
})

test_that("ims detects identical data (100% match)", {
  X <- data.frame(
    age = c(25, 30, 35, 40),
    gender = c("M", "F", "M", "F")
  )
  Y <- X  # Identical

  result <- ims(X, Y)

  expect_equal(result$ims, 1)
  expect_equal(result$ims_pct, 100)
  expect_equal(result$n_identical, 4)
  expect_false(result$privacy_pass)
})

test_that("ims detects no matches for completely different data", {
  X <- data.frame(
    age = c(20, 25, 30),
    gender = c("M", "M", "M")
  )
  Y <- data.frame(
    age = c(50, 55, 60),
    gender = c("F", "F", "F")
  )

  result <- ims(X, Y)

  expect_equal(result$ims, 0)
  expect_equal(result$n_identical, 0)
  expect_true(result$privacy_pass)
})

test_that("ims handles partial copies", {
  X <- data.frame(
    age = c(25, 30, 35, 40, 45),
    gender = c("M", "F", "M", "F", "M")
  )
  # Y contains 2 copies from X and 2 new records
  Y <- data.frame(
    age = c(25, 30, 99, 99),
    gender = c("M", "F", "X", "Y")
  )

  result <- ims(X, Y)

  expect_equal(result$n_identical, 2)
  expect_equal(result$ims, 0.5)
})

test_that("ims handles non-dataframe input", {
  expect_error(ims(1:10, data.frame(a = 1:10)),
               "X must be a data frame")
  expect_error(ims(data.frame(a = 1:10), 1:10),
               "Y must be a data frame")
})

test_that("ims handles missing common variables", {
  X <- data.frame(a = 1:10, b = 1:10)
  Y <- data.frame(c = 1:10, d = 1:10)

  expect_error(ims(X, Y), "No common variables found")
})

test_that("ims print method works", {
  set.seed(123)
  X <- data.frame(
    age = sample(20:60, 30, replace = TRUE),
    gender = sample(c("M", "F"), 30, replace = TRUE)
  )
  Y <- data.frame(
    age = sample(20:60, 30, replace = TRUE),
    gender = sample(c("M", "F"), 30, replace = TRUE)
  )

  result <- ims(X, Y)

  expect_output(print(result), "Identical Match Share")
  expect_output(print(result), "IMS")
})

test_that("ims summary method works", {
  X <- data.frame(
    age = c(25, 30, 35, 40),
    gender = c("M", "F", "M", "F")
  )
  Y <- X[c(1, 2, 1, 2), ]  # Copies with repetition

  result <- ims(X, Y)
  summ <- summary(result)

  expect_s3_class(summ, "summary.ims")
  expect_true("match_distribution" %in% names(summ))
})

test_that("ims works with selected variables", {
  X <- data.frame(
    age = c(25, 30, 35),
    gender = c("M", "F", "M"),
    extra = c("A", "B", "C")
  )
  Y <- data.frame(
    age = c(25, 30, 99),
    gender = c("M", "F", "X"),
    extra = c("X", "Y", "Z")  # Different
  )

  # With all variables
  result_all <- ims(X, Y)
  expect_equal(result_all$n_identical, 0)

  # With selected variables only
  result_subset <- ims(X, Y, vars = c("age", "gender"))
  expect_equal(result_subset$n_identical, 2)
})

# Tests for repu function
test_that("repu returns correct S3 class", {
  set.seed(123)
  X <- data.frame(
    age = sample(20:60, 50, replace = TRUE),
    gender = sample(c("M", "F"), 50, replace = TRUE)
  )
  Y <- X[sample(nrow(X), 30, replace = TRUE), ]

  result <- repu(X, Y)

  expect_s3_class(result, "ims")
})

test_that("repu focuses on unique records", {
  # Create data with some duplicates and some uniques
  X <- data.frame(
    age = c(25, 25, 30, 35, 40),  # age 25 is duplicated
    gender = c("M", "M", "F", "M", "F")
  )
  # Y copies one duplicate (age=25,M) and one unique (age=30,F)
  Y <- data.frame(
    age = c(25, 30),
    gender = c("M", "F")
  )

  # Regular IMS counts both matches
  result_ims <- ims(X, Y)
  expect_equal(result_ims$n_identical, 2)

  # RepU only counts match to unique record (age=30,F)
  result_repu <- repu(X, Y, uniques_only = TRUE)
  expect_equal(result_repu$n_identical, 1)
})

test_that("repu with uniques_only=FALSE equals ims", {
  set.seed(123)
  X <- data.frame(
    age = sample(20:40, 30, replace = TRUE),
    gender = sample(c("M", "F"), 30, replace = TRUE)
  )
  Y <- X[sample(nrow(X), 20, replace = TRUE), ]

  result_ims <- ims(X, Y)
  result_repu <- repu(X, Y, uniques_only = FALSE)

  expect_equal(result_ims$ims, result_repu$ims)
  expect_equal(result_ims$n_identical, result_repu$n_identical)
})

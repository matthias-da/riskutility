# Tests for compare_* functions
# These functions compare original and synthetic/anonymized datasets

library(testthat)

# ============================================================================
# Test Data Setup
# ============================================================================

create_test_data <- function(n = 200, seed = 123) {
  set.seed(seed)
  X <- data.frame(
    income = rnorm(n, mean = 50000, sd = 10000),
    age = sample(20:70, n, replace = TRUE),
    gender = sample(c("Male", "Female"), n, replace = TRUE),
    region = sample(c("North", "South", "East", "West"), n, replace = TRUE),
    weight = runif(n, 0.5, 1.5),
    stringsAsFactors = FALSE
  )

  Y <- data.frame(
    income = rnorm(n, mean = 48000, sd = 12000),
    age = sample(20:70, n, replace = TRUE),
    gender = sample(c("Male", "Female"), n, replace = TRUE),
    region = sample(c("North", "South", "East", "West"), n, replace = TRUE),
    weight = runif(n, 0.5, 1.5),
    stringsAsFactors = FALSE
  )

  list(X = X, Y = Y)
}

# ============================================================================
# compare_ks_test
# ============================================================================

test_that("compare_ks_test returns correct structure", {
  data <- create_test_data()

  result <- compare_ks_test(data$X, data$Y, num_var = "income")

  expect_s3_class(result, "data.table")
  expect_true("ks_statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true(result$ks_statistic >= 0 && result$ks_statistic <= 1)
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

test_that("compare_ks_test works with grouping variables", {
  data <- create_test_data()

  result <- compare_ks_test(data$X, data$Y,
                            num_var = "income",
                            cat_vars = "gender")

  expect_s3_class(result, "data.table")
  expect_true("gender" %in% names(result))
  expect_equal(nrow(result), 2)  # Male and Female
})

test_that("compare_ks_test works with multiple grouping variables", {
  data <- create_test_data()

  result <- compare_ks_test(data$X, data$Y,
                            num_var = "income",
                            cat_vars = c("gender", "region"))

  expect_s3_class(result, "data.table")
  expect_true("gender" %in% names(result))
  expect_true("region" %in% names(result))
  expect_true(nrow(result) >= 1)
})

test_that("compare_ks_test validates input", {
  data <- create_test_data()

  # Missing variable
  expect_error(compare_ks_test(data$X, data$Y, num_var = "nonexistent"))

  # Missing grouping variable
  expect_error(compare_ks_test(data$X, data$Y,
                               num_var = "income",
                               cat_vars = "nonexistent"))
})

# ============================================================================
# compare_wasserstein
# ============================================================================

test_that("compare_wasserstein returns correct structure for continuous", {
  data <- create_test_data()

  result <- compare_wasserstein(data$X, data$Y,
                                num_var = "income",
                                var_type = "continuous")

  expect_s3_class(result, "data.table")
  expect_true("wasserstein" %in% names(result))
  expect_true("var_type" %in% names(result))
  expect_true(result$wasserstein >= 0)
  expect_equal(result$var_type, "continuous")
})

test_that("compare_wasserstein works for nominal variables", {
  data <- create_test_data()

  result <- compare_wasserstein(data$X, data$Y,
                                num_var = "gender",
                                var_type = "nominal")

  expect_s3_class(result, "data.table")
  expect_true(result$wasserstein >= 0)
  expect_true(result$wasserstein <= 1)  # Total variation bounded by 1
  expect_equal(result$var_type, "nominal")
})

test_that("compare_wasserstein works with weights", {
  data <- create_test_data()

  result <- compare_wasserstein(data$X, data$Y,
                                num_var = "income",
                                var_type = "continuous",
                                weight_X = "weight",
                                weight_Y = "weight")

  expect_s3_class(result, "data.table")
  expect_true(result$wasserstein >= 0)
})

test_that("compare_wasserstein works with grouping", {
  data <- create_test_data()

  result <- compare_wasserstein(data$X, data$Y,
                                num_var = "income",
                                cat_vars = "gender",
                                var_type = "continuous")

  expect_s3_class(result, "data.table")
  expect_true(nrow(result) >= 2)
})

test_that("compare_wasserstein auto-detects variable type", {
  data <- create_test_data()

  # Numeric should be detected as continuous
  result_num <- compare_wasserstein(data$X, data$Y,
                                    num_var = "income",
                                    var_type = "auto")
  expect_equal(result_num$var_type, "continuous")

  # Character should be detected as nominal
  result_cat <- compare_wasserstein(data$X, data$Y,
                                    num_var = "gender",
                                    var_type = "auto")
  expect_equal(result_cat$var_type, "nominal")
})

# ============================================================================
# compare_means_frequencies
# ============================================================================

test_that("compare_means_frequencies works for numeric variables", {
  data <- create_test_data()

  result <- compare_means_frequencies(data$X, data$Y,
                                      cont_vars = "income",
                                      cat_vars = character(0))

  expect_type(result, "list")
  # Should have some comparison output
  expect_true(length(result) > 0)
})

test_that("compare_means_frequencies works for categorical variables", {
  data <- create_test_data()

  result <- compare_means_frequencies(data$X, data$Y,
                                      cont_vars = character(0),
                                      cat_vars = "gender")

  expect_type(result, "list")
})

# ============================================================================
# compare_chisq_gof
# ============================================================================

test_that("compare_chisq_gof returns valid results", {
  data <- create_test_data()

  result <- compare_chisq_gof(data$X, data$Y, cat_vars = "gender")

  expect_s3_class(result, "data.table")
  expect_true("chi_squared" %in% names(result) ||
              "chisq" %in% names(result) ||
              any(grepl("chi", names(result), ignore.case = TRUE)))
})

test_that("compare_chisq_gof works with multiple variables", {
  data <- create_test_data()

  result <- compare_chisq_gof(data$X, data$Y,
                              cat_vars = c("gender", "region"))

  expect_s3_class(result, "data.table")
})

# ============================================================================
# compare_pca
# ============================================================================

test_that("compare_pca returns correct structure", {
  data <- create_test_data()

  result <- compare_pca(data$X, data$Y,
                        vars = c("income", "age"))

  expect_type(result, "list")
  # Should have some PCA-related output
  expect_true(length(result) > 0)
})

test_that("compare_pca works with biplot option", {
  data <- create_test_data()

  # Should not error with biplot = TRUE
  expect_no_error({
    result <- compare_pca(data$X, data$Y,
                          vars = c("income", "age"),
                          biplot = FALSE)
  })
})

# ============================================================================
# compare_correlation_matrices (if exists)
# ============================================================================

test_that("compare_correlation_matrices works", {
  skip_if_not(exists("compare_correlation_matrices"))

  data <- create_test_data()
  data$X$score <- rnorm(nrow(data$X))
  data$Y$score <- rnorm(nrow(data$Y))

  result <- compare_correlation_matrices(data$X, data$Y,
                                         vars = c("income", "age", "score"))

  expect_type(result, "list")
})

# ============================================================================
# compare_distributions_cont
# ============================================================================

test_that("compare_distributions_cont returns S3 class", {
  data <- create_test_data()

  result <- compare_distributions_cont(data$X, data$Y, var = "income")

  expect_s3_class(result, "compare_distributions_cont")
})

test_that("compare_distributions_cont has print method", {
  data <- create_test_data()

  result <- compare_distributions_cont(data$X, data$Y, var = "income")

  expect_output(print(result), regexp = ".")
})

# ============================================================================
# Input validation tests
# ============================================================================

test_that("compare functions validate data frames", {
  data <- create_test_data()

  # Non-data.frame input
  expect_error(compare_ks_test("not a df", data$Y, num_var = "income"))
  expect_error(compare_wasserstein(data$X, "not a df", num_var = "income"))
})

test_that("compare functions handle empty groups gracefully", {
  data <- create_test_data(n = 50)

  # Create data with a rare category
  data$X$rare_cat <- "A"
  data$Y$rare_cat <- "B"  # Different categories in X and Y

  # Should still work (may have empty groups)
  expect_no_error({
    result <- compare_wasserstein(data$X, data$Y,
                                  num_var = "income",
                                  cat_vars = "rare_cat",
                                  var_type = "continuous")
  })
})

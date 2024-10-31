test_that("compare function", {
  # Create sample data
  set.seed(123)
  X <- data.frame(
    age = sample(20:80, 100, replace = TRUE),
    income = rnorm(100, mean = 50000, sd = 10000),
    weights = runif(100, 0.5, 2)
  )
  Y <- data.frame(
    age = sample(20:80, 10000, replace = TRUE),
    income = rnorm(10000, mean = 50000, sd = 10000)
  )

  # Test the compare function
  result <- compare(X, Y, variables = c("age", "income"), weights = "weights")
  expect_is(result, "compare")
  expect_equal(names(result), c("formula", "ecdf", "kind"))
#  expect_equal(dim(result$ecdf), c(20000, 3))
  expect_equal(result$kind, "numeric")
})

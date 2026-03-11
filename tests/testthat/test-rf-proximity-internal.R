test_that(".rf_proximity returns correct structure", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(50), y = rnorm(50))
  d2 <- data.frame(x = rnorm(50), y = rnorm(50))
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 1)

  expect_type(res, "list")
  expect_s3_class(res$forest, "ranger")
  expect_true(is.matrix(res$terminal_nodes))
  expect_equal(nrow(res$terminal_nodes), 100)  # n1 + n2
  expect_equal(ncol(res$terminal_nodes), 50)   # n_trees
  expect_equal(res$n1, 50)
  expect_equal(res$n2, 50)
  expect_true(is.numeric(res$oob_error))
  expect_true(res$oob_error >= 0 && res$oob_error <= 1)
})

test_that(".rf_proximity returns importance when requested", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(50), y = rnorm(50))
  d2 <- data.frame(x = rnorm(50) + 2, y = rnorm(50))
  res_imp <- riskutility:::.rf_proximity(d1, d2, importance = TRUE,
                                          n_trees = 50, seed = 1)
  res_no  <- riskutility:::.rf_proximity(d1, d2, importance = FALSE,
                                          n_trees = 50, seed = 1)

  expect_true(is.numeric(res_imp$importance))
  expect_equal(length(res_imp$importance), 2)  # x, y
  expect_null(res_no$importance)
})

test_that(".rf_proximity validates n_trees >= 10", {
  skip_if_not_installed("ranger")
  d1 <- data.frame(x = 1:20)
  d2 <- data.frame(x = 21:40)
  expect_error(riskutility:::.rf_proximity(d1, d2, n_trees = 5),
               "n_trees")
})

test_that(".rf_proximity detects .rf_label collision", {
  skip_if_not_installed("ranger")
  d1 <- data.frame(x = 1:20, .rf_label = 1:20)
  d2 <- data.frame(x = 21:40, .rf_label = 21:40)
  expect_error(riskutility:::.rf_proximity(d1, d2, n_trees = 10),
               "rf_label")
})

test_that(".rf_proximity uses vars subset", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(30), y = rnorm(30), z = rnorm(30))
  d2 <- data.frame(x = rnorm(30), y = rnorm(30), z = rnorm(30))
  res <- riskutility:::.rf_proximity(d1, d2, vars = c("x", "y"),
                                      n_trees = 50, seed = 1)
  # importance should only have x, y (not z)
  expect_equal(length(res$importance), 2)
  expect_true(all(names(res$importance) %in% c("x", "y")))
})

test_that(".rf_proximity handles high-cardinality factors", {
  skip_if_not_installed("ranger")
  set.seed(1)
  # Create factor with > 53 levels to trigger order mode
  d1 <- data.frame(x = factor(sample(paste0("cat", 1:60), 100, TRUE)),
                   y = rnorm(100))
  d2 <- data.frame(x = factor(sample(paste0("cat", 1:60), 100, TRUE)),
                   y = rnorm(100))
  expect_message(
    res <- riskutility:::.rf_proximity(d1, d2, n_trees = 20, seed = 1),
    "order"
  )
  expect_s3_class(res$forest, "ranger")
})

test_that(".rf_proximity OOB NA handling for small data", {
  skip_if_not_installed("ranger")
  set.seed(1)
  # Very small dataset where OOB predictions may produce NAs
  d1 <- data.frame(x = 1:5, y = rnorm(5))
  d2 <- data.frame(x = 6:10, y = rnorm(5))
  # Should run without error (warning about NAs is acceptable)
  res <- suppressWarnings(
    riskutility:::.rf_proximity(d1, d2, n_trees = 10, seed = 1)
  )
  expect_false(any(is.na(res$terminal_nodes)))
})

test_that(".rf_proximity uses modifyList for user overrides", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(50), y = rnorm(50))
  d2 <- data.frame(x = rnorm(50), y = rnorm(50))
  # Override importance via ... should work without collision
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 1,
                                      importance = TRUE,
                                      min.node.size = 5)
  expect_s3_class(res$forest, "ranger")
})

test_that(".rf_proximity is reproducible with seed", {
  skip_if_not_installed("ranger")
  d1 <- data.frame(x = rnorm(50), y = rnorm(50))
  d2 <- data.frame(x = rnorm(50), y = rnorm(50))
  r1 <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 42)
  r2 <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 42)
  expect_identical(r1$terminal_nodes, r2$terminal_nodes)
})

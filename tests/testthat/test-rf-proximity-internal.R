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

# ── Task 3: .proximity_from_nodes ─────────────────────────────────────────────

test_that(".proximity_from_nodes returns correct dimensions", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(20), y = rnorm(20))
  d2 <- data.frame(x = rnorm(30), y = rnorm(30))
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 1)

  idx1 <- 1:20     # data1 indices
  idx2 <- 21:50    # data2 indices
  prox <- riskutility:::.proximity_from_nodes(res$terminal_nodes, idx1, idx2)

  expect_true(is.matrix(prox))
  expect_equal(nrow(prox), 30)  # length(idx2)
  expect_equal(ncol(prox), 20)  # length(idx1)
})

test_that(".proximity_from_nodes values are in [0, 1]", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(20))
  d2 <- data.frame(x = rnorm(20))
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 100, seed = 1)

  prox <- riskutility:::.proximity_from_nodes(res$terminal_nodes, 1:20, 21:40)
  expect_true(all(prox >= 0 & prox <= 1))
})

test_that(".proximity_from_nodes: identical data has high self-proximity", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d <- data.frame(x = rnorm(30), y = rnorm(30))
  # data2 = copy of data1 (memorized)
  res <- riskutility:::.rf_proximity(d, d, n_trees = 200, seed = 1)

  # Self-proximity (record i in d1 vs record i in d2)
  prox <- riskutility:::.proximity_from_nodes(res$terminal_nodes, 1:30, 31:60)
  self_prox <- diag(prox)
  other_prox <- prox[row(prox) != col(prox)]

  # Self-proximity should generally be higher than cross-proximity
  expect_true(mean(self_prox) > mean(other_prox))
})

test_that(".proximity_from_nodes is symmetric", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(15), y = rnorm(15))
  d2 <- data.frame(x = rnorm(15), y = rnorm(15))
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 1)

  # prox(idx1->idx2) should be transpose of prox(idx2->idx1)
  prox_ab <- riskutility:::.proximity_from_nodes(res$terminal_nodes, 1:15, 16:30)
  prox_ba <- riskutility:::.proximity_from_nodes(res$terminal_nodes, 16:30, 1:15)
  expect_equal(prox_ab, t(prox_ba))
})

test_that(".proximity_from_nodes: tie correction produces rational values", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(10))
  d2 <- data.frame(x = rnorm(10))
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 1)

  prox <- riskutility:::.proximity_from_nodes(res$terminal_nodes, 1:10, 11:20)
  # Values should be multiples of 1/n_trees (allow for floating-point tolerance)
  expect_true(all(abs(prox * 50 - round(prox * 50)) < 1e-9))
})

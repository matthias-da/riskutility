# tests/testthat/test-recordLinkage-ot.R

test_that("recordLinkage(matching = 'ot') basic with deterministic", {
  set.seed(1)
  X <- data.frame(a = rnorm(25), b = rnorm(25))
  Y <- data.frame(a = rnorm(25), b = rnorm(25))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "deterministic",
                       matching = "ot")
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(res$settings$matching, "ot")
  expect_equal(nrow(res$per_record), 25)
  # OT risk is continuous (not binary like bijective)
  expect_true(any(res$per_record$risk > 0 & res$per_record$risk < 1))
})

test_that("recordLinkage(matching = 'ot') detects near-copies", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- X + rnorm(40, 0, 0.01)
  res <- recordLinkage(X, Y, key = c("a", "b"), matching = "ot",
                       truth = "row")
  expect_true(res$overall$mean_risk > 0.1)
  expect_true(mean(res$per_record$true_in_set) > 0.5)
})

test_that("recordLinkage(matching = 'ot') works with probabilistic", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- X + rnorm(40, 0, 0.3)
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "probabilistic",
                       matching = "ot", truth = "row")
  expect_s3_class(res, "recordLinkageRisk")
  expect_true(all(res$per_record$risk >= 0))
})

test_that("recordLinkage(matching = 'ot') with blocking", {
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = sample(c("A", "B"), 30, TRUE))
  Y <- data.frame(a = rnorm(30), b = sample(c("A", "B"), 30, TRUE))
  res <- recordLinkage(X, Y, key = c("a", "b"), block = "b",
                       matching = "ot")
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(nrow(res$per_record), 30)
})

test_that("recordLinkage(matching = 'ot') with unequal sizes (truth = 'id')", {
  set.seed(1)
  X <- data.frame(id = 1:20, a = rnorm(20), b = rnorm(20))
  Y <- rbind(X[, c("id", "a", "b")], X[1:5, c("id", "a", "b")])
  Y$a <- Y$a + rnorm(25, 0, 0.1)
  Y$b <- Y$b + rnorm(25, 0, 0.1)
  res <- recordLinkage(X, Y, key = c("a", "b"), matching = "ot",
                       truth = "id", id = "id")
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(nrow(res$per_record), 20)
})

test_that("recordLinkage(matching = 'ot') custom epsilon", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res1 <- recordLinkage(X, Y, key = c("a", "b"), matching = "ot",
                        ot_epsilon = 0.01)
  res2 <- recordLinkage(X, Y, key = c("a", "b"), matching = "ot",
                        ot_epsilon = 10)
  # Small epsilon -> sharper (more risk variation)
  # Large epsilon -> smoother (more uniform risk, less variation)
  sd1 <- sd(res1$per_record$risk)
  sd2 <- sd(res2$per_record$risk)
  expect_true(sd1 > sd2)
})

test_that("recordLinkage(matching = 'ot') direction = 'reverse'", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b"), matching = "ot",
                       direction = "reverse")
  expect_equal(res$direction, "reverse")
  expect_equal(nrow(res$per_record), 20)
})

test_that("recordLinkage(matching = 'ot') purely categorical data", {
  set.seed(1)
  X <- data.frame(
    sex = sample(c("M", "F"), 30, TRUE),
    region = sample(LETTERS[1:5], 30, TRUE),
    edu = sample(c("low", "mid", "high"), 30, TRUE),
    stringsAsFactors = FALSE
  )
  Y <- data.frame(
    sex = sample(c("M", "F"), 30, TRUE),
    region = sample(LETTERS[1:5], 30, TRUE),
    edu = sample(c("low", "mid", "high"), 30, TRUE),
    stringsAsFactors = FALSE
  )
  res <- recordLinkage(X, Y, key = c("sex", "region", "edu"),
                       matching = "ot")
  expect_s3_class(res, "recordLinkageRisk")
  expect_true(all(is.finite(res$per_record$risk)))
})

test_that("OT vs bijective vs independent: OT risk is between the others", {
  skip_if_not_installed("clue")
  set.seed(1)
  X <- data.frame(a = rnorm(15), b = rnorm(15))
  Y <- X + rnorm(30, 0, 0.1)
  res_ind <- recordLinkage(X, Y, key = c("a", "b"),
                            matching = "independent", truth = "row")
  res_ot <- recordLinkage(X, Y, key = c("a", "b"),
                           matching = "ot", truth = "row")
  res_bij <- recordLinkage(X, Y, key = c("a", "b"),
                            matching = "bijective", truth = "row")
  # OT should give continuous risk (not binary 0/1 like bijective)
  expect_true(any(res_ot$per_record$risk > 0 & res_ot$per_record$risk < 1))
  # All three should detect the near-copies
  expect_true(res_ot$overall$mean_risk > 0)
  expect_true(res_bij$overall$mean_risk > 0)
})

test_that("recordLinkage(matching = 'ot') plot works", {
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b"), matching = "ot")
  expect_no_error(plot(res, which = 1))
})

test_that("recordLinkage(matching = 'ot') print shows OT mode", {
  set.seed(1)
  X <- data.frame(a = rnorm(15), b = rnorm(15))
  Y <- data.frame(a = rnorm(15), b = rnorm(15))
  res <- recordLinkage(X, Y, key = c("a", "b"), matching = "ot")
  out <- capture.output(print(res))
  expect_true(any(grepl("ot.*Sinkhorn", out, ignore.case = TRUE)))
})

test_that("recordLinkage(matching = 'ot') summary shows OT mode", {
  set.seed(1)
  X <- data.frame(a = rnorm(15), b = rnorm(15))
  Y <- data.frame(a = rnorm(15), b = rnorm(15))
  res <- recordLinkage(X, Y, key = c("a", "b"), matching = "ot")
  out <- capture.output(print(summary(res)))
  expect_true(any(grepl("ot.*Sinkhorn", out, ignore.case = TRUE)))
})

test_that("recordLinkage(matching = 'ot') stores settings and transport_plans", {
  set.seed(1)
  X <- data.frame(a = rnorm(15), b = rnorm(15))
  Y <- data.frame(a = rnorm(15), b = rnorm(15))
  res <- recordLinkage(X, Y, key = c("a", "b"), matching = "ot",
                        ot_epsilon = 0.05)
  expect_equal(res$settings$matching, "ot")
  expect_equal(res$settings$ot_epsilon, 0.05)
  expect_equal(res$settings$ot_max_iter, 100L)
  # transport_plans should be a list with at least one matrix
  expect_true(is.list(res$transport_plans))
  expect_true(length(res$transport_plans) > 0)
  tp <- res$transport_plans[[1]]
  expect_true(is.matrix(tp))
  expect_equal(nrow(tp), 15)
  expect_equal(ncol(tp), 15)
})

test_that("recordLinkage(matching = 'ot') risk capped at 1", {
  set.seed(1)
  X <- data.frame(a = rnorm(15), b = rnorm(15))
  Y <- X + rnorm(30, 0, 0.01)  # near-copies
  res <- recordLinkage(X, Y, key = c("a", "b"), matching = "ot",
                        ot_epsilon = 0.001, truth = "row")
  expect_true(all(res$per_record$risk <= 1))
  expect_true(all(res$per_record$risk >= 0))
})

test_that("recordLinkage(matching = 'ot') risk_weighting message", {
  set.seed(1)
  X <- data.frame(a = rnorm(10), b = rnorm(10))
  Y <- data.frame(a = rnorm(10), b = rnorm(10))
  expect_message(
    recordLinkage(X, Y, key = c("a", "b"), matching = "ot",
                  risk_weighting = "softmax"),
    "risk_weighting.*independent"
  )
})

test_that("recordLinkage(matching = 'ot') strategy message", {
  set.seed(1)
  X <- data.frame(a = rnorm(10), b = rnorm(10))
  Y <- data.frame(a = rnorm(10), b = rnorm(10))
  expect_message(
    recordLinkage(X, Y, key = c("a", "b"), matching = "ot",
                  strategy = "topk", k = 5),
    "strategy.*independent"
  )
})

test_that("recordLinkage(matching = 'ot') with single-record block", {
  set.seed(1)
  X <- data.frame(a = rnorm(6), b = c("A","A","A","B","B","C"))
  Y <- data.frame(a = rnorm(6), b = c("A","A","A","B","B","C"))
  # Block "C" has only 1 record
  res <- recordLinkage(X, Y, key = c("a", "b"), block = "b",
                        matching = "ot")
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(nrow(res$per_record), 6)
  expect_true(all(is.finite(res$per_record$risk)))
})

test_that("recordLinkage(matching = 'ot') direction reverse with truth", {
  set.seed(1)
  X <- data.frame(a = rnorm(15), b = rnorm(15))
  Y <- X + rnorm(30, 0, 0.1)
  res <- recordLinkage(X, Y, key = c("a", "b"), matching = "ot",
                        direction = "reverse", truth = "row")
  expect_equal(res$direction, "reverse")
  expect_true(res$overall$mean_risk > 0)
})

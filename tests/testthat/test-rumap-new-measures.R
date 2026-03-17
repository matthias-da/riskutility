test_that("rumap accepts mmd as utility measure", {
  set.seed(42)
  n <- 100
  orig <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                     key = factor(sample(letters[1:3], n, TRUE)),
                     target = factor(sample(c("A", "B"), n, TRUE)))
  synth <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                      key = factor(sample(letters[1:3], n, TRUE)),
                      target = factor(sample(c("A", "B"), n, TRUE)))
  result <- rumap(orig, list(synth1 = synth),
                  risk_measures = "ims",
                  utility_measures = "mmd",
                  key_vars = "key", target_var = "target",
                  holdout_fraction = 0.3)
  expect_s3_class(result, "rumap")
  expect_true("mmd" %in% names(result$utility))
  expect_false(is.na(result$utility$mmd[1]))
})

test_that("rumap accepts copula as utility measure", {
  set.seed(42)
  n <- 100
  orig <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                     target = factor(sample(c("A", "B"), n, TRUE)))
  synth <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                      target = factor(sample(c("A", "B"), n, TRUE)))
  result <- rumap(orig, list(synth1 = synth),
                  risk_measures = character(0),
                  utility_measures = "copula",
                  key_vars = "x1", target_var = "target")
  expect_s3_class(result, "rumap")
  expect_true("copula" %in% names(result$utility))
})

test_that("rumap accepts tail as utility measure", {
  set.seed(42)
  n <- 200
  orig <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                     target = factor(sample(c("A", "B"), n, TRUE)))
  synth <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                      target = factor(sample(c("A", "B"), n, TRUE)))
  result <- rumap(orig, list(synth1 = synth),
                  risk_measures = character(0),
                  utility_measures = "tail",
                  key_vars = "x1", target_var = "target")
  expect_s3_class(result, "rumap")
  expect_true("tail" %in% names(result$utility))
})

test_that("rumap accepts tstr as utility measure", {
  set.seed(42)
  n <- 200
  orig <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  orig$target <- orig$x1 + rnorm(n, sd = 0.5)
  synth <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  synth$target <- synth$x1 + rnorm(n, sd = 0.5)
  result <- rumap(orig, list(synth1 = synth),
                  risk_measures = character(0),
                  utility_measures = "tstr",
                  key_vars = "x1", target_var = "target")
  expect_s3_class(result, "rumap")
  expect_true("tstr" %in% names(result$utility))
})

test_that("rumap with all new utility measures together", {
  set.seed(42)
  n <- 200
  orig <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  orig$target <- orig$x1 + rnorm(n, sd = 0.5)
  synth <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  synth$target <- synth$x1 + rnorm(n, sd = 0.5)
  result <- rumap(orig, list(synth1 = synth),
                  risk_measures = character(0),
                  utility_measures = c("mmd", "tstr", "copula", "tail"),
                  key_vars = "x1", target_var = "target")
  expect_s3_class(result, "rumap")
  # All should be present
  expect_true(all(c("mmd", "tstr", "copula", "tail") %in% names(result$utility)))
})

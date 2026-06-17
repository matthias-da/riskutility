# tests/testthat/test-recordLinkage.R

# Helper data used across tests
.make_test_data <- function(n = 100, seed = 1) {
  set.seed(seed)
  x <- data.frame(
    age    = sample(18:80, n, TRUE),
    sex    = factor(sample(c("f", "m"), n, TRUE)),
    region = factor(sample(paste0("R", 1:5), n, TRUE)),
    stringsAsFactors = FALSE
  )
  # Perturbation: swap age within region
  x_anon <- x
  for (r in levels(x$region)) {
    idx <- which(x$region == r)
    if (length(idx) > 1L) {
      x_anon$age[idx] <- sample(x$age[idx])
    }
  }
  list(x = x, x_anon = x_anon)
}


# ── Class structure tests ──────────────────────────────────────────────

test_that("recordLinkage returns correct class and expected fields", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))

  expect_s3_class(res, "recordLinkageRisk")
  expect_true(is.data.frame(res$per_record))
  expect_true(all(c("risk", "cand_n", "true_in_set", "d_true", "d_min") %in%
                    names(res$per_record)))
  expect_true(is.list(res$overall))
  expect_true(is.list(res$settings))
})

test_that("recordLinkage stores n_original, n_synthetic, key_vars, method", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex"))

  expect_equal(res$n_original, 50)
  expect_equal(res$n_synthetic, 50)
  expect_equal(res$key_vars, c("age", "sex"))
  expect_equal(res$method, "deterministic")
})

test_that("per_record has correct number of rows", {
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_equal(nrow(res$per_record), 30)
})

test_that("overall has all expected summary fields", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_true(all(c("mean_risk", "max_risk", "n_high_risk", "pct_high_risk",
                     "risk_quantiles", "pct_risk_gt0", "pct_true_in_set",
                     "mean_candidate_size") %in% names(res$overall)))
})

test_that("privacy_pass is logical", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_true(is.logical(res$privacy_pass))
})


# ── Deterministic method tests ──────────────────────────────────────────

test_that("exact copy gives risk = 1 for all records (nearest)", {
  set.seed(42)
  x <- data.frame(
    a = c(1, 2, 3, 4, 5),
    b = factor(c("x", "y", "x", "y", "x"))
  )
  res <- recordLinkage(x, x, key = c("a", "b"), strategy = "nearest")
  expect_equal(res$per_record$risk, rep(1, 5))
  expect_equal(res$overall$mean_risk, 1)
})

test_that("independent random data gives low risk", {
  set.seed(123)
  n <- 200
  x <- data.frame(a = sample(1:1000, n), b = factor(sample(letters[1:10], n, TRUE)))
  x_anon <- data.frame(a = sample(1:1000, n), b = factor(sample(letters[1:10], n, TRUE)))
  res <- recordLinkage(x, x_anon, key = c("a", "b"), strategy = "nearest")
  # With 200 records and 1000*10 = 10000 possible combos, risk should be low
  expect_lt(res$overall$mean_risk, 0.3)
})

test_that("threshold strategy excludes distant records", {
  set.seed(1)
  x <- data.frame(a = c(1, 50, 100), b = factor(c("x", "y", "z")))
  x_anon <- data.frame(a = c(2, 51, 99), b = factor(c("x", "y", "z")))
  res <- recordLinkage(x, x_anon, key = c("a", "b"),
                       strategy = "threshold", threshold = 0.005)
  # Very tight threshold - some records should have 0 candidates
  expect_true(any(res$per_record$cand_n == 0))
})

test_that("topk strategy limits candidates", {
  set.seed(1)
  x <- data.frame(a = c(1, 2, 3), b = factor(c("x", "x", "x")))
  x_anon <- data.frame(a = c(1, 2, 3), b = factor(c("x", "x", "x")))
  res <- recordLinkage(x, x_anon, key = c("a", "b"),
                       strategy = "topk", k = 1)
  expect_true(all(res$per_record$cand_n >= 1))
})

test_that("topk_threshold strategy works", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       strategy = "topk_threshold", k = 3, threshold = 0.2)
  expect_s3_class(res, "recordLinkageRisk")
})

test_that("nearest_threshold strategy works", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       strategy = "nearest_threshold", threshold = 0.1)
  expect_s3_class(res, "recordLinkageRisk")
})

test_that("blocking reduces computation and gives valid results", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       block = "region")
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(nrow(res$per_record), 50)
})

test_that("return_matches stores candidate indices", {
  d <- .make_test_data(20)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       return_matches = TRUE)
  expect_true(!is.null(res$matches))
  expect_equal(length(res$matches), 20)
})


# ── Softmax weighting tests ────────────────────────────────────────────

test_that("softmax weights sum to 1", {
  w <- riskutility:::.softmax_risk(c(0.1, 0.3, 0.5))
  expect_equal(sum(w), 1, tolerance = 1e-10)
})

test_that("softmax gives highest weight to closest record", {
  w <- riskutility:::.softmax_risk(c(0.1, 0.3, 0.5))
  expect_true(w[1] > w[2])
  expect_true(w[2] > w[3])
})

test_that("softmax with equal distances gives uniform weights", {
  w <- riskutility:::.softmax_risk(c(0.5, 0.5, 0.5))
  expect_equal(w, rep(1/3, 3), tolerance = 1e-10)
})

test_that("softmax with user-supplied kappa works", {
  w <- riskutility:::.softmax_risk(c(0.1, 0.5), kappa = 10)
  expect_equal(sum(w), 1, tolerance = 1e-10)
  expect_true(w[1] > w[2])
})

test_that("softmax weighting produces valid risk values", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       risk_weighting = "softmax")
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
})


# ── Kernel weighting tests ──────────────────────────────────────────

test_that("kernel (gaussian) weights sum to 1", {
  w <- riskutility:::.kernel_risk(c(0.1, 0.3, 0.5))
  expect_equal(sum(w), 1, tolerance = 1e-10)
})

test_that("kernel gives highest weight to closest record", {
  w <- riskutility:::.kernel_risk(c(0.05, 0.3, 0.8))
  expect_true(w[1] > w[2])
  expect_true(w[2] > w[3])
})

test_that("kernel with equal distances gives uniform weights", {
  w <- riskutility:::.kernel_risk(c(0.5, 0.5, 0.5))
  expect_equal(w, rep(1/3, 3), tolerance = 1e-10)
})

test_that("kernel with user-supplied bandwidth works", {
  w <- riskutility:::.kernel_risk(c(0.1, 0.5), bandwidth = 0.2)
  expect_equal(sum(w), 1, tolerance = 1e-10)
  expect_true(w[1] > w[2])
})

test_that("kernel epanechnikov gives zero weight beyond bandwidth", {
  # bandwidth = 0.3, so d=0.5 > bandwidth -> weight 0
  w <- riskutility:::.kernel_risk(c(0.1, 0.5), bandwidth = 0.3,
                                      kernel = "epanechnikov")
  expect_equal(w[2], 0, tolerance = 1e-10)
  expect_equal(w[1], 1, tolerance = 1e-10)
})

test_that("kernel tricube gives zero weight beyond bandwidth", {
  w <- riskutility:::.kernel_risk(c(0.1, 0.5), bandwidth = 0.3,
                                      kernel = "tricube")
  expect_equal(w[2], 0, tolerance = 1e-10)
  expect_equal(w[1], 1, tolerance = 1e-10)
})

test_that("kernel single distance returns 1", {
  w <- riskutility:::.kernel_risk(0.5)
  expect_equal(w, 1)
})

test_that("kernel empty distance returns empty", {
  w <- riskutility:::.kernel_risk(numeric(0))
  expect_equal(length(w), 0)
})

test_that("kernel weighting produces valid risk values", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       risk_weighting = "kernel")
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
  expect_equal(res$settings$risk_weighting, "kernel")
  expect_equal(res$settings$kernel, "gaussian")
})

test_that("kernel differs from uniform with multi-candidate sets", {
  d <- .make_test_data(50)
  ru <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                      strategy = "threshold", threshold = 0.3,
                      risk_weighting = "uniform")
  rm <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                      strategy = "threshold", threshold = 0.3,
                      risk_weighting = "kernel")
  # With many candidates, kernel should produce different risks
  multi <- which(ru$per_record$cand_n > 2)
  if (length(multi) > 0) {
    expect_false(all(ru$per_record$risk[multi] == rm$per_record$risk[multi]))
  }
})

test_that("kernel + probabilistic method works", {
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       method = "probabilistic",
                       risk_weighting = "kernel")
  expect_s3_class(res, "recordLinkageRisk")
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
})

test_that("all three kernels produce valid output", {
  d <- .make_test_data(30)
  for (k in c("gaussian", "epanechnikov", "tricube")) {
    res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                         risk_weighting = "kernel", kernel = k)
    expect_s3_class(res, "recordLinkageRisk")
    expect_true(all(res$per_record$risk >= 0))
    expect_true(all(res$per_record$risk <= 1))
  }
})

test_that("print and summary show kernel kernel info", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       risk_weighting = "kernel",
                       kernel = "epanechnikov")
  expect_output(print(res), "epanechnikov")
  s <- summary(res)
  expect_output(print(s), "epanechnikov")
})


# ── Probabilistic method tests ─────────────────────────────────────────

test_that("probabilistic method returns fs_params", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       method = "probabilistic")
  expect_true(!is.null(res$fs_params))
  expect_true(all(c("m_probs", "u_probs", "threshold_used") %in%
                    names(res$fs_params)))
})

test_that("probabilistic: m > u for matching data", {
  # When true matches mostly agree, m should be > u for most variables
  set.seed(1)
  n <- 200
  x <- data.frame(
    a = sample(1:10, n, TRUE),
    b = factor(sample(c("x", "y", "z"), n, TRUE))
  )
  # Light perturbation: only swap a few values
  x_anon <- x
  swap_idx <- sample(n, 20)
  x_anon$a[swap_idx] <- sample(1:10, 20, TRUE)

  res <- recordLinkage(x, x_anon, key = c("a", "b"),
                       method = "probabilistic")
  # For variable b (unperturbed), m should be >= u
  expect_true(res$fs_params$m_probs["b"] >= res$fs_params$u_probs["b"])
})

test_that("probabilistic: exact copy gives high risk", {
  set.seed(42)
  x <- data.frame(
    a = c(1, 2, 3, 4, 5),
    b = factor(c("x", "y", "x", "y", "x"))
  )
  res <- recordLinkage(x, x, key = c("a", "b"),
                       method = "probabilistic")
  expect_gt(res$overall$mean_risk, 0)
})

test_that("probabilistic: user-supplied m/u works", {
  d <- .make_test_data(30)
  m <- c(age = 0.9, sex = 0.95, region = 0.95)
  u <- c(age = 0.1, sex = 0.5, region = 0.2)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       method = "probabilistic",
                       m_probs = m, u_probs = u)
  expect_equal(res$fs_params$m_probs, m)
  expect_equal(res$fs_params$u_probs, u)
})

test_that("probabilistic with softmax weighting works", {
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       method = "probabilistic",
                       risk_weighting = "softmax")
  expect_s3_class(res, "recordLinkageRisk")
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
})



# ── synth_pair method tests ─────────────────────────────────────────────

test_that("synth_pair method works correctly", {
  d <- .make_test_data(50)
  pair <- synth_pair(d$x, d$x_anon, key_vars = c("age", "sex", "region"))
  res <- recordLinkage(pair)
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(res$n_original, 50)
})

test_that("synth_pair errors without key_vars", {
  d <- .make_test_data(10)
  pair <- synth_pair(d$x, d$x_anon)
  expect_error(recordLinkage(pair), "key_vars")
})

test_that("synth_pair passes additional arguments", {
  d <- .make_test_data(50)
  pair <- synth_pair(d$x, d$x_anon, key_vars = c("age", "sex", "region"))
  res <- recordLinkage(pair, risk_weighting = "softmax")
  expect_equal(res$settings$risk_weighting, "softmax")
})


# ── S3 methods tests ────────────────────────────────────────────────────

test_that("print method runs without error", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_output(print(res), "Record Linkage Risk")
})

test_that("summary returns correct class", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  s <- summary(res)
  expect_s3_class(s, "summary.recordLinkageRisk")
})

test_that("print.summary runs without error", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  s <- summary(res)
  expect_output(print(s), "Summary: Record Linkage Risk")
})

test_that("summary for probabilistic shows fs_params", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       method = "probabilistic")
  s <- summary(res)
  expect_output(print(s), "Fellegi-Sunter")
})

test_that("plot which=1 runs without error", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_no_error(plot(res, which = 1))
})

test_that("plot which=2 runs without error", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_no_error(plot(res, which = 2))
})

test_that("plot which=3 shows message for non-probabilistic", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  # Should not error, just show placeholder text
  expect_no_error(plot(res, which = 3))
})

test_that("plot which=3 works for probabilistic method", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       method = "probabilistic")
  expect_no_error(plot(res, which = 3))
})


# ── Edge cases ──────────────────────────────────────────────────────────

test_that("single record works", {
  x <- data.frame(a = 1, b = factor("x"))
  x_anon <- data.frame(a = 1, b = factor("x"))
  res <- recordLinkage(x, x_anon, key = c("a", "b"))
  expect_equal(nrow(res$per_record), 1)
  expect_equal(res$per_record$risk, 1)
})

test_that("all NAs in anonymized data handled with na_anon='ignore'", {
  x <- data.frame(a = c(1, 2, 3), b = c("x", "y", "z"),
                  stringsAsFactors = FALSE)
  x_anon <- data.frame(a = c(NA_real_, NA_real_, NA_real_),
                        b = c(NA_character_, NA_character_, NA_character_),
                        stringsAsFactors = FALSE)
  # With ignore, all distances collapse when no valid vars remain
  expect_no_error(
    recordLinkage(x, x_anon, key = c("a", "b"), na_anon = "ignore")
  )
})

test_that("non-data.frame input is rejected", {
  expect_error(
    recordLinkage(matrix(1:4, 2, 2), data.frame(a = 1:2), key = "a"),
    "data.frame"
  )
})

test_that("missing key vars are caught", {
  x <- data.frame(a = 1:3, b = factor(c("x", "y", "z")))
  x_anon <- data.frame(a = 1:3, b = factor(c("x", "y", "z")))
  expect_error(
    recordLinkage(x, x_anon, key = c("a", "nonexistent")),
    "key"
  )
})

test_that("empty block produces zero risk", {
  x <- data.frame(
    a = c(1, 2), b = factor(c("x", "y")),
    region = factor(c("A", "B"))
  )
  x_anon <- data.frame(
    a = c(3, 4), b = factor(c("z", "w")),
    region = factor(c("C", "D"))  # different blocks
  )
  res <- recordLinkage(x, x_anon, key = c("a", "b", "region"),
                       block = "region")
  # No matching blocks, so all risks should be 0
  expect_equal(res$per_record$risk, c(0, 0))
})


# ── Property-based tests ────────────────────────────────────────────────

test_that("risk values are in [0, 1]", {
  d <- .make_test_data(100)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
})

test_that("mean risk equals mean of per-record risks", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_equal(res$overall$mean_risk, mean(res$per_record$risk),
               tolerance = 1e-10)
})

test_that("identical data gives maximum risk", {
  set.seed(1)
  n <- 20
  x <- data.frame(
    a = 1:n,
    b = factor(sample(c("x", "y"), n, TRUE))
  )
  res <- recordLinkage(x, x, key = c("a", "b"), strategy = "nearest")
  expect_equal(res$overall$mean_risk, 1)
})

test_that("deterministic method is deterministic (no randomness)", {
  d <- .make_test_data(50)
  r1 <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  r2 <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_equal(r1$per_record$risk, r2$per_record$risk)
})


# ── Backward compatibility ──────────────────────────────────────────────

test_that("default method=deterministic, risk_weighting=uniform is backward compatible", {
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_equal(res$method, "deterministic")
  expect_equal(res$settings$risk_weighting, "uniform")
})

test_that("truth='id' mode works correctly", {
  set.seed(1)
  x <- data.frame(id = paste0("id", 1:5), a = c(1, 2, 3, 4, 5),
                  b = factor(c("x", "y", "x", "y", "x")))
  x_anon <- data.frame(id = paste0("id", 1:5), a = c(1, 2, 3, 4, 5),
                        b = factor(c("x", "y", "x", "y", "x")))
  res <- recordLinkage(x, x_anon, key = c("a", "b"),
                       truth = "id", id = "id")
  expect_equal(res$overall$mean_risk, 1)
})

test_that("quantile threshold works", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       strategy = "threshold",
                       threshold = list(type = "quantile", p = 0.1))
  expect_s3_class(res, "recordLinkageRisk")
})



# ── Reverse direction tests ──────────────────────────────────────────────

test_that("reverse returns correct dimensions and fields", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       direction = "reverse")
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(nrow(res$per_record), nrow(d$x_anon))
  expect_equal(res$n_query, nrow(d$x_anon))
  expect_equal(res$direction, "reverse")
})

test_that("reverse truth='row' works with equal sizes", {
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       direction = "reverse", truth = "row")
  expect_equal(nrow(res$per_record), 30)
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
})

test_that("reverse truth='id' maps anon IDs to original IDs", {
  set.seed(1)
  x <- data.frame(id = paste0("id", 1:5), a = c(1, 2, 3, 4, 5),
                  b = factor(c("x", "y", "x", "y", "x")))
  x_anon <- data.frame(id = paste0("id", 1:5), a = c(1, 2, 3, 4, 5),
                        b = factor(c("x", "y", "x", "y", "x")))
  res <- recordLinkage(x, x_anon, key = c("a", "b"),
                       direction = "reverse", truth = "id", id = "id")
  expect_equal(res$overall$mean_risk, 1)
  expect_equal(res$direction, "reverse")
})

test_that("reverse exact copy gives risk = 1", {
  set.seed(42)
  x <- data.frame(
    a = c(1, 2, 3, 4, 5),
    b = factor(c("x", "y", "x", "y", "x"))
  )
  res <- recordLinkage(x, x, key = c("a", "b"),
                       direction = "reverse", strategy = "nearest")
  expect_equal(res$per_record$risk, rep(1, 5))
  expect_equal(res$overall$mean_risk, 1)
})

test_that("reverse deterministic differs from forward", {
  d <- .make_test_data(50)
  fwd <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  rev <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       direction = "reverse")
  # Forward has one row per original, reverse per synthetic
  expect_equal(nrow(fwd$per_record), nrow(d$x))
  expect_equal(nrow(rev$per_record), nrow(d$x_anon))
  expect_equal(fwd$direction, "forward")
  expect_equal(rev$direction, "reverse")
})

test_that("reverse probabilistic works", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       method = "probabilistic", direction = "reverse")
  expect_s3_class(res, "recordLinkageRisk")
  expect_true(!is.null(res$fs_params))
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
})

test_that("reverse blocking works", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       direction = "reverse", block = "region")
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(nrow(res$per_record), nrow(d$x_anon))
})

test_that("reverse softmax weighting gives valid risk", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       direction = "reverse", risk_weighting = "softmax")
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
})

test_that("reverse kernel weighting gives valid risk", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       direction = "reverse", risk_weighting = "kernel")
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
})

test_that("reverse print shows direction", {
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       direction = "reverse")
  expect_output(print(res), "reverse")
  expect_output(print(res), "risk per synthetic record")
})

test_that("reverse summary shows direction", {
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       direction = "reverse")
  s <- summary(res)
  expect_equal(s$direction, "reverse")
  expect_output(print(s), "reverse")
})

test_that("reverse plot works for which=1,2", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       direction = "reverse")
  expect_no_error(plot(res, which = 1))
  expect_no_error(plot(res, which = 2))
})

test_that("reverse with synth_pair works", {
  d <- .make_test_data(50)
  pair <- synth_pair(d$x, d$x_anon, key_vars = c("age", "sex", "region"))
  res <- recordLinkage(pair, direction = "reverse")
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(res$direction, "reverse")
  expect_equal(nrow(res$per_record), nrow(d$x_anon))
})

test_that("reverse with return_matches stores correct indices", {
  d <- .make_test_data(20)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       direction = "reverse", return_matches = TRUE)
  expect_true(!is.null(res$matches))
  expect_equal(length(res$matches), nrow(d$x_anon))
  # Matches should be indices into X (the search data in reverse)
  all_idx <- unlist(res$matches)
  expect_true(all(all_idx >= 1L & all_idx <= nrow(d$x)))
})

test_that("forward default backward compat: direction='forward', n_query set", {
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_equal(res$direction, "forward")
  expect_equal(res$n_query, nrow(d$x))
  expect_equal(res$settings$direction, "forward")
})


# ── Per-record enhancement tests ─────────────────────────────────────────

test_that("risk_band column exists and has correct factor levels", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_true("risk_band" %in% names(res$per_record))
  expect_s3_class(res$per_record$risk_band, "factor")
  expect_equal(levels(res$per_record$risk_band),
               c("very_low", "low", "moderate", "high", "very_high",
                 "unique_match"))
})

test_that("risk_band is consistent with risk values", {
  set.seed(42)
  x <- data.frame(a = 1:5, b = factor(c("x", "y", "x", "y", "x")))
  res <- recordLinkage(x, x, key = c("a", "b"), strategy = "nearest")
  # All risk=1 → should be "unique_match"
  expect_true(all(res$per_record$risk_band == "unique_match"))
})

test_that("d_rank column exists and is integer", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_true("d_rank" %in% names(res$per_record))
  expect_type(res$per_record$d_rank, "integer")
})

test_that("d_rank >= 1 when true_in_set is TRUE", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  in_set <- res$per_record$true_in_set
  if (any(in_set)) {
    expect_true(all(res$per_record$d_rank[in_set] >= 1L))
  }
})

test_that("d_rank is 1 for exact copy with nearest strategy", {
  x <- data.frame(a = 1:5, b = factor(c("x", "y", "z", "w", "v")))
  res <- recordLinkage(x, x, key = c("a", "b"), strategy = "nearest")
  # True match is the closest → rank 1
  expect_true(all(res$per_record$d_rank == 1L))
})

test_that("d_rank handles ties with ties.method='min'", {
  # Two identical records in x_anon -> tied distances
  x <- data.frame(a = 1, b = factor("x"))
  x_anon <- data.frame(a = c(1, 1, 2), b = factor(c("x", "x", "y")))
  res <- recordLinkage(x, x_anon, key = c("a", "b"), strategy = "nearest",
                       truth = "id", id = "a")
  # With ties.method = "min", rank should be 1 (best case)
  expect_equal(res$per_record$d_rank, 1L)
})

test_that("d_rank computed for probabilistic method", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       method = "probabilistic")
  expect_true("d_rank" %in% names(res$per_record))
  expect_true(any(!is.na(res$per_record$d_rank)))
})


# ── var_importance tests ─────────────────────────────────────────────────

test_that("var_importance is named numeric for deterministic method", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_true(!is.null(res$var_importance))
  expect_type(res$var_importance, "double")
  expect_equal(length(res$var_importance), 3)
  expect_equal(names(res$var_importance), c("age", "sex", "region"))
})

test_that("var_importance is named numeric for probabilistic method", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       method = "probabilistic")
  expect_type(res$var_importance, "double")
  expect_equal(names(res$var_importance), c("age", "sex", "region"))
})


# ── risk_gini tests ──────────────────────────────────────────────────────

test_that("risk_gini is in [0, 1]", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_true(res$overall$risk_gini >= 0)
  expect_true(res$overall$risk_gini <= 1)
})

test_that("risk_gini = 0 for uniform risk", {
  # All records have risk=1
  x <- data.frame(a = 1:5, b = factor(c("x", "y", "z", "w", "v")))
  res <- recordLinkage(x, x, key = c("a", "b"), strategy = "nearest")
  expect_equal(res$overall$risk_gini, 0, tolerance = 1e-10)
})

test_that("risk_gini in summary output", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  s <- summary(res)
  expect_true(!is.null(s$risk_gini))
  expect_output(print(s), "Gini")
})


# ── top_at_risk tests ────────────────────────────────────────────────────

test_that("top_at_risk returns correct number of rows", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  top <- top_at_risk(res, n = 5)
  expect_equal(nrow(top), 5)
  expect_true("record_id" %in% names(top))
})

test_that("top_at_risk is sorted by risk descending", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  top <- top_at_risk(res, n = 10)
  expect_true(all(diff(top$risk) <= 0))
})

test_that("top_at_risk with data appends QI columns", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  top <- top_at_risk(res, n = 5, data = d$x)
  expect_true(all(c("age", "sex", "region") %in% names(top)))
})

test_that("top_at_risk handles n > nrow", {
  d <- .make_test_data(10)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  top <- top_at_risk(res, n = 100)
  expect_equal(nrow(top), 10)
})


# ── risk_by_group tests ──────────────────────────────────────────────────

test_that("risk_by_group produces expected columns", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  rg <- risk_by_group(res, group = d$x$region)
  expect_true(all(c("mean_risk", "max_risk", "n", "n_high", "pct_high") %in%
                    names(rg)))
})

test_that("risk_by_group with column name works", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  rg <- risk_by_group(res, group = "region", data = d$x)
  expect_true("region" %in% names(rg))
  expect_equal(nrow(rg), length(unique(d$x$region)))
})

test_that("risk_by_group sorted by mean_risk descending", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  rg <- risk_by_group(res, group = d$x$region)
  expect_true(all(diff(rg$mean_risk) <= 0))
})


# ── merge_per_record tests ───────────────────────────────────────────────

test_that("merge_per_record produces correct dimensions", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  merged <- merge_per_record(res, data = d$x)
  expect_equal(nrow(merged), nrow(d$x))
  expect_true(all(c("age", "sex", "region", "risk", "d_rank", "risk_band") %in%
                    names(merged)))
})

test_that("merge_per_record errors on dimension mismatch", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_error(merge_per_record(res, data = d$x[1:5, ]))
})


# ── inspect_record tests ─────────────────────────────────────────────────

test_that("inspect_record errors without return_matches", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_error(inspect_record(res, i = 1), "return_matches")
})

test_that("inspect_record works with return_matches", {
  d <- .make_test_data(20)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       return_matches = TRUE)
  ir <- inspect_record(res, i = 1)
  expect_s3_class(ir, "inspect_record")
  expect_equal(ir$record_id, 1)
  expect_true(is.numeric(ir$risk))
})

test_that("inspect_record prints without error", {
  d <- .make_test_data(20)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       return_matches = TRUE)
  ir <- inspect_record(res, i = 1, data_orig = d$x, data_anon = d$x_anon)
  expect_output(print(ir), "Record Linkage Inspection")
  expect_true(!is.null(ir$query_record))
  expect_true(!is.null(ir$candidate_records) || ir$n_candidates == 0)
})


# ── plot which=3 for non-probabilistic methods ───────────────────────────

test_that("plot which=3 works for deterministic method", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_no_error(plot(res, which = 3))
})


# ── summary uses per_record$risk_band ────────────────────────────────────

test_that("summary uses risk_band from per_record", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  s <- summary(res)
  # Verify bands are consistent: sum should equal n_query
  expect_equal(sum(s$risk_bands), res$n_query)
})

test_that("summary includes var_importance", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  s <- summary(res)
  expect_true(!is.null(s$var_importance))
  expect_output(print(s), "Variable Importance")
})


# ── New plot types (which = 5:8) ─────────────────────────────────────────

test_that("plot which=5 (risk band barplot) works", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_no_error(plot(res, which = 5))
})

test_that("plot which=6 (rank distribution) works", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_no_error(plot(res, which = 6))
})

test_that("plot which=6 handles no valid ranks gracefully", {
  # All records in different blocks -> no true matches found
  x <- data.frame(a = c(1, 2), b = factor(c("x", "y")),
                   region = factor(c("A", "B")))
  x_anon <- data.frame(a = c(3, 4), b = factor(c("z", "w")),
                         region = factor(c("C", "D")))
  res <- recordLinkage(x, x_anon, key = c("a", "b", "region"),
                       block = "region")
  expect_no_error(plot(res, which = 6))
})

test_that("plot which=7 (risk by group) works with group vector", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_no_error(plot(res, which = 7, group = d$x$region))
})

test_that("plot which=7 works with column name and data", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_no_error(plot(res, which = 7, group = "region", data = d$x))
})

test_that("plot which=7 shows placeholder without group", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_no_error(plot(res, which = 7))
})

test_that("plot which=8 (Lorenz curve) works", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_no_error(plot(res, which = 8))
})

test_that("plot which=8 works for uniform risk (Gini=0)", {
  x <- data.frame(a = 1:5, b = factor(c("x", "y", "z", "w", "v")))
  res <- recordLinkage(x, x, key = c("a", "b"), strategy = "nearest")
  expect_no_error(plot(res, which = 8))
})

test_that("plot which=1:8 all together works", {
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  expect_no_error(plot(res, which = 1:8, group = d$x$region))
})


# ── Bijective matching (GDBRL) tests ────────────────────────────────────

test_that("bijective: exact copy yields all risk = 1", {
  skip_if_not_installed("clue")
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x, key = c("age", "sex", "region"),
                       matching = "bijective")
  expect_true(all(res$per_record$risk == 1))
})

test_that("bijective: risk is binary {0, 1}", {
  skip_if_not_installed("clue")
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       matching = "bijective")
  expect_true(all(res$per_record$risk %in% c(0, 1)))
})

test_that("bijective: assignments are unique within blocks", {
  skip_if_not_installed("clue")
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       matching = "bijective")
  assigned <- res$per_record$bijective_assigned
  nonzero <- assigned[assigned > 0]
  expect_equal(length(nonzero), length(unique(nonzero)))
})

test_that("bijective: bijective_assigned column exists", {
  skip_if_not_installed("clue")
  d <- .make_test_data(30)
  res_bij <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                           matching = "bijective")
  expect_true("bijective_assigned" %in% names(res_bij$per_record))
})

test_that("independent: bijective_assigned column absent", {
  d <- .make_test_data(30)
  res_ind <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                           matching = "independent")
  expect_false("bijective_assigned" %in% names(res_ind$per_record))
})

test_that("bijective: cand_n is 0 or 1", {
  skip_if_not_installed("clue")
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       matching = "bijective")
  expect_true(all(res$per_record$cand_n %in% c(0L, 1L)))
})

test_that("bijective: works with probabilistic method", {
  skip_if_not_installed("clue")
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       method = "probabilistic", matching = "bijective")
  expect_s3_class(res, "recordLinkageRisk")
  expect_true(all(res$per_record$risk %in% c(0, 1)))
})

test_that("bijective: works with blocking", {
  skip_if_not_installed("clue")
  d <- .make_test_data()
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       block = "region", matching = "bijective")
  expect_s3_class(res, "recordLinkageRisk")
  # Assignments should be unique within each block
  assigned <- res$per_record$bijective_assigned
  for (blk in unique(d$x$region)) {
    idx <- which(d$x$region == blk)
    blk_assigned <- assigned[idx]
    nonzero <- blk_assigned[blk_assigned > 0]
    expect_equal(length(nonzero), length(unique(nonzero)))
  }
})

test_that("bijective: works with reverse direction", {
  skip_if_not_installed("clue")
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       direction = "reverse", matching = "bijective")
  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(nrow(res$per_record), nrow(d$x_anon))
  expect_true(all(res$per_record$risk %in% c(0, 1)))
})

test_that("bijective: return_matches gives single-element lists", {
  skip_if_not_installed("clue")
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       matching = "bijective", return_matches = TRUE)
  lens <- vapply(res$matches, length, integer(1))
  expect_true(all(lens %in% c(0L, 1L)))
})

test_that("bijective: print shows 'bijective'", {
  skip_if_not_installed("clue")
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       matching = "bijective")
  out <- capture.output(print(res))
  expect_true(any(grepl("bijective", out, ignore.case = TRUE)))
})

test_that("bijective: summary includes matching field", {
  skip_if_not_installed("clue")
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       matching = "bijective")
  s <- summary(res)
  expect_equal(s$matching, "bijective")
  out <- capture.output(print(s))
  expect_true(any(grepl("bijective", out, ignore.case = TRUE)))
})

test_that("bijective: settings stores matching", {
  skip_if_not_installed("clue")
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       matching = "bijective")
  expect_equal(res$settings$matching, "bijective")
})

test_that("bijective: handles rectangular blocks via blocking (nq > ns in block)", {
  skip_if_not_installed("clue")
  # Create data where blocking produces unequal-size blocks
  set.seed(99)
  n <- 40
  x <- data.frame(
    age = sample(20:60, n, TRUE),
    sex = factor(c(rep("f", 30), rep("m", 10))),
    stringsAsFactors = FALSE
  )
  x_anon <- x
  x_anon$age <- x_anon$age + sample(-3:3, n, TRUE)
  # Blocking on sex gives blocks of size 30 and 10
  res <- recordLinkage(x, x_anon, key = c("age", "sex"),
                       block = "sex", matching = "bijective")
  expect_s3_class(res, "recordLinkageRisk")
  expect_true(all(res$per_record$risk %in% c(0, 1)))
})

## -- GDBRL correctness: hand-traced toy examples -------------------------

test_that("bijective: known permutation yields correct assignment and risk", {
  skip_if_not_installed("clue")
  # Anon rows swapped: s1=50, s2=20, s3=80 vs orig q1=20, q2=50, q3=80
  # Optimal: q1->s2(cost=0), q2->s1(cost=0), q3->s3(cost=0)
  # true_idx = [1,2,3], assigned = [2,1,3] => risk = [0,0,1]
  x <- data.frame(age = c(20, 50, 80))
  x_anon <- data.frame(age = c(50, 20, 80))
  res <- recordLinkage(x, x_anon, key = "age", matching = "bijective")
  expect_identical(res$per_record$bijective_assigned, c(2L, 1L, 3L))
  expect_identical(res$per_record$risk, c(0, 0, 1))
})

test_that("bijective: swapped anon yields all risk = 0", {
  skip_if_not_installed("clue")
  # q1(10), q2(11) vs s1(11), s2(10) -- rows swapped
  # Optimal: q1->s2(cost=0), q2->s1(cost=0). Both assigned to wrong => risk=0
  x <- data.frame(age = c(10, 11))
  x_anon <- data.frame(age = c(11, 10))
  res <- recordLinkage(x, x_anon, key = "age", matching = "bijective")
  expect_identical(res$per_record$bijective_assigned, c(2L, 1L))
  expect_true(all(res$per_record$risk == 0))
})

test_that("bijective: competing nearest neighbors resolved correctly", {
  skip_if_not_installed("clue")
  # q1(10), q2(12) both closest to s1(11), but s2(30) is far
  # Bijective: q1->s1 (cost=1), q2->s2 (cost=18) is optimal over
  #            q1->s2 (cost=20), q2->s1 (cost=1). Both get true match.
  x <- data.frame(age = c(10, 12))
  x_anon <- data.frame(age = c(11, 30))
  res <- recordLinkage(x, x_anon, key = "age", matching = "bijective")
  expect_identical(res$per_record$bijective_assigned, c(1L, 2L))
  expect_true(all(res$per_record$risk == 1))
})

test_that("bijective: matches brute-force LSAP on 4x4 problem", {
  skip_if_not_installed("clue")
  # Brute-force all 4! = 24 permutations to find optimal assignment
  brute_lsap <- function(cost) {
    n <- nrow(cost)
    all_perms <- function(v) {
      if (length(v) == 1) return(list(v))
      out <- list()
      for (i in seq_along(v)) {
        rest <- all_perms(v[-i])
        for (p in rest) out <- c(out, list(c(v[i], p)))
      }
      out
    }
    perms <- all_perms(seq_len(n))
    best_cost <- Inf; best_perm <- NULL
    for (p in perms) {
      tc <- sum(cost[cbind(seq_len(n), p)])
      if (tc < best_cost) { best_cost <- tc; best_perm <- p }
    }
    best_perm
  }
  x <- data.frame(age = c(10, 30, 50, 70))
  x_anon <- data.frame(age = c(12, 28, 53, 68))
  rng <- diff(range(c(x$age, x_anon$age)))
  cost <- outer(x$age, x_anon$age, function(a, b) abs(a - b)) / rng
  bf_assign <- brute_lsap(cost)
  res <- recordLinkage(x, x_anon, key = "age", matching = "bijective")
  expect_identical(res$per_record$bijective_assigned, bf_assign)
})

test_that("bijective: matches brute-force LSAP on 5x5 problem", {
  skip_if_not_installed("clue")
  brute_lsap <- function(cost) {
    n <- nrow(cost)
    all_perms <- function(v) {
      if (length(v) == 1) return(list(v))
      out <- list()
      for (i in seq_along(v)) {
        rest <- all_perms(v[-i])
        for (p in rest) out <- c(out, list(c(v[i], p)))
      }
      out
    }
    perms <- all_perms(seq_len(n))
    best_cost <- Inf; best_perm <- NULL
    for (p in perms) {
      tc <- sum(cost[cbind(seq_len(n), p)])
      if (tc < best_cost) { best_cost <- tc; best_perm <- p }
    }
    best_perm
  }
  x <- data.frame(age = c(10, 25, 40, 55, 70))
  x_anon <- data.frame(age = c(13, 22, 44, 52, 73))
  rng <- diff(range(c(x$age, x_anon$age)))
  cost <- outer(x$age, x_anon$age, function(a, b) abs(a - b)) / rng
  bf_assign <- brute_lsap(cost)
  res <- recordLinkage(x, x_anon, key = "age", matching = "bijective")
  expect_identical(res$per_record$bijective_assigned, bf_assign)
})

test_that("bijective: mean risk >= independent mean risk", {
  skip_if_not_installed("clue")
  set.seed(42)
  n <- 100
  x <- data.frame(
    age = sample(18:80, n, TRUE),
    sex = factor(sample(c("f", "m"), n, TRUE)),
    region = factor(sample(paste0("R", 1:5), n, TRUE))
  )
  x_anon <- x
  for (r in levels(x$region)) {
    idx <- which(x$region == r)
    x_anon$age[idx] <- sample(x$age[idx])
  }
  res_ind <- recordLinkage(x, x_anon, key = c("age", "sex", "region"),
                           matching = "independent")
  res_bij <- recordLinkage(x, x_anon, key = c("age", "sex", "region"),
                           matching = "bijective")
  expect_gte(res_bij$overall$mean_risk, res_ind$overall$mean_risk)
})

test_that("bijective: blocking constrains assignments within blocks", {
  skip_if_not_installed("clue")
  x <- data.frame(age = c(10, 20, 50, 60),
                   grp = factor(c("A", "A", "B", "B")))
  x_anon <- data.frame(age = c(12, 18, 52, 58),
                        grp = factor(c("A", "A", "B", "B")))
  res <- recordLinkage(x, x_anon, key = c("age", "grp"), block = "grp",
                       matching = "bijective")
  a <- res$per_record$bijective_assigned
  expect_true(all(a[1:2] %in% 1:2))
  expect_true(all(a[3:4] %in% 3:4))
  expect_true(all(res$per_record$risk == 1))
})

test_that("independent matching is default and unchanged", {
  d <- .make_test_data(30)
  res1 <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"))
  res2 <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                        matching = "independent")
  expect_equal(res1$per_record$risk, res2$per_record$risk)
  expect_equal(res1$settings$matching, "independent")
})


## -- Hand-computed correctness tests ------------------------------------------

test_that("deterministic: permuted anon gives correct risk and distances", {
  # X=(20,40,60), X_anon=(41,19,61). Range=42.
  # rec1: nearest=a2(1/42), true=a1(21/42). Risk=0.
  # rec2: nearest=a1(1/42), true=a2(21/42). Risk=0.
  # rec3: nearest=a3(1/42), true=a3. Risk=1.
  x <- data.frame(age = c(20, 40, 60))
  xa <- data.frame(age = c(41, 19, 61))
  res <- recordLinkage(x, xa, key = "age")
  expect_equal(res$per_record$risk, c(0, 0, 1))
  expect_equal(res$per_record$d_true[1], 21/42, tolerance = 1e-10)
  expect_equal(res$per_record$d_min[1], 1/42, tolerance = 1e-10)
  expect_equal(res$per_record$d_min[3], 1/42, tolerance = 1e-10)
})

test_that("deterministic: mixed types give correct Gower distance", {
  # d(r1,a1) = (|20-22|/10 + 0)/2 = 0.1
  # d(r1,a2) = (|20-28|/10 + 1)/2 = 0.9
  x <- data.frame(age = c(20, 30), sex = factor(c("M", "F")))
  xa <- data.frame(age = c(22, 28), sex = factor(c("M", "F")))
  res <- recordLinkage(x, xa, key = c("age", "sex"))
  expect_equal(res$per_record$risk, c(1, 1))
  expect_equal(res$per_record$d_true, c(0.1, 0.1), tolerance = 1e-10)
})

test_that("deterministic: tied distances give risk = 1/n_tied", {
  # X=(20,40), X_anon=(30,30). All d=0.5. Both tied at min.
  x <- data.frame(age = c(20, 40))
  xa <- data.frame(age = c(30, 30))
  res <- recordLinkage(x, xa, key = "age")
  expect_equal(res$per_record$risk, c(0.5, 0.5))
  expect_equal(res$per_record$cand_n, c(2L, 2L))
  expect_true(all(res$per_record$true_in_set))
})

test_that("deterministic: threshold strategy excludes/includes candidates", {
  # X=(20,80), X_anon=(21,79). Range=60. d≈0.017.
  x <- data.frame(age = c(20, 80))
  xa <- data.frame(age = c(21, 79))
  # Strict threshold: no candidates
  res1 <- recordLinkage(x, xa, key = "age", strategy = "threshold",
                         threshold = 0.01)
  expect_equal(res1$per_record$risk, c(0, 0))
  expect_equal(res1$per_record$cand_n, c(0L, 0L))
  # Generous threshold: includes true match
  res2 <- recordLinkage(x, xa, key = "age", strategy = "threshold",
                         threshold = 0.5)
  expect_equal(res2$per_record$risk, c(1, 1))
})

test_that("deterministic: topk strategy controls candidate set size", {
  # rec1(30): d=[0.2, 0.12, 0.8]. True=a1. Nearest=a2(0.12).
  x <- data.frame(age = c(30, 30, 50))
  xa <- data.frame(age = c(25, 33, 50))
  res1 <- recordLinkage(x, xa, key = "age", strategy = "topk", k = 1)
  expect_equal(res1$per_record$risk[1], 0)  # true not in top-1
  res2 <- recordLinkage(x, xa, key = "age", strategy = "topk", k = 2)
  expect_equal(res2$per_record$risk[1], 0.5)  # true in top-2
})

test_that("deterministic: nearest_threshold strategy", {
  x <- data.frame(age = c(30, 50))
  xa <- data.frame(age = c(29, 50))
  # d(30,29)=1/21≈0.048 < 0.05 → include
  res1 <- recordLinkage(x, xa, key = "age",
                         strategy = "nearest_threshold", threshold = 0.05)
  expect_equal(res1$per_record$risk[1], 1)
  # d_min > 0.01 → empty set → risk=0
  res2 <- recordLinkage(x, xa, key = "age",
                         strategy = "nearest_threshold", threshold = 0.01)
  expect_equal(res2$per_record$risk[1], 0)
})

test_that("deterministic: softmax weighting hand-computed", {
  # rec1(25): d=[0.2, 0.32] within threshold. kappa=2/0.12=16.667.
  # w_true = exp(-kappa*0.2) / (exp(-kappa*0.2) + exp(-kappa*0.32))
  x <- data.frame(age = c(25, 33, 50))
  xa <- data.frame(age = c(30, 33, 50))
  res <- recordLinkage(x, xa, key = "age", strategy = "threshold",
                        threshold = 0.5, risk_weighting = "softmax")
  kappa <- 2 / (8/25 - 5/25)
  neg_kd <- -kappa * c(5/25, 8/25)
  neg_kd <- neg_kd - max(neg_kd)
  w <- exp(neg_kd) / sum(exp(neg_kd))
  expect_equal(res$per_record$risk[1], w[1], tolerance = 1e-4)
})

test_that("deterministic: kernel (gaussian) weighting hand-computed", {
  x <- data.frame(age = c(25, 33, 50))
  xa <- data.frame(age = c(30, 33, 50))
  res <- recordLinkage(x, xa, key = "age", strategy = "threshold",
                        threshold = 0.5, risk_weighting = "kernel",
                        kernel = "gaussian", bandwidth = 0.1)
  # rec1: d=[0.2, 0.32]. u=d/0.1=[2, 3.2]. K=exp(-u^2/2).
  k1 <- exp(-2^2 / 2); k2 <- exp(-3.2^2 / 2)
  expect_equal(res$per_record$risk[1], k1 / (k1 + k2), tolerance = 1e-4)
})

test_that("deterministic: epanechnikov kernel zero beyond bandwidth", {
  # rec1(50): d to a1(10)=1, d to a2(49)=0.025. bw=0.1.
  # K(1/0.1=10) = 0 (beyond bandwidth). True=a1 → risk=0.
  x <- data.frame(age = c(50, 10))
  xa <- data.frame(age = c(10, 49))
  res <- recordLinkage(x, xa, key = "age",
                        strategy = "threshold", threshold = 1.0,
                        risk_weighting = "kernel", kernel = "epanechnikov",
                        bandwidth = 0.1)
  expect_equal(res$per_record$risk[1], 0)
})

test_that("probabilistic: hand-computed risk with user m/u (categorical)", {
  # LR_agree = log(0.9/0.3) = log(3). LR_disagree = log(0.1/0.7).
  x <- data.frame(sex = factor(c("M", "F", "M")))
  xa <- data.frame(sex = factor(c("M", "F", "F")))
  res <- recordLinkage(x, xa, key = "sex", method = "probabilistic",
                        m_probs = c(sex = 0.9), u_probs = c(sex = 0.3),
                        fs_threshold = 0)
  # rec1(M): above 0 = {a1}. True=a1. Risk=1.
  # rec2(F): above 0 = {a2,a3}. True=a2. Risk=0.5.
  # rec3(M): above 0 = {a1}. True=a3. Risk=0.
  expect_equal(res$per_record$risk, c(1, 0.5, 0))
})

test_that("probabilistic: 2-var all-agree and all-disagree patterns", {
  x <- data.frame(sex = factor(c("M","F","M")), job = factor(c("A","B","A")))
  xa <- data.frame(sex = factor(c("M","M","F")), job = factor(c("A","B","A")))
  res <- recordLinkage(x, xa, key = c("sex","job"), method = "probabilistic",
                        m_probs = c(sex = 0.9, job = 0.9),
                        u_probs = c(sex = 0.3, job = 0.3),
                        fs_threshold = 0)
  # rec1(M,A) vs a1(M,A): both agree. LR >> 0. Only a1 above 0. True=a1. Risk=1.
  # rec2(F,B): no candidates with both agree. All LR < 0. Risk=0.
  expect_equal(res$per_record$risk, c(1, 0, 0))
})

test_that("probabilistic: numeric variable with MAD-based tolerance", {
  x <- data.frame(age = c(20, 40, 60))
  xa <- data.frame(age = c(21, 39, 62))
  # MAD(c(20,40,60)) = 1.4826 * 20 ≈ 29.65
  res <- recordLinkage(x, xa, key = "age", method = "probabilistic",
                        m_probs = c(age = 0.95), u_probs = c(age = 0.2),
                        fs_threshold = 0)
  # rec1: {a1,a2} agree (within tol 29.65). Risk=1/2.
  # rec2: all 3 agree. Risk=1/3.
  # rec3: {a2,a3} agree. Risk=1/2.
  expect_equal(res$per_record$risk, c(0.5, 1/3, 0.5), tolerance = 1e-4)
})


## -- Edge case tests ----------------------------------------------------------

test_that("zero-range numeric variable: all distances = 0", {
  x <- data.frame(age = c(50, 50, 50))
  res <- recordLinkage(x, x, key = "age")
  expect_equal(res$per_record$risk, rep(1/3, 3), tolerance = 1e-10)
  expect_true(all(res$per_record$d_min == 0))
})

test_that("single candidate in block: forced match risk = 1", {
  x <- data.frame(age = c(20, 40), grp = factor(c("A", "B")))
  xa <- data.frame(age = c(25, 35), grp = factor(c("A", "B")))
  res <- recordLinkage(x, xa, key = c("age", "grp"), block = "grp")
  expect_equal(res$per_record$risk, c(1, 1))
  expect_equal(res$per_record$cand_n, c(1L, 1L))
})

test_that("all records identical: risk = 1/n", {
  x <- data.frame(age = rep(30, 5), sex = factor(rep("M", 5)))
  res <- recordLinkage(x, x, key = c("age", "sex"))
  expect_equal(res$per_record$risk, rep(0.2, 5), tolerance = 1e-10)
  expect_true(all(res$per_record$cand_n == 5L))
})

test_that("only categorical variables: correct Gower distances", {
  x <- data.frame(sex = factor(c("M","F","M")), job = factor(c("A","B","B")))
  xa <- data.frame(sex = factor(c("M","F","M")), job = factor(c("A","A","B")))
  # rec2(F,B): d to a2(F,A)=0.5, d to a3(M,B)=0.5. Tied → risk=0.5.
  res <- recordLinkage(x, xa, key = c("sex", "job"))
  expect_equal(res$per_record$risk, c(1, 0.5, 1))
})

test_that("NA handling: na_anon modes give correct distances", {
  x <- data.frame(age = c(20, 30), sex = factor(c("M", "F")))
  xa <- data.frame(age = c(NA_real_, NA_real_), sex = factor(c("M", "F")))

  # match: NA → d=0
  r1 <- recordLinkage(x, xa, key = c("age", "sex"), na_anon = "match")
  expect_equal(r1$per_record$d_true, c(0, 0))
  expect_equal(r1$per_record$risk, c(1, 1))

  # mismatch: NA → d=1 for that variable
  r2 <- recordLinkage(x, xa, key = c("age", "sex"), na_anon = "mismatch")
  expect_equal(r2$per_record$d_true, c(0.5, 0.5))
  expect_equal(r2$per_record$risk, c(1, 1))

  # ignore: NA variable excluded from distance
  r3 <- recordLinkage(x, xa, key = c("age", "sex"), na_anon = "ignore")
  expect_equal(r3$per_record$d_true, c(0, 0))
  expect_equal(r3$per_record$risk, c(1, 1))
})

test_that("ordinal factors auto-detected and scaled by range", {
  x <- data.frame(edu = ordered(c("low","high"), levels = c("low","mid","high")))
  xa <- data.frame(edu = ordered(c("low","mid"), levels = c("low","mid","high")))
  res <- recordLinkage(x, xa, key = "edu")
  # Ordinal: low=1, mid=2, high=3. Range=2.
  # rec2(high=3): d(a2=mid=2) = 1/2 = 0.5. Nearest=a2=true. Risk=1.
  expect_equal(res$per_record$d_true, c(0, 0.5), tolerance = 1e-10)
  expect_equal(res$per_record$risk, c(1, 1))
})

test_that("custom variable weights change distances correctly", {
  x <- data.frame(age = c(20, 20), sex = factor(c("M", "M")))
  xa <- data.frame(age = c(25, 20), sex = factor(c("F", "M")))
  # weights=(3,1): rec1 d(a1)=(3*1+1*1)/4=1, d(a2)=0. nearest=a2. true=a1. risk=0.
  res <- recordLinkage(x, xa, key = c("age", "sex"),
                        weights = c(age = 3, sex = 1))
  expect_equal(res$per_record$risk, c(0, 1))
})

test_that("truth='id' maps records by ID column", {
  x <- data.frame(id = c(1, 2), age = c(20, 40))
  xa <- data.frame(id = c(2, 1), age = c(38, 22))
  res <- recordLinkage(x, xa, key = "age", truth = "id", id = "id")
  # rec1(20) true=a2(22,id=1). Range=20. d_true=2/20=0.1. Nearest=a2. Risk=1.
  expect_equal(res$per_record$risk, c(1, 1))
  expect_equal(res$per_record$d_true, c(0.1, 0.1), tolerance = 1e-10)
})


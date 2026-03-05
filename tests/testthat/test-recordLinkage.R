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


# ── Midastouch weighting tests ──────────────────────────────────────────

test_that("midastouch (gaussian) weights sum to 1", {
  w <- riskutility:::.midastouch_risk(c(0.1, 0.3, 0.5))
  expect_equal(sum(w), 1, tolerance = 1e-10)
})

test_that("midastouch gives highest weight to closest record", {
  w <- riskutility:::.midastouch_risk(c(0.05, 0.3, 0.8))
  expect_true(w[1] > w[2])
  expect_true(w[2] > w[3])
})

test_that("midastouch with equal distances gives uniform weights", {
  w <- riskutility:::.midastouch_risk(c(0.5, 0.5, 0.5))
  expect_equal(w, rep(1/3, 3), tolerance = 1e-10)
})

test_that("midastouch with user-supplied bandwidth works", {
  w <- riskutility:::.midastouch_risk(c(0.1, 0.5), bandwidth = 0.2)
  expect_equal(sum(w), 1, tolerance = 1e-10)
  expect_true(w[1] > w[2])
})

test_that("midastouch epanechnikov gives zero weight beyond bandwidth", {
  # bandwidth = 0.3, so d=0.5 > bandwidth -> weight 0
  w <- riskutility:::.midastouch_risk(c(0.1, 0.5), bandwidth = 0.3,
                                      kernel = "epanechnikov")
  expect_equal(w[2], 0, tolerance = 1e-10)
  expect_equal(w[1], 1, tolerance = 1e-10)
})

test_that("midastouch tricube gives zero weight beyond bandwidth", {
  w <- riskutility:::.midastouch_risk(c(0.1, 0.5), bandwidth = 0.3,
                                      kernel = "tricube")
  expect_equal(w[2], 0, tolerance = 1e-10)
  expect_equal(w[1], 1, tolerance = 1e-10)
})

test_that("midastouch single distance returns 1", {
  w <- riskutility:::.midastouch_risk(0.5)
  expect_equal(w, 1)
})

test_that("midastouch empty distance returns empty", {
  w <- riskutility:::.midastouch_risk(numeric(0))
  expect_equal(length(w), 0)
})

test_that("midastouch weighting produces valid risk values", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       risk_weighting = "midastouch")
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
  expect_equal(res$settings$risk_weighting, "midastouch")
  expect_equal(res$settings$kernel, "gaussian")
})

test_that("midastouch differs from uniform with multi-candidate sets", {
  d <- .make_test_data(50)
  ru <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                      strategy = "threshold", threshold = 0.3,
                      risk_weighting = "uniform")
  rm <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                      strategy = "threshold", threshold = 0.3,
                      risk_weighting = "midastouch")
  # With many candidates, midastouch should produce different risks
  multi <- which(ru$per_record$cand_n > 2)
  if (length(multi) > 0) {
    expect_false(all(ru$per_record$risk[multi] == rm$per_record$risk[multi]))
  }
})

test_that("midastouch + probabilistic method works", {
  d <- .make_test_data(30)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       method = "probabilistic",
                       risk_weighting = "midastouch")
  expect_s3_class(res, "recordLinkageRisk")
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
})

test_that("all three kernels produce valid output", {
  d <- .make_test_data(30)
  for (k in c("gaussian", "epanechnikov", "tricube")) {
    res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                         risk_weighting = "midastouch", kernel = k)
    expect_s3_class(res, "recordLinkageRisk")
    expect_true(all(res$per_record$risk >= 0))
    expect_true(all(res$per_record$risk <= 1))
  }
})

test_that("print and summary show midastouch kernel info", {
  d <- .make_test_data(50)
  res <- recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                       risk_weighting = "midastouch",
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


# ── PRAM method tests ──────────────────────────────────────────────────

test_that("PRAM method returns pram_info", {
  set.seed(1)
  x <- data.frame(a = factor(c("x", "y", "x")), b = factor(c("1", "2", "1")))
  x_anon <- x  # no perturbation
  pram_m <- list(
    a = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2,
               dimnames = list(c("x", "y"), c("x", "y"))),
    b = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2,
               dimnames = list(c("1", "2"), c("1", "2")))
  )
  res <- recordLinkage(x, x_anon, key = c("a", "b"),
                       method = "pram", pram_matrix = pram_m)
  expect_true(!is.null(res$pram_info))
  expect_equal(res$pram_info$variables_used, c("a", "b"))
})

test_that("PRAM: diagonal-dominant matrix gives high risk for unperturbed data", {
  set.seed(1)
  x <- data.frame(a = factor(c("x", "y")), b = factor(c("1", "2")))
  x_anon <- x  # identical
  pram_m <- list(
    a = matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2,
               dimnames = list(c("x", "y"), c("x", "y"))),
    b = matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2,
               dimnames = list(c("1", "2"), c("1", "2")))
  )
  res <- recordLinkage(x, x_anon, key = c("a", "b"),
                       method = "pram", pram_matrix = pram_m)
  # With 99% diagonal, risk should be very high

  expect_gt(res$overall$mean_risk, 0.5)
})

test_that("PRAM errors without pram_matrix", {
  d <- .make_test_data(10)
  expect_error(
    recordLinkage(d$x, d$x_anon, key = c("age", "sex", "region"),
                  method = "pram"),
    "pram_matrix"
  )
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
  expect_equal(res$settings$method, "deterministic")
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

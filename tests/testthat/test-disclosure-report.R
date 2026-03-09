# Tests for disclosure_report (Comprehensive Disclosure Risk Report)

# Shared test data
make_report_data <- function(n = 150, seed = 42) {
  set.seed(seed)
  original <- data.frame(
    age_group = factor(sample(c("18-30", "31-45", "46-60", "60+"), n, replace = TRUE)),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE)),
    income = factor(sample(c("low", "medium", "high"), n, replace = TRUE))
  )
  synthetic <- data.frame(
    age_group = factor(sample(c("18-30", "31-45", "46-60", "60+"), n, replace = TRUE)),
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),
    region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE)),
    income = factor(sample(c("low", "medium", "high"), n, replace = TRUE))
  )
  list(original = original, synthetic = synthetic,
       key_vars = c("age_group", "gender", "region"),
       target_var = "income")
}

test_that("disclosure_report returns correct S3 class", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              key_vars = d$key_vars,
                              target_var = d$target_var,
                              compute = c("dcap", "ims"),
                              seed = 42, verbose = FALSE)

  expect_s3_class(report, "disclosure_report")
  expect_true("results" %in% names(report))
  expect_true("summary" %in% names(report))
  expect_true("overall_risk" %in% names(report))
  expect_true("n_pass" %in% names(report))
  expect_true("n_warn" %in% names(report))
  expect_true("parameters" %in% names(report))
})

test_that("disclosure_report computes attribution metrics", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              key_vars = d$key_vars,
                              target_var = d$target_var,
                              compute = "attribution",
                              seed = 42, verbose = FALSE)

  expect_true("dcap" %in% names(report$results))
  expect_true("tcap" %in% names(report$results))
  expect_true("weap" %in% names(report$results))
  expect_true("disco" %in% names(report$results))
})

test_that("disclosure_report computes distance metrics", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              compute = "distance",
                              seed = 42, verbose = FALSE)

  # At minimum IMS should work (dcr/nndr may need numeric data)
  expect_true("ims" %in% names(report$results))
})

test_that("disclosure_report computes privacy models", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              key_vars = d$key_vars,
                              target_var = d$target_var,
                              compute = "privacy",
                              seed = 42, verbose = FALSE)

  expect_true("kanonymity" %in% names(report$results))
  expect_true("ldiversity" %in% names(report$results))
  expect_true("tcloseness" %in% names(report$results))
  expect_true("individual_risk" %in% names(report$results))
})

test_that("disclosure_report computes membership metrics", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              compute = c("nnaa"),
                              seed = 42, verbose = FALSE)

  expect_true("nnaa" %in% names(report$results))
})

test_that("disclosure_report skips attribution without key_vars", {
  d <- make_report_data()
  expect_message(
    report <- disclosure_report(d$original, d$synthetic,
                                key_vars = NULL,
                                compute = "all",
                                seed = 42, verbose = TRUE),
    "Attribution metrics skipped"
  )

  # Should not contain attribution results
  expect_false("dcap" %in% names(report$results))
})

test_that("disclosure_report skips privacy without key_vars", {
  d <- make_report_data()
  expect_message(
    report <- disclosure_report(d$original, d$synthetic,
                                key_vars = NULL,
                                compute = "privacy",
                                seed = 42, verbose = TRUE),
    "Privacy model metrics skipped"
  )
})

test_that("disclosure_report accepts specific metric names", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              key_vars = d$key_vars,
                              target_var = d$target_var,
                              compute = c("dcap", "kanonymity"),
                              seed = 42, verbose = FALSE)

  expect_true("dcap" %in% names(report$results))
  expect_true("kanonymity" %in% names(report$results))
  expect_false("tcap" %in% names(report$results))
  expect_false("ims" %in% names(report$results))
})

test_that("disclosure_report overall_risk is valid", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              key_vars = d$key_vars,
                              target_var = d$target_var,
                              compute = c("dcap", "ims"),
                              seed = 42, verbose = FALSE)

  expect_true(report$overall_risk %in% c("LOW", "MEDIUM", "HIGH"))
})

test_that("disclosure_report summary table has expected columns", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              key_vars = d$key_vars,
                              target_var = d$target_var,
                              compute = c("dcap", "ims"),
                              seed = 42, verbose = FALSE)

  expect_true(all(c("metric", "value", "reference", "ratio", "status") %in%
                    names(report$summary)))
})

test_that("disclosure_report print method works", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              key_vars = d$key_vars,
                              target_var = d$target_var,
                              compute = c("dcap", "ims"),
                              seed = 42, verbose = FALSE)

  expect_output(print(report), "DISCLOSURE RISK REPORT")
  expect_output(print(report), "OVERALL")
})

test_that("disclosure_report summary method works", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              key_vars = d$key_vars,
                              target_var = d$target_var,
                              compute = c("dcap", "kanonymity"),
                              seed = 42, verbose = FALSE)

  s <- summary(report)
  expect_s3_class(s, "summary.disclosure_report")
  expect_true(length(s$attribution_results) > 0 || length(s$privacy_results) > 0)
})

test_that("disclosure_report print.summary method works", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              key_vars = d$key_vars,
                              target_var = d$target_var,
                              compute = c("dcap", "kanonymity"),
                              seed = 42, verbose = FALSE)

  s <- summary(report)
  expect_output(print(s), "Disclosure Risk Report Summary")
})

test_that("disclosure_report plot method runs without error", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              key_vars = d$key_vars,
                              target_var = d$target_var,
                              compute = c("dcap", "ims"),
                              seed = 42, verbose = FALSE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(report, which = 1))
})

test_that("disclosure_report synth_pair dispatch works", {
  d <- make_report_data()
  pair <- synth_pair(d$original, d$synthetic,
                     key_vars = d$key_vars,
                     target_var = d$target_var)

  report <- disclosure_report(pair,
                              compute = c("dcap", "ims"),
                              seed = 42, verbose = FALSE)

  expect_s3_class(report, "disclosure_report")
})

test_that("disclosure_report errors on invalid input", {
  expect_error(disclosure_report("not_df", data.frame(a = 1)), "data frame")
  expect_error(disclosure_report(data.frame(a = 1), "not_df"), "data frame")
})

test_that("disclosure_report parameters are stored", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              key_vars = d$key_vars,
                              target_var = d$target_var,
                              compute = c("dcap"),
                              seed = 42, verbose = FALSE)

  expect_equal(report$parameters$key_vars, d$key_vars)
  expect_equal(report$parameters$target_var, d$target_var)
  expect_equal(report$parameters$n_original, nrow(d$original))
  expect_equal(report$parameters$n_synthetic, nrow(d$synthetic))
})

test_that("disclosure_report counts pass/warn correctly", {
  d <- make_report_data()
  report <- disclosure_report(d$original, d$synthetic,
                              key_vars = d$key_vars,
                              target_var = d$target_var,
                              compute = c("dcap", "ims"),
                              seed = 42, verbose = FALSE)

  expect_equal(report$n_pass + report$n_warn + report$n_error + report$n_na,
               report$n_total)
})

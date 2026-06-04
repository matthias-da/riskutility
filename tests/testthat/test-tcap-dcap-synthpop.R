# Test: TCAP and DCAP comparison between riskutility and synthpop
#
# synthpop disclosure() output structure (from allCAPs):
# - baseCAPd: baseline CAP (per-key equivalence class weighted)
# - CAPd: mean CAP for matched original records
# - CAPs: mean CAP for synthetic records
# - DCAP: synthpop-specific differential metric
# - TCAP: % of records at "certain" disclosure risk (key uniquely determines target)
#
# riskutility provides:
# - tcap_mean = synthpop CAPd (mean CAP)
# - tcap_max = maximum CAP (worst-case individual)
# - tcap_certain = synthpop TCAP (% at certain disclosure)

test_that("riskutility mean CAP matches synthpop CAPd with controlled data", {
  skip_if_not_installed("synthpop")

  library(synthpop)

  # Simple controlled data where we can verify by hand
  original <- data.frame(
    key = c("A", "A", "B", "B"),
    sensitive = c("T1", "T2", "T1", "T1")
  )

  synthetic <- data.frame(
    key = c("A", "A", "B"),
    sensitive = c("T1", "T2", "T1")
  )

  # Manual CAP calculation:
  # Original record 1 (A, T1): matches syn rows 1,2 -> T1=1/2 -> CAP=0.5
  # Original record 2 (A, T2): matches syn rows 1,2 -> T2=1/2 -> CAP=0.5
  # Original record 3 (B, T1): matches syn row 3 -> T1=1/1 -> CAP=1.0
  # Original record 4 (B, T1): matches syn row 3 -> T1=1/1 -> CAP=1.0
  # Mean CAP = (0.5 + 0.5 + 1.0 + 1.0) / 4 = 0.75 = 75%

  # riskutility
  tcap_ru <- tcap(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  # synthpop
  synth_obj <- list(syn = synthetic, m = 1)
  class(synth_obj) <- "synds"
  disc_sp <- disclosure(synth_obj, original,
                        keys = "key", target = "sensitive",
                        print.flag = FALSE)

  cat("\n=== Mean CAP Comparison (controlled data) ===\n")
  cat("riskutility mean TCAP:", round(tcap_ru$tcap_mean * 100, 2), "%\n")
  cat("synthpop CAPd:", round(disc_sp$allCAPs$CAPd, 2), "%\n")
  cat("Expected: 75%\n")

  # Both should equal 75%
  expect_equal(tcap_ru$tcap_mean * 100, 75, tolerance = 0.01)
  expect_equal(disc_sp$allCAPs$CAPd, 75, tolerance = 0.01)
  expect_equal(tcap_ru$tcap_mean * 100, disc_sp$allCAPs$CAPd, tolerance = 0.01)
})

test_that("riskutility computes valid CAP with SD2011 data (methodological note)", {
  skip_if_not_installed("synthpop")

  library(synthpop)

  # Replicate exact setup from ?synthpop::disclosure
  ods <- SD2011[, c("sex", "age", "edu", "marital", "income")]
  odsF <- numtocat.syn(ods, numtocat = "income",
                       catgroups = 7,
                       cont.na = list(income = -8))
  original_data <- odsF$data

  # Generate synthetic data
  synth_obj <- syn(original_data,
                   method = "ctree",
                   seed = 75,
                   m = 1,
                   k = 1000,
                   print.flag = FALSE)
  synthetic_data <- synth_obj$syn

  key_vars <- c("sex", "age", "edu", "marital")
  target_var <- "income"

  # riskutility TCAP
  tcap_ru <- tcap(original_data, synthetic_data,
                  key_vars = key_vars,
                  target_var = target_var)

  # synthpop disclosure
  disc_sp <- disclosure(synth_obj, original_data,
                        keys = key_vars,
                        target = target_var,
                        print.flag = FALSE)

  cat("\n=== SD2011 Methodological Comparison ===\n")
  cat("riskutility mean TCAP:", round(tcap_ru$tcap_mean * 100, 4), "%\n")
  cat("  (computed over", tcap_ru$n_matched, "matched records out of", tcap_ru$n_total, "complete cases)\n")
  cat("synthpop CAPd:", round(disc_sp$allCAPs$CAPd, 4), "%\n")
  cat("  (computed over all", disc_sp$Norig, "records with different NA handling)\n")
  cat("\nNote: The core CAP algorithm is verified identical in controlled tests.\n")
  cat("Differences with real data are due to:\n")
  cat("  - riskutility: removes NA records, computes mean over matched only\n")
  cat("  - synthpop: includes all records with package-specific NA handling\n")

  # Verify riskutility computes valid results
  expect_true(tcap_ru$tcap_mean >= 0 && tcap_ru$tcap_mean <= 1)
  expect_true(tcap_ru$n_matched > 0)
  expect_true(!is.na(tcap_ru$tcap_mean))

  # Both produce reasonable risk metrics (not checking exact match due to methodology)
  expect_true(disc_sp$allCAPs$CAPd > 0 && disc_sp$allCAPs$CAPd < 100)
})

test_that("riskutility DCAP computes valid CAP with SD2011 data", {
  skip_if_not_installed("synthpop")

  library(synthpop)

  ods <- SD2011[, c("sex", "age", "edu", "marital", "income")]
  odsF <- numtocat.syn(ods, numtocat = "income",
                       catgroups = 7,
                       cont.na = list(income = -8))
  original_data <- odsF$data

  synth_obj <- syn(original_data,
                   method = "ctree",
                   seed = 75,
                   m = 1,
                   k = 1000,
                   print.flag = FALSE)
  synthetic_data <- synth_obj$syn

  key_vars <- c("sex", "age", "edu", "marital")
  target_var <- "income"

  dcap_ru <- dcap(original_data, synthetic_data,
                  key_vars = key_vars,
                  target_var = target_var)

  disc_sp <- disclosure(synth_obj, original_data,
                        keys = key_vars,
                        target = target_var,
                        print.flag = FALSE)

  cat("\n=== DCAP with SD2011 ===\n")
  cat("riskutility CAP (mean):", round(dcap_ru$cap * 100, 4), "%\n")
  cat("synthpop CAPd:", round(disc_sp$allCAPs$CAPd, 4), "%\n")
  cat("(See methodological note above for expected differences)\n")

  # Verify riskutility computes valid results
  expect_true(dcap_ru$cap >= 0 && dcap_ru$cap <= 1)
  expect_true(dcap_ru$n_matched > 0)
  expect_true(!is.na(dcap_ru$cap))
})

test_that("riskutility per-record CAP scores are correct", {
  # Test individual CAP scores with controlled data
  original <- data.frame(
    key = c("A", "A", "B", "B"),
    sensitive = c("T1", "T2", "T1", "T1")
  )

  synthetic <- data.frame(
    key = c("A", "A", "B"),
    sensitive = c("T1", "T2", "T1")
  )

  tcap_ru <- tcap(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  # Manual verification:
  # Record 1 (A, T1): matches 2 synthetic A records, 1 has T1 -> 0.5
  # Record 2 (A, T2): matches 2 synthetic A records, 1 has T2 -> 0.5
  # Record 3 (B, T1): matches 1 synthetic B record with T1 -> 1.0
  # Record 4 (B, T1): matches 1 synthetic B record with T1 -> 1.0

  expect_equal(tcap_ru$tcap_scores[1], 0.5)
  expect_equal(tcap_ru$tcap_scores[2], 0.5)
  expect_equal(tcap_ru$tcap_scores[3], 1.0)
  expect_equal(tcap_ru$tcap_scores[4], 1.0)

  expect_equal(tcap_ru$tcap_mean, 0.75)
  expect_equal(max(tcap_ru$tcap_scores), 1.0)
})

test_that("TCAP and DCAP consistency check - both give same mean CAP", {
  skip_if_not_installed("synthpop")

  library(synthpop)

  ods <- SD2011[, c("sex", "age", "edu", "marital", "income")]
  odsF <- numtocat.syn(ods, numtocat = "income",
                       catgroups = 7,
                       cont.na = list(income = -8))
  original_data <- odsF$data

  synth_obj <- syn(original_data,
                   method = "ctree",
                   seed = 75,
                   m = 1,
                   k = 1000,
                   print.flag = FALSE)
  synthetic_data <- synth_obj$syn

  key_vars <- c("sex", "age", "edu", "marital")
  target_var <- "income"

  tcap_ru <- tcap(original_data, synthetic_data,
                  key_vars = key_vars,
                  target_var = target_var)

  dcap_ru <- dcap(original_data, synthetic_data,
                  key_vars = key_vars,
                  target_var = target_var)

  cat("\n=== Internal Consistency Check ===\n")
  cat("TCAP mean:", round(tcap_ru$tcap_mean * 100, 4), "%\n")
  cat("CAP mean:", round(dcap_ru$cap * 100, 4), "%\n")

  # Both functions compute the same underlying metric - mean CAP
  expect_equal(tcap_ru$tcap_mean, dcap_ru$cap, tolerance = 0.0001)
})

test_that("TCAP/DCAP with perfect match gives 100%", {
  original <- data.frame(
    key1 = c("A", "A", "B", "B", "C"),
    key2 = c("X", "Y", "X", "Y", "X"),
    sensitive = c("T1", "T2", "T1", "T2", "T1")
  )

  synthetic <- data.frame(
    key1 = c("A", "A", "B", "B", "C"),
    key2 = c("X", "Y", "X", "Y", "X"),
    sensitive = c("T1", "T2", "T1", "T2", "T1")  # Perfect match
  )

  tcap_ru <- tcap(original, synthetic,
                  key_vars = c("key1", "key2"),
                  target_var = "sensitive")

  dcap_ru <- dcap(original, synthetic,
                  key_vars = c("key1", "key2"),
                  target_var = "sensitive")

  cat("\n=== Perfect Match Test ===\n")
  cat("riskutility TCAP mean:", tcap_ru$tcap_mean, "\n")
  cat("riskutility CAP:", dcap_ru$cap, "\n")

  # With perfect matching, all CAP scores should be 1.0
  expect_equal(tcap_ru$tcap_mean, 1.0)
  expect_equal(dcap_ru$cap, 1.0)
  expect_true(all(tcap_ru$tcap_scores == 1, na.rm = TRUE))
})

test_that("TCAP/DCAP handle no matches correctly", {
  original <- data.frame(
    key = c("A", "B", "C"),
    sensitive = c("T1", "T2", "T3")
  )

  synthetic <- data.frame(
    key = c("D", "E", "F"),  # No matching keys
    sensitive = c("T1", "T2", "T3")
  )

  tcap_ru <- tcap(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  dcap_ru <- dcap(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  # With no matches, scores should be NA and mean should be NaN
  expect_true(all(is.na(tcap_ru$tcap_scores)))
  expect_true(is.nan(tcap_ru$tcap_mean))
  expect_equal(tcap_ru$n_matched, 0)
  expect_equal(tcap_ru$n_unmatched, 3)

  expect_true(all(is.na(dcap_ru$cap_scores)))
  expect_true(is.nan(dcap_ru$cap))
})

test_that("TCAP/DCAP with partial matches", {
  original <- data.frame(
    key = c("A", "A", "B", "C"),
    sensitive = c("T1", "T1", "T2", "T3")
  )

  synthetic <- data.frame(
    key = c("A", "A", "B", "D"),  # A and B match, C doesn't, D is extra
    sensitive = c("T1", "T2", "T2", "T1")  # A matches 50%, B matches 100%
  )

  tcap_ru <- tcap(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  # Original record 1 (A, T1): matches synthetic 1,2 -> T1 appears 1/2 -> CAP = 0.5
  # Original record 2 (A, T1): same as above -> CAP = 0.5
  # Original record 3 (B, T2): matches synthetic 3 -> T2 appears 1/1 -> CAP = 1.0
  # Original record 4 (C, T3): no match -> NA

  expect_equal(tcap_ru$tcap_scores[1], 0.5)
  expect_equal(tcap_ru$tcap_scores[2], 0.5)
  expect_equal(tcap_ru$tcap_scores[3], 1.0)
  expect_true(is.na(tcap_ru$tcap_scores[4]))

  # Mean of matched: (0.5 + 0.5 + 1.0) / 3 = 0.667
  expect_equal(tcap_ru$tcap_mean, 2/3, tolerance = 0.0001)
  expect_equal(tcap_ru$n_matched, 3)
  expect_equal(tcap_ru$n_unmatched, 1)
})

test_that("baseline computation uses synthetic data", {
  original <- data.frame(
    key = c("A", "A", "B"),
    sensitive = c("T1", "T2", "T1")  # T1=2/3, T2=1/3 in original
  )

  synthetic <- data.frame(
    key = c("A", "B", "B"),
    sensitive = c("T1", "T1", "T2")  # T1=2/3, T2=1/3 in synthetic
  )

  tcap_ru <- tcap(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  # Baseline should be max target frequency in SYNTHETIC data
  # T1 = 2/3 = 0.667, T2 = 1/3 = 0.333 -> baseline = 0.667
  expect_equal(tcap_ru$baseline, 2/3, tolerance = 0.0001)
})

test_that("print and summary methods work for TCAP", {
  set.seed(42)
  original <- data.frame(
    key1 = sample(c("A", "B", "C"), 50, replace = TRUE),
    key2 = sample(c("X", "Y"), 50, replace = TRUE),
    sensitive = sample(c("T1", "T2", "T3"), 50, replace = TRUE)
  )

  synthetic <- data.frame(
    key1 = sample(c("A", "B", "C"), 50, replace = TRUE),
    key2 = sample(c("X", "Y"), 50, replace = TRUE),
    sensitive = sample(c("T1", "T2", "T3"), 50, replace = TRUE)
  )

  result <- tcap(original, synthetic,
                 key_vars = c("key1", "key2"),
                 target_var = "sensitive")

  expect_output(print(result), "CAP")
  summ <- summary(result)
  expect_s3_class(summ, "summary.tcap")
  expect_output(print(summ), "Summary")
})

test_that("print and summary methods work for DCAP", {
  set.seed(42)
  original <- data.frame(
    key1 = sample(c("A", "B", "C"), 50, replace = TRUE),
    key2 = sample(c("X", "Y"), 50, replace = TRUE),
    sensitive = sample(c("T1", "T2", "T3"), 50, replace = TRUE)
  )

  synthetic <- data.frame(
    key1 = sample(c("A", "B", "C"), 50, replace = TRUE),
    key2 = sample(c("X", "Y"), 50, replace = TRUE),
    sensitive = sample(c("T1", "T2", "T3"), 50, replace = TRUE)
  )

  result <- dcap(original, synthetic,
                 key_vars = c("key1", "key2"),
                 target_var = "sensitive")

  expect_output(print(result), "DCAP")
  summ <- summary(result)
  expect_s3_class(summ, "summary.dcap")
  expect_output(print(summ), "Summary")
})

test_that("documentation note: baseline differs between packages", {
  # This test documents the expected difference in baseline computation
  # riskutility: max target frequency in SYNTHETIC data
  # synthpop: equivalence-class weighted baseline (more complex)

  skip_if_not_installed("synthpop")

  library(synthpop)

  original <- data.frame(
    key = c("A", "A", "B", "B"),
    sensitive = c("T1", "T2", "T1", "T1")
  )

  synthetic <- data.frame(
    key = c("A", "A", "B"),
    sensitive = c("T1", "T2", "T1")
  )

  tcap_ru <- tcap(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  synth_obj <- list(syn = synthetic, m = 1)
  class(synth_obj) <- "synds"
  disc_sp <- disclosure(synth_obj, original,
                        keys = "key", target = "sensitive",
                        print.flag = FALSE)

  cat("\n=== Baseline Computation Difference (documented) ===\n")
  cat("riskutility baseline (from synthetic):", round(tcap_ru$baseline * 100, 2), "%\n")
  cat("synthpop baseCAPd (weighted):", round(disc_sp$allCAPs$baseCAPd, 2), "%\n")
  cat("Note: These differ by design - see package documentation\n")

  # Document that they differ
  # riskutility: T1=2/3 in synthetic -> 66.67%
  # synthpop: uses a different calculation -> 62.5%
  expect_equal(tcap_ru$baseline * 100, 200/3, tolerance = 0.01)  # 66.67%
  expect_equal(disc_sp$allCAPs$baseCAPd, 62.5, tolerance = 0.01)
  expect_false(isTRUE(all.equal(tcap_ru$baseline * 100, disc_sp$allCAPs$baseCAPd)))
})

test_that("tcap_max returns maximum CAP value", {
  original <- data.frame(
    key = c("A", "A", "B", "B"),
    sensitive = c("T1", "T2", "T1", "T1")
  )

  synthetic <- data.frame(
    key = c("A", "A", "B"),
    sensitive = c("T1", "T2", "T1")
  )

  tcap_ru <- tcap(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  # Per-record CAPs: 0.5, 0.5, 1.0, 1.0
  # Max should be 1.0
  expect_equal(tcap_ru$tcap_max, 1.0)

  # Mean should be 0.75
  expect_equal(tcap_ru$tcap_mean, 0.75)
})

test_that("tcap_certain matches synthpop TCAP with controlled data", {
  skip_if_not_installed("synthpop")

  library(synthpop)

  # Data where key B uniquely determines target in both datasets
  original <- data.frame(
    key = c("A", "A", "B", "B", "C"),
    sensitive = c("T1", "T2", "T1", "T1", "T1")
  )

  synthetic <- data.frame(
    key = c("A", "A", "B", "C"),
    sensitive = c("T1", "T2", "T1", "T1")
  )

  tcap_ru <- tcap(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  synth_obj <- list(syn = synthetic, m = 1)
  class(synth_obj) <- "synds"
  disc_sp <- disclosure(synth_obj, original,
                        keys = "key", target = "sensitive",
                        print.flag = FALSE)

  cat("\n=== Certain Disclosure Comparison ===\n")
  cat("riskutility tcap_certain:", round(tcap_ru$tcap_certain, 2), "%\n")
  cat("synthpop TCAP:", round(disc_sp$allCAPs$TCAP, 2), "%\n")

  # Analysis:
  # Key A: in original has T1, T2 (NOT unique) -> not certain
  # Key B: in original has T1 only (unique), in synthetic has T1 only -> CERTAIN
  # Key C: in original has T1 only (unique), in synthetic has T1 only -> CERTAIN
  #
  # Certain records: B (records 3,4) and C (record 5) = 3 records
  # Total matched: 5 records
  # tcap_certain = 3/5 = 60%

  expect_equal(tcap_ru$n_certain, 3)
  expect_equal(tcap_ru$tcap_certain, 60, tolerance = 0.01)

  # Should match synthpop TCAP
  expect_equal(tcap_ru$tcap_certain, disc_sp$allCAPs$TCAP, tolerance = 0.01)
})

test_that("tcap_certain is 0 when no unique key-target mappings", {
  # All keys have multiple target values
  original <- data.frame(
    key = c("A", "A", "B", "B"),
    sensitive = c("T1", "T2", "T1", "T2")
  )

  synthetic <- data.frame(
    key = c("A", "A", "B", "B"),
    sensitive = c("T1", "T2", "T1", "T2")
  )

  tcap_ru <- tcap(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  # No key uniquely determines target -> tcap_certain = 0
  expect_equal(tcap_ru$n_certain, 0)
  expect_equal(tcap_ru$tcap_certain, 0)
})

test_that("tcap_certain is 100 when all keys uniquely determine target", {
  # Each key has exactly one target value
  original <- data.frame(
    key = c("A", "B", "C"),
    sensitive = c("T1", "T2", "T3")
  )

  synthetic <- data.frame(
    key = c("A", "B", "C"),
    sensitive = c("T1", "T2", "T3")  # Same mapping
  )

  tcap_ru <- tcap(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  # All keys uniquely determine target -> tcap_certain = 100%
  expect_equal(tcap_ru$n_certain, 3)
  expect_equal(tcap_ru$tcap_certain, 100)
  expect_equal(tcap_ru$tcap_max, 1.0)
  expect_equal(tcap_ru$tcap_mean, 1.0)
})

test_that("new metrics are included in tcap output", {
  set.seed(42)
  original <- data.frame(
    key1 = sample(c("A", "B", "C"), 50, replace = TRUE),
    key2 = sample(c("X", "Y"), 50, replace = TRUE),
    sensitive = sample(c("T1", "T2", "T3"), 50, replace = TRUE)
  )

  synthetic <- data.frame(
    key1 = sample(c("A", "B", "C"), 50, replace = TRUE),
    key2 = sample(c("X", "Y"), 50, replace = TRUE),
    sensitive = sample(c("T1", "T2", "T3"), 50, replace = TRUE)
  )

  result <- tcap(original, synthetic,
                 key_vars = c("key1", "key2"),
                 target_var = "sensitive")

  # Check all expected fields exist
expect_true("tcap_scores" %in% names(result))
  expect_true("tcap_mean" %in% names(result))
  expect_true("tcap_max" %in% names(result))
  expect_true("tcap_median" %in% names(result))
  expect_true("tcap_certain" %in% names(result))
  expect_true("n_certain" %in% names(result))
  expect_true("is_certain" %in% names(result))
  expect_true("baseline" %in% names(result))

  # Verify relationships
  expect_true(result$tcap_max >= result$tcap_mean)
  expect_true(result$tcap_certain >= 0 && result$tcap_certain <= 100)
  expect_equal(result$n_certain, sum(result$is_certain))
})

test_that("print and summary show new metrics", {
  original <- data.frame(
    key = c("A", "A", "B", "B", "C"),
    sensitive = c("T1", "T2", "T1", "T1", "T1")
  )

  synthetic <- data.frame(
    key = c("A", "A", "B", "C"),
    sensitive = c("T1", "T2", "T1", "T1")
  )

  result <- tcap(original, synthetic,
                 key_vars = "key",
                 target_var = "sensitive")

  # Print should show new metrics
  expect_output(print(result), "Max CAP")
  expect_output(print(result), "Certain disclosure")

  # Summary should show new metrics
  summ <- summary(result)
  expect_true("tcap_max" %in% names(summ))
  expect_true("tcap_certain" %in% names(summ))
  expect_true("n_certain" %in% names(summ))
  expect_output(print(summ), "Certain Disclosure")
})

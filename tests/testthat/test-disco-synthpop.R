# Test: DiSCO implementation with both disclosure types
#
# riskutility provides two disclosure types:
#
# - "potential": Counts any synthetic record that matches an original record
#   on key+target (broader measure - any match counts)
#
# - "certain": Only counts disclosures where the key combination uniquely
#   determines the target in the original data (compatible with synthpop's DiSCO)

test_that("disco correctly counts key+target matches with potential disclosure", {
  original <- data.frame(
    key1 = c("A", "A", "B", "B", "C"),
    key2 = c("X", "Y", "X", "Y", "X"),
    sensitive = c("T1", "T2", "T1", "T2", "T1")
  )

  synthetic <- data.frame(
    key1 = c("A", "A", "B", "C", "D"),
    key2 = c("X", "Y", "X", "X", "X"),
    sensitive = c("T1", "T2", "T1", "T1", "T1")
  )

  result <- disco(original, synthetic,
                  key_vars = c("key1", "key2"),
                  target_var = "sensitive",
                  disclosure_type = "potential")

  # Manual verification:
  # Syn 1: (A,X,T1) matches Orig 1 - DISCO
  # Syn 2: (A,Y,T2) matches Orig 2 - DISCO
  # Syn 3: (B,X,T1) matches Orig 3 - DISCO
  # Syn 4: (C,X,T1) matches Orig 5 - DISCO
  # Syn 5: (D,X,T1) no key match - NOT DISCO

  expect_equal(result$n_disco, 4)
  expect_equal(result$n_disco_potential, 4)
  expect_equal(sort(result$disco_idx), c(1, 2, 3, 4))
  expect_equal(result$disclosure_type, "potential")
})

test_that("disco certain disclosure only counts unique key->target mappings", {
  # Key A maps to multiple targets (T1, T2) - NOT certain
  # Key B maps to single target (T1) - certain
  original <- data.frame(
    key = c("A", "A", "A", "B", "B"),
    sensitive = c("T1", "T1", "T2", "T1", "T1")
  )

  synthetic <- data.frame(
    key = c("A", "A", "B"),
    sensitive = c("T1", "T2", "T1")
  )

  result_potential <- disco(original, synthetic,
                            key_vars = "key",
                            target_var = "sensitive",
                            disclosure_type = "potential")

  result_certain <- disco(original, synthetic,
                          key_vars = "key",
                          target_var = "sensitive",
                          disclosure_type = "certain")

  # Potential: all 3 synthetic records match something
  expect_equal(result_potential$n_disco, 3)
  expect_equal(result_potential$n_disco_potential, 3)

  # Certain: only synthetic (B, T1) counts because B uniquely determines T1
  expect_equal(result_certain$n_disco, 1)
  expect_equal(result_certain$n_disco_certain, 1)
  expect_equal(result_certain$disco_idx, 3)  # Only the B record

  # Both results should have both counts
  expect_equal(result_potential$n_disco_certain, 1)
  expect_equal(result_certain$n_disco_potential, 3)
})

test_that("disco handles duplicates correctly", {
  original <- data.frame(
    key = c("A", "A", "A", "B"),
    sensitive = c("T1", "T1", "T2", "T1")
  )

  synthetic <- data.frame(
    key = c("A", "B"),
    sensitive = c("T1", "T1")
  )

  result <- disco(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive",
                  disclosure_type = "potential")

  # Syn 1: (A, T1) matches Orig 1 and 2 - DISCO
  # Syn 2: (B, T1) matches Orig 4 - DISCO
  expect_equal(result$n_disco, 2)

  # First synthetic record should match 2 original records
  expect_equal(length(result$matched_original_idx[[1]]), 2)

  # Second synthetic record should match 1 original record
  expect_equal(length(result$matched_original_idx[[2]]), 1)
})

test_that("disco returns zero when no matches exist", {
  original <- data.frame(
    key = c("A", "B"),
    sensitive = c("T1", "T2")
  )

  synthetic <- data.frame(
    key = c("A", "B"),
    sensitive = c("T2", "T1")  # Swapped targets
  )

  result <- disco(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  expect_equal(result$n_disco, 0)
  expect_equal(result$n_disco_potential, 0)
  expect_equal(result$n_disco_certain, 0)
  expect_equal(result$pct_disco, 0)
})

test_that("disco matches synthpop when using certain disclosure", {
  skip_on_cran()
  skip_if_not_installed("synthpop")

  library(synthpop)

  # Case where key DOES uniquely determine target
  original <- data.frame(
    key = c("A", "A", "B", "B"),
    sensitive = c("S1", "S1", "S2", "S2")  # A always S1, B always S2
  )

  synthetic <- data.frame(
    key = c("A", "B"),
    sensitive = c("S1", "S2")
  )

  # riskutility with "certain" disclosure
  disc_ru <- disco(original, synthetic,
                   key_vars = "key",
                   target_var = "sensitive",
                   disclosure_type = "certain")

  # synthpop
  disc_sp <- disclosure(synthetic, original,
                        keys = "key",
                        target = "sensitive",
                        print.flag = FALSE)

  cat("\n=== Equivalence Test (unique key->target mapping) ===\n")
  cat("riskutility (certain) n_disco:", disc_ru$n_disco, "\n")
  cat("riskutility pct_original_disclosed:", disc_ru$pct_original_disclosed, "%\n")
  cat("synthpop DiSCO:", disc_sp$attrib$DiSCO, "%\n")

  # Both should show 100% disclosure
  expect_equal(disc_ru$pct_original_disclosed, 100)
  expect_equal(disc_sp$attrib$DiSCO, 100)
})

test_that("disco certain disclosure matches synthpop when keys don't uniquely determine target", {
  skip_on_cran()
  skip_if_not_installed("synthpop")

  library(synthpop)

  # Key A -> S1 (2x), S2 (1x) - NOT unique
  # Key B -> S1 (2x) - unique
  original <- data.frame(
    key = c("A", "A", "A", "B", "B"),
    sensitive = c("S1", "S1", "S2", "S1", "S1")
  )

  synthetic <- data.frame(
    key = c("A", "A", "B"),
    sensitive = c("S1", "S2", "S1")
  )

  # riskutility with "certain" disclosure
  disc_ru <- disco(original, synthetic,
                   key_vars = "key",
                   target_var = "sensitive",
                   disclosure_type = "certain")

  # synthpop
  disc_sp <- disclosure(synthetic, original,
                        keys = "key",
                        target = "sensitive",
                        print.flag = FALSE)

  cat("\n=== synthpop Compatibility Test ===\n")
  cat("Original: Key A -> S1 (2x), S2 (1x); Key B -> S1 (2x)\n")
  cat("riskutility (certain) pct_original_disclosed:", disc_ru$pct_original_disclosed, "%\n")
  cat("synthpop DiSCO:", disc_sp$attrib$DiSCO, "%\n")

  # Both should report 40% (2/5 original records disclosed - the B records)
  expect_equal(disc_ru$pct_original_disclosed, 40)
  expect_equal(disc_sp$attrib$DiSCO, 40)
})

test_that("disco computes baseline correctly", {
  set.seed(42)

  original <- data.frame(
    key1 = sample(c("A", "B", "C"), 100, replace = TRUE),
    key2 = sample(c("X", "Y"), 100, replace = TRUE),
    sensitive = sample(c("T1", "T2", "T3"), 100, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  )

  synthetic <- data.frame(
    key1 = sample(c("A", "B", "C"), 100, replace = TRUE),
    key2 = sample(c("X", "Y"), 100, replace = TRUE),
    sensitive = sample(c("T1", "T2", "T3"), 100, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  )

  result <- disco(original, synthetic,
                  key_vars = c("key1", "key2"),
                  target_var = "sensitive")

  # Baseline should be positive and less than or equal to key matches
  expect_true(result$baseline_disco >= 0)
  expect_true(result$baseline_disco <= result$n_key_matches)

  # DiSCO ratio should be computable
  expect_true(!is.na(result$disco_ratio) || result$baseline_disco == 0)
})

test_that("disco with larger realistic dataset", {
  skip_on_cran()
  skip_if_not_installed("synthpop")

  library(synthpop)

  set.seed(123)
  n <- 500

  original <- data.frame(
    sex = sample(c("M", "F"), n, replace = TRUE),
    age_group = sample(c("18-30", "31-45", "46-60", "60+"), n, replace = TRUE),
    education = sample(c("Low", "Medium", "High"), n, replace = TRUE),
    income = sample(c("Low", "Medium", "High"), n, replace = TRUE)
  )

  synth_obj <- syn(original, seed = 456, m = 1, print.flag = FALSE)
  synthetic <- synth_obj$syn

  key_vars <- c("sex", "age_group", "education")
  target_var <- "income"

  result_potential <- disco(original, synthetic,
                            key_vars = key_vars,
                            target_var = target_var,
                            disclosure_type = "potential")

  result_certain <- disco(original, synthetic,
                          key_vars = key_vars,
                          target_var = target_var,
                          disclosure_type = "certain")

  cat("\n=== Realistic Dataset Test ===\n")
  cat("Original records:", result_potential$n_original, "\n")
  cat("Synthetic records:", result_potential$n_synthetic, "\n")
  cat("Potential DiSCO:", result_potential$n_disco_potential,
      "(", round(result_potential$pct_disco_potential, 1), "%)\n")
  cat("Certain DiSCO:", result_potential$n_disco_certain,
      "(", round(result_potential$pct_disco_certain, 1), "%)\n")
  cat("Original records disclosed:", result_certain$n_originals_disclosed,
      "(", round(result_certain$pct_original_disclosed, 1), "%)\n")

  # Basic sanity checks
  expect_true(result_potential$n_disco >= 0)
  expect_true(result_potential$n_disco <= result_potential$n_synthetic)
  expect_true(result_certain$n_disco_certain <= result_certain$n_disco_potential)
  expect_true(result_certain$pct_original_disclosed >= 0)
  expect_true(result_certain$pct_original_disclosed <= 100)
})

test_that("disco output structure is correct", {
  original <- data.frame(
    key = c("A", "B", "C"),
    sensitive = c("T1", "T2", "T1")
  )

  synthetic <- data.frame(
    key = c("A", "B"),
    sensitive = c("T1", "T2")
  )

  result <- disco(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  # Check all expected fields exist

  expect_true("disco_idx" %in% names(result))
  expect_true("disco_records" %in% names(result))
  expect_true("n_disco" %in% names(result))
  expect_true("pct_disco" %in% names(result))
  expect_true("n_disco_potential" %in% names(result))
  expect_true("n_disco_certain" %in% names(result))
  expect_true("pct_disco_potential" %in% names(result))
  expect_true("pct_disco_certain" %in% names(result))
  expect_true("n_originals_disclosed" %in% names(result))
  expect_true("pct_original_disclosed" %in% names(result))
  expect_true("disclosure_type" %in% names(result))
  expect_true("n_key_matches" %in% names(result))
  expect_true("baseline_disco" %in% names(result))

  # Check class

  expect_s3_class(result, "disco")
})

test_that("disco print and summary methods work", {
  original <- data.frame(
    key = c("A", "A", "B"),
    sensitive = c("T1", "T2", "T1")
  )

  synthetic <- data.frame(
    key = c("A", "B"),
    sensitive = c("T1", "T1")
  )

  result <- disco(original, synthetic,
                  key_vars = "key",
                  target_var = "sensitive")

  # Print should not error
  expect_output(print(result), "DiSCO")

  # Summary should not error
  summ <- summary(result)
  expect_s3_class(summ, "summary.disco")
  expect_output(print(summ), "Summary")
})

# Test: the three TCAP variants ("kinds") and their exact definitions
#
# tcap() always computes and returns all three; `kind` only selects the
# headline that print()/summary() emphasize.
#
# - tcap_certain     (kind "certain", default): % of matched original records
#   whose key uniquely determines the target in BOTH original and synthetic
#   data (and CAP = 1). riskutility's own, strictest variant.
# - tcap_matched     (kind "matched"): % of matched original records with
#   CAP = 1, i.e. disclosive in the synthetic data and correctly attributed.
#   Reproduces synthpop's TCAP as defined up to synthpop 1.9-2
#   (sum(tab_DiSCO) / sum(tab_kd[tab_ks > 0])).
# - tcap_conditional (kind "conditional"): same numerator, divided by the
#   number of original records matched to a *disclosive* synthetic key class.
#   The definition of Little, Allmendinger & Elliot (2025, JOS 41(1)),
#   adopted by synthpop from version 1.9-3. Their reference implementation
#   returns 0 when no record is matched to a disclosive class.
#
# On any data: tcap_certain <= tcap_matched <= tcap_conditional.

test_that("all three kinds are hand-verified on controlled data", {
  # Syn classes: A -> {T1, T2} (not disclosive); B -> {T1}, C -> {T1} (disclosive)
  original <- data.frame(
    key = c("A", "A", "B", "B", "C"),
    sensitive = c("T1", "T2", "T1", "T1", "T1")
  )
  synthetic <- data.frame(
    key = c("A", "A", "B", "C"),
    sensitive = c("T1", "T2", "T1", "T1")
  )

  r <- tcap(original, synthetic, key_vars = "key", target_var = "sensitive")

  # CAP scores: 0.5, 0.5, 1, 1, 1 -> matched 5, CAP=1 count 3
  # certain: keys B, C unique in both -> records 3,4,5
  # disclosive-matched: records in classes B, C -> 3
  expect_equal(r$n_matched, 5)
  expect_equal(r$n_certain, 3)
  expect_equal(r$n_disclosive_correct, 3)
  expect_equal(r$n_disclosive, 3)

  expect_equal(r$tcap_certain, 60)          # 3 / 5
  expect_equal(r$tcap_matched, 60)          # 3 / 5  (synthpop <= 1.9-2: 60)
  expect_equal(r$tcap_conditional, 100)     # 3 / 3  (synthpop >= 1.9-3: 100)
})

test_that("the three kinds separate on data where the definitions differ", {
  # Syn classes: A -> {T1} disclosive, B -> {T2} disclosive, D -> {T3,T4} not.
  # Key A is NOT unique in the original ({T1, T1, T2}) -> certain < matched.
  # Class D is matched but not disclosive -> matched < conditional.
  original <- data.frame(
    key = c("A", "A", "A", "B", "D", "D"),
    sensitive = c("T1", "T1", "T2", "T2", "T3", "T4")
  )
  synthetic <- data.frame(
    key = c("A", "B", "D", "D"),
    sensitive = c("T1", "T2", "T3", "T4")
  )

  r <- tcap(original, synthetic, key_vars = "key", target_var = "sensitive")

  # CAP scores: 1, 1, 0, 1, 0.5, 0.5
  expect_equal(r$tcap_scores, c(1, 1, 0, 1, 0.5, 0.5))
  expect_equal(r$n_matched, 6)
  expect_equal(r$n_certain, 1)              # only record 4 (key B unique in both)
  expect_equal(r$n_disclosive_correct, 3)   # records 1, 2, 4 (CAP = 1)
  expect_equal(r$n_disclosive, 4)           # records 1, 2, 3 (class A), 4 (class B)

  expect_equal(r$tcap_certain, 100 / 6)     # 16.67
  expect_equal(r$tcap_matched, 50)          # 3 / 6
  expect_equal(r$tcap_conditional, 75)      # 3 / 4

  # per-record disclosive flag is consistent with the count
  expect_length(r$is_disclosive, nrow(original))
  expect_equal(sum(r$is_disclosive), r$n_disclosive)
  expect_equal(which(r$is_disclosive), 1:4)
})

test_that("kinds always order as certain <= matched <= conditional", {
  for (seed in c(1, 7, 42)) {
    set.seed(seed)
    X <- data.frame(
      k1 = sample(LETTERS[1:6], 80, replace = TRUE),
      k2 = sample(c("X", "Y"), 80, replace = TRUE),
      t = sample(paste0("T", 1:3), 80, replace = TRUE)
    )
    Y <- data.frame(
      k1 = sample(LETTERS[1:6], 60, replace = TRUE),
      k2 = sample(c("X", "Y"), 60, replace = TRUE),
      t = sample(paste0("T", 1:3), 60, replace = TRUE)
    )
    r <- tcap(X, Y, key_vars = c("k1", "k2"), target_var = "t")
    expect_true(r$n_matched > 0)
    expect_lte(r$tcap_certain, r$tcap_matched + 1e-9)
    if (r$n_disclosive > 0) {
      expect_lte(r$tcap_matched, r$tcap_conditional + 1e-9)
    }
  }
})

test_that("no key matches: matched is NA, conditional is 0", {
  original <- data.frame(key = c("A", "B", "C"), sensitive = c("T1", "T2", "T3"))
  synthetic <- data.frame(key = c("D", "E", "F"), sensitive = c("T1", "T2", "T3"))

  r <- tcap(original, synthetic, key_vars = "key", target_var = "sensitive")

  expect_true(is.na(r$tcap_matched))
  expect_equal(r$tcap_conditional, 0)   # reference-implementation convention
  expect_equal(r$n_disclosive, 0)
  expect_equal(r$n_disclosive_correct, 0)

  # print() must not error on a zero-match object
  expect_output(print(r), "undefined")
})

test_that("matches but no disclosive synthetic class: all kinds are 0", {
  original <- data.frame(
    key = c("A", "A", "B", "B"),
    sensitive = c("T1", "T2", "T1", "T2")
  )
  synthetic <- original  # every synthetic class has two target values

  r <- tcap(original, synthetic, key_vars = "key", target_var = "sensitive")

  expect_equal(r$n_matched, 4)
  expect_equal(r$n_disclosive, 0)
  expect_equal(r$tcap_certain, 0)
  expect_equal(r$tcap_matched, 0)
  expect_equal(r$tcap_conditional, 0)
})

test_that("perfect unique mapping gives 100 for all kinds", {
  original <- data.frame(key = c("A", "B", "C"), sensitive = c("T1", "T2", "T3"))
  synthetic <- original

  r <- tcap(original, synthetic, key_vars = "key", target_var = "sensitive")

  expect_equal(r$tcap_certain, 100)
  expect_equal(r$tcap_matched, 100)
  expect_equal(r$tcap_conditional, 100)
  expect_equal(r$n_disclosive, 3)
})

test_that("kind selects the headline but never changes the values", {
  original <- data.frame(
    key = c("A", "A", "A", "B", "D", "D"),
    sensitive = c("T1", "T1", "T2", "T2", "T3", "T4")
  )
  synthetic <- data.frame(
    key = c("A", "B", "D", "D"),
    sensitive = c("T1", "T2", "T3", "T4")
  )

  r_def <- tcap(original, synthetic, key_vars = "key", target_var = "sensitive")
  r_cond <- tcap(original, synthetic, key_vars = "key", target_var = "sensitive",
                 kind = "conditional")

  expect_equal(r_def$kind, "certain")
  expect_equal(r_cond$kind, "conditional")

  expect_equal(r_cond$tcap_certain, r_def$tcap_certain)
  expect_equal(r_cond$tcap_matched, r_def$tcap_matched)
  expect_equal(r_cond$tcap_conditional, r_def$tcap_conditional)

  expect_error(
    tcap(original, synthetic, key_vars = "key", target_var = "sensitive",
         kind = "nonsense"),
    "arg"
  )
})

test_that("kind passes through the synth_pair method", {
  original <- data.frame(key = c("A", "A", "B"), sensitive = c("T1", "T2", "T1"))
  synthetic <- data.frame(key = c("A", "B"), sensitive = c("T1", "T1"))

  pair <- synth_pair(original, synthetic, key_vars = "key", target_var = "sensitive")
  r_sp <- tcap(pair, kind = "matched")
  r_df <- tcap(original, synthetic, key_vars = "key", target_var = "sensitive")

  expect_equal(r_sp$kind, "matched")
  expect_equal(r_sp$tcap_matched, r_df$tcap_matched)
  expect_equal(r_sp$tcap_conditional, r_df$tcap_conditional)
})

test_that("print shows all three variants and marks the selected kind", {
  original <- data.frame(
    key = c("A", "A", "B", "B", "C"),
    sensitive = c("T1", "T2", "T1", "T1", "T1")
  )
  synthetic <- data.frame(
    key = c("A", "A", "B", "C"),
    sensitive = c("T1", "T2", "T1", "T1")
  )

  r_def <- tcap(original, synthetic, key_vars = "key", target_var = "sensitive")
  out_def <- capture.output(print(r_def))
  expect_true(any(grepl("Certain disclosure.*selected", out_def)))
  expect_true(any(grepl("TCAP matched", out_def)))
  expect_true(any(grepl("TCAP conditional", out_def)))

  r_cond <- tcap(original, synthetic, key_vars = "key", target_var = "sensitive",
                 kind = "conditional")
  out_cond <- capture.output(print(r_cond))
  expect_true(any(grepl("TCAP conditional.*selected", out_cond)))
  expect_false(any(grepl("Certain disclosure.*selected", out_cond)))
})

test_that("summary carries all variants and prints them", {
  original <- data.frame(
    key = c("A", "A", "A", "B", "D", "D"),
    sensitive = c("T1", "T1", "T2", "T2", "T3", "T4")
  )
  synthetic <- data.frame(
    key = c("A", "B", "D", "D"),
    sensitive = c("T3", "T2", "T3", "T4")
  )

  r <- tcap(original, synthetic, key_vars = "key", target_var = "sensitive")
  s <- summary(r)

  expect_equal(s$tcap_matched, r$tcap_matched)
  expect_equal(s$tcap_conditional, r$tcap_conditional)
  expect_equal(s$n_disclosive, r$n_disclosive)
  expect_equal(s$n_disclosive_correct, r$n_disclosive_correct)
  expect_equal(s$kind, r$kind)

  expect_output(print(s), "conditional")
  expect_output(print(s), "Little")
})

test_that("gower method: certain equals matched, conditional uses homogeneous match sets", {
  # Factor key with tiny threshold makes Gower matching identical to exact
  # matching, so the exact-method hand trace carries over. For gower,
  # is_certain is the simplified CAP = 1 criterion, hence certain == matched.
  original <- data.frame(
    key = factor(c("A", "A", "A", "B", "D", "D"), levels = c("A", "B", "D")),
    sensitive = c("T1", "T1", "T2", "T2", "T3", "T4")
  )
  synthetic <- data.frame(
    key = factor(c("A", "B", "D", "D"), levels = c("A", "B", "D")),
    sensitive = c("T1", "T2", "T3", "T4")
  )

  r <- tcap(original, synthetic, key_vars = "key", target_var = "sensitive",
            method = "gower", gower_threshold = 0.01)

  expect_equal(r$n_matched, 6)
  expect_equal(r$n_disclosive_correct, 3)
  expect_equal(r$tcap_matched, 50)
  expect_equal(r$tcap_certain, r$tcap_matched)
  # records 5, 6 match syn class D with mixed targets {T3, T4} -> not homogeneous
  expect_equal(r$n_disclosive, 4)
  expect_equal(r$tcap_conditional, 75)
})

test_that("continuous target: kinds are computed on the binned target", {
  set.seed(99)
  original <- data.frame(
    key = sample(c("A", "B", "C"), 60, replace = TRUE),
    income = rlnorm(60, 10, 0.5)
  )
  synthetic <- data.frame(
    key = sample(c("A", "B", "C"), 60, replace = TRUE),
    income = rlnorm(60, 10, 0.5)
  )

  r <- tcap(original, synthetic, key_vars = "key", target_var = "income",
            cont_bins = 4)

  expect_false(is.null(r$tcap_matched))
  expect_false(is.null(r$tcap_conditional))
  expect_true(r$tcap_matched >= 0 && r$tcap_matched <= 100)
  expect_true(r$tcap_conditional >= 0 && r$tcap_conditional <= 100)
})

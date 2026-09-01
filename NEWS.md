# riskutility 0.2.0

Bug fixes and one API addition, all surfaced while preparing the JSS software
paper. Three of these change reported values, so results produced with 0.1.0
are not directly comparable.

## New

* `tcap()` gains `cont_bins` (default 10), matching `dcap()`. A continuous
  target is now discretised into quantile bins before matching. Previously an
  exact match on a continuous value essentially never occurred and every TCAP
  score collapsed to zero.
* `tcap()` now reports all three TCAP variants in use in the literature and
  gains a `kind` argument (`"certain"`, the default; `"matched"`;
  `"conditional"`) that selects which one `print()`/`summary()` highlight --
  reported values are unchanged. `tcap_certain` keeps riskutility's strict
  definition (the key uniquely determines the target in both the original and
  the synthetic data). The new `tcap_matched` reproduces synthpop's TCAP as
  defined up to synthpop 1.9-2; the new `tcap_conditional` implements the
  definition of Little, Allmendinger & Elliot (2025, Journal of Official
  Statistics 41(1), 255-308), which synthpop adopts from version 1.9-3. The
  previous documentation claim that `tcap_certain` equals synthpop's TCAP held
  only on data where every disclosive synthetic key class is also unambiguous
  in the original data; the documentation now states the exact mapping.

## Bug fixes (change reported values)

* `propscore()`: the synthetic-side kernel density was indexed with
  `length(ps)` -- a scalar -- rather than `nrow(p)`, so the index ran backwards
  and every density diagnostic reported by `summary()` (`kl`, `density_ratio`,
  `mean_ratio`, `sd_ratio` and their Bayes-space counterparts) was computed on
  the wrong records.
* `tstr()`: R-squared is now clamped at 0, and the TSTR/TRTR ratio is returned
  as `NA` with a warning when the train-on-real model has no predictive power.
  Previously two negative R-squared values could divide to a large positive
  ratio and be reported as excellent utility.
* `tail_fidelity()`: the Jensen-Shannon divergence is now computed in density
  space. It previously used Bayes-space (clr) ordinates, which are
  sign-indefinite, and could return values far outside `[0, log 2]`.

## Other

* `dcr()` now warns when the training and holdout sets differ in size. The
  `dcr_share` reference value of 0.5 and the 0.55 rule used by `privacy_pass`
  assume equal sizes; with an unequal split the expected share under no
  memorisation is `n_train / (n_train + n_holdout)`.
* Tests that compare riskutility against synthpop output are now version-aware
  (synthpop 1.9-3 redefines its TCAP measure) and are skipped on CRAN, as
  requested by the synthpop maintainers and the CRAN team. They continue to
  run in local `devtools::test()` runs.
* `print.tcap()` no longer errors on a result with zero matched records.
* `compare_feature_importance(importance_type = "permutation")` no longer
  requires the vip package (archived from CRAN on 2026-07-08 and therefore
  flagged by the CRAN incoming checks): permutation importance is computed
  internally, and the documented caret-style metric names (`"Accuracy"`,
  `"Kappa"`, `"ROC AUC"`, `"RMSE"`, `"MAE"`, `"MAPE"`, `"R-squared"`) as well
  as custom metric functions are handled directly, with error metrics
  sign-flipped so that larger importance always means a more influential
  feature. vip was removed from `Suggests`.

# riskutility 0.1.0

Initial release: a comprehensive framework for measuring disclosure risk and
data utility of anonymized and synthetic data. All measures share a consistent
S3 API (`print()`, `summary()`, `plot()`) and feed a multivariate Risk-Utility
(R-U) map.

## Disclosure risk

* **Attribution-based (CAP family):** `dcap()` (reports both the raw mean CAP and
  the differential CAP = mean CAP minus baseline), `tcap()`, `weap()`, `disco()`.
* **ML-based:** `rapid()` (Risk of Attribute Prediction-Induced Disclosure;
  random-forest default, also `lm`/`cart`/`gbm`/`logit`) with `confint()`,
  permutation test, threshold selection, synthesizer cross-validation, and six
  plot types.
* **Distance-based (holdout):** `dcr()`, `nndr()`, `ims()`, `repu()`, including the
  DCR-Delusion caveat and null-distribution diagnostics.
* **Membership inference:** `domias()`, `nnaa()`, `mia_classifier()`.
* **Classical SDC privacy models:** `kanonymity()`, `ldiversity()`
  (distinct/entropy/recursive), `tcloseness()` (EMD), `suda()`,
  `individual_risk()`, `population_uniqueness()` (Pitman/Zayatz/SNB),
  `epsilon_identifiability()`, `delta_presence()`, `hitting_rate()`,
  `singling_out()`, `linkability()`, `attacker_risk()`
  (prosecutor/journalist/marketer), `drisk()`.
* **Record linkage:** `recordLinkage()` with deterministic, probabilistic
  (Fellegi-Sunter), PRAM, predictive, random-forest, RBRL, robust-Mahalanobis,
  and embedding (autoencoder) methods; independent, bijective (Hungarian / GDBRL),
  and optimal-transport (Sinkhorn) matching; blocking and per-record accessors.
  All eight methods share a single re-identification-risk definition — the
  probability of identifying the *true* match within the attacker's candidate
  set. For the random-forest and embedding methods, the nearest-neighbour
  similarity (their former `risk` value) is now retained in an `nn_similarity`
  diagnostic column. `na_anon` (`ignore`/`match`/`mismatch`) is honored
  consistently across all methods (PRAM no longer reports an artificial zero
  risk for records with a missing key). New options: `compute_baseline = TRUE`
  reports the no-perturbation reference risk (with `risk_reduction`), and
  `expected_risk = TRUE` reports a perturbation-aware expected PRAM risk over
  the transition distribution. User-supplied `m_probs`/`u_probs` are validated
  and clamped to the open interval (0,1).
* **Reporting:** `disclosure_report()` produces a comprehensive multi-metric report.

## Data utility

* Propensity-score utility: `propscore()`, `pMSE()`, `specks()`.
* Global / interval: `gower()`, `mqs()`, `ci_overlap()`, `ci_proximity()`.
* Distributional and structural: `compare_wasserstein()`, `compare_ks_test()`,
  `compare_chisq_gof()`, `compare_pca()`, `compare_embedding()`,
  `compare_correlation_matrices()`, `hellinger()`, `energy_distance()`, `mmd()`,
  `copula_fidelity()`, `tail_fidelity()`, `contingency_fidelity()`.
* Downstream / model-based: `tstr()` (train on synthetic, test on real),
  `compare_feature_importance()`, `compare_model_performance()`,
  `regression_fidelity()`, `subgroup_utility()`.
* Information-theoretic: `KLDiv()`, `JSDiv()`, `CrossEntropy()`, entropy and
  mutual-information helpers, `privacy_score()`.

## Multivariate Risk-Utility map

* `rumap()`: normalized multivariate R-U evaluation with Pareto-frontier
  identification, internal-consistency metrics, and seven visualizations
  (scatter, heatmap, dot plot, parallel coordinates, radial, PCA biplot,
  blockwise PCA).

## Integration

* `synth_pair()` container plus `from_synthpop()` and `from_simPop()` converters;
  most measures dispatch on `synth_pair` objects as well as plain data frames.

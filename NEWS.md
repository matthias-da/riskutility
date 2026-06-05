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

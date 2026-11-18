# riskutility 0.3.0

RAPID (`rapid()` and its companion functions) is rebuilt on the measure of the RAPID paper
(Templ, Thees and Müller, 2026, under review). It measures a gain over a marginal baseline for
both types of target and reports the chance level in every result. The new vignette
`vignette("rapid", package = "riskutility")` explains the measure and the result object.
Several defaults change reported values. Each item below gives the call that reproduces the
0.2.0 rule, as in the section "Changes in 0.3.0" of `?rapid`.

## Value-changing

* For a continuous target the default rule is the band rule (`num_eval_method = "band"`). A
  record is at risk when the attacker's predictive distribution puts clearly more probability
  into a tolerance band around the true value than the marginal distribution of the original
  values does (gain `r = (g - b) / (1 - b) > num_tau`, default `num_tau = 0.3`). The 0.2.0 rule
  is `rapid(X, Y, kv, tv, num_eval_method = "point", num_epsilon = 0.1)`.
* For a categorical target the default rule is the gain over the marginal baseline
  (`cat_eval_method = "RCS_marginal"`, `cat_tau = 0.3`), as in the paper. The 0.2.0 rule is
  `rapid(X, Y, kv, tv, cat_eval_method = "RCS_conditional")` (its default `cat_tau` stays 1).
* A threshold on a gain (`num_tau` under the band rule, `cat_tau` under `RCS_marginal`, also
  through `update()`) must be at least 0 and below 1. Another value is an error of class
  `rapid_error_input`. A gain is at most 1, so a threshold of 1 or more flags no record, and a
  negative threshold flags records on which the attacker does worse than the marginal
  distribution. A 0.2.0 call with `cat_tau = 1` is therefore an error that names its rule, `rapid(X,
  Y, kv, tv, cat_tau = 1, cat_eval_method = "RCS_conditional")`.
* The tolerance `num_epsilon` of the relative error metrics is a proportion (default 0.1, that
  is 10 %), as in the paper. `num_epsilon = 10` is now `num_epsilon = 0.1`. A relative tolerance
  of 1 or more is still read as percent, with a warning.
* `rapid_threshold_select()` under the point rule takes its grid `epsilon_range` as proportions
  (default `seq(0.01, 0.5, by = 0.02)`) and recommends the smallest significant tolerance, the
  strictest one. 0.2.0 recommended the largest, the least strict. The 0.2.0 choice is
  `max(res$results$threshold[res$results$significant])`.
* The linear model chooses its working scale (`num_scale = "auto"`). It takes the log scale when
  the error metric is relative, all synthetic target values are positive and the log model has
  the higher likelihood on the original scale. Imputed synthetic values take no part in this
  choice. The 0.2.0 model is `model_type = "lm", num_scale = "identity"`.
* The tree of `model_type = "cart"` is grown with `cp = 0`, `minbucket = 5`, `xval = 10` and
  pruned at the complexity with the smallest cross-validated error. The 0.2.0 tree, with the
  rpart defaults, is
  `model_type = "cart", prune = FALSE, cp = 0.01, minsplit = 20, minbucket = 7`. `minsplit` must
  be given, because with `minbucket = 7` alone rpart sets `minsplit = 21`.
* `rapid_test()`, `rapid_threshold_select()`, `rapid_qi_contribution()` and
  `rapid_synthesizer_cv()` fit through `rapid()` and use its defaults. Their numeric results
  change with the band rule as above.
* `rumap()` reports for RAPID the excess of a random forest attacker over its chance level,
  `max(0, rapid - chance)`. 0.2.0 reported the raw RAPID of a linear model under the point rule,
  which is now
  `rapid(X, Y, kv, tv, model_type = "lm", num_eval_method = "point", num_scale = "identity",
  num_epsilon = 0.1)$rapid`.
* `disclosure_report()` passes RAPID when the p-value of `rapid_test()` (199 permutations) is
  above the new argument `rapid_alpha` (default 0.05). Its row reports the chance level and the
  excess. 0.2.0 passed RAPID when it was below 0.15 under the 0.2.0 rule. With the default
  attacker `rapid_model = "rf"` that verdict is now
  `rapid(X, Y, kv, tv, model_type = "rf", num_eval_method = "point",
  num_epsilon = 0.1)$rapid < 0.15`
  for a continuous target and
  `rapid(X, Y, kv, tv, model_type = "rf", cat_eval_method = "RCS_conditional")$rapid < 0.15`
  for a categorical one. For `rapid_model = "lm"` or `"cart"` add the settings of the
  attacker items above (`num_scale = "identity"`, or `prune = FALSE, cp = 0.01, minsplit = 20,
  minbucket = 7`).
* `summary()` of a `rapid` object no longer labels the result LOW, MEDIUM or HIGH (cut-offs 0.05
  and 0.15). Read the excess over the chance level and `rapid_test()` instead.
* `rapid_synthesizer_cv()` has no default seed and is silent by default. The 0.2.0 call is
  `rapid_synthesizer_cv(..., seed = 2025, verbose = TRUE)`. No RAPID function sets a seed unless
  one is given.
* The predicted class of a record (column `predicted_class`) takes the first class with the
  highest probability. 0.2.0 broke ties at random.

## New

* Every `rapid()` result reports the chance level `chance`, the share of records flagged when
  the attacker's predictions are scored against permuted true values (`n_perm`, default 10), the
  shares of the single permutations (`chance_rates`) and the excess `excess = rapid - chance`.
  The print says when RAPID is not distinguishable from the chance level. The functions that
  build on `rapid()` report the three values as well.
* `threshold_rule = "permutation"` replaces the threshold by the `perm_level` quantile (default
  0.95) of all permuted gains.
* `num_scale` (`"auto"`, `"identity"`, `"log"`) sets the working scale of the linear model.
* New lookup attacker `model_type = "lookup"`, one fully grown ranger tree that reads off the
  synthetic records with the same key values. It measures whether the release reproduces the
  records' values and is read by its excess. The print adds a line that says so.
* The print of a `rapid` object starts with the attack scenario and an attacker line with the
  settings that drive the result (for example `rf: 500 trees, min.node.size 5` or
  `cart: 12 leaves, pruned by cross-validation`), stored in the field `attacker`.
* Continuous targets get a predictive distribution per attacker. It is normal or lognormal for
  lm, the leaf values for cart and lookup, a quantile regression forest for rf, a quantile model
  for gbm. `num_dist` chooses between the exact probability (`"exact"`), `n_quantiles` quantiles
  (`"grid"`, default 100) and the exact form where the attacker has one (`"auto"`, default).
* With `num_dist = "exact"` the forest uses the weights of quantile regression forests
  (Meinshausen, 2006) instead of ranger's random draws, and lm and cart use their own
  distribution.
* The record table holds the band (`band_lo`, `band_hi`), the confidence `g`, the baseline `b`
  and the gain `r` of every record.
* Arguments in `...` reach the attacker model (ranger, lm, rpart, glm, xgboost), and `prune`
  sets whether the tree is pruned.
* `rapid_suite()` fits several attackers on the same data (default `rf`, `cart`, `gbm`,
  `lookup`) or combines fitted `rapid` objects, and reports the worst case, the largest RAPID
  and the largest excess with the model that has it. It has `print()`, `summary()`, `plot()`,
  `as.data.frame()` and `[[`.
* `rapid_learner()` builds a user-defined attacker from a fitting function and prediction
  functions. It is accepted wherever `model_type` is.
* `rapid(..., holdout = )` (or the holdout of a `synth_pair`) scores records that were not used
  for the release with the same attacker and rule. The result holds their RAPID, chance level
  and excess (`holdout`) and `release_specific`, RAPID minus holdout RAPID.
* `update()` changes arguments of a `rapid` object. It rescores the stored predictions without a
  refit when only the rule changes, and refits otherwise.
* `predict()` applies the attacker of a `rapid` object fitted with `store_model = TRUE` to new
  records, and scores them when they hold the target.
* `as.data.frame()` for `rapid` objects (the record table, `at_risk_only = TRUE` for the flagged
  records) and for the objects of `rapid_test()`, `rapid_threshold_select()`,
  `rapid_qi_contribution()`, `rapid_synthesizer_cv()` and `rapid_suite()`. `nobs()` and
  `format()` for `rapid` objects.
* `rapid_test()` adds the Wald test against the chance level, with
  `z = (R - pi0) / sqrt(pi0 (1 - pi0) / n)` and `pi0` the mean of the permutation distribution.
* `confint()` for `rapid` objects gains `method = "clopper-pearson"`.
* `rapid_test()`, `rapid_threshold_select()` and `rapid_suite()` accept a fitted `rapid` object
  and work on its stored predictions without a refit. `rapid_test()` and
  `rapid_synthesizer_cv()` gain a `plot()` method.
* All RAPID functions take the data as `X, Y, key_vars, target_var`. The old names
  `original_data`, `synthetic_data`, `quasi_identifiers` and `sensitive_attribute` are accepted
  as aliases. An argument whose name is a near miss of a `rapid()` argument is an error.
* Errors and warnings of the RAPID functions carry classes (`rapid_error_input`,
  `rapid_warning_deprecated` and others), so that they can be caught by class.
* Missing values of the target in `X`, `Y`, the holdout or the `newdata` of `predict()` give a
  warning of class `rapid_warning_missing_target` with their number per data set and their
  treatment. 0.2.0 replaced them without a message. Under the default `na_strategy =
  "constant"` the records of `X` stay in the denominator of RAPID and are scored against
  `na_constant_value`, so RAPID tends to be lower than on the observed records. Use
  `na_strategy = "drop"` to compute RAPID on the records with an observed target.
* `plot()` of a `rapid` object draws the gain rules (band, `RCS_marginal`) with ggplot2. Type 1
  is the histogram of the gains with the permutation null as an outline, type 2 the confidence
  `g` against the baseline `b`, and type 3 RAPID against the threshold with the chance level as
  a curve. Each ggplot type is returned invisibly, a list for several types.
* `disclosure_report()` gains `rapid_alpha` (see above) and stores the `rapid_test()` object in
  `results$rapid_test`.
* `rapid_qi_contribution()` measures the risk contribution of each quasi-identifier under a
  declared attack scenario. Each quasi-identifier (or a named group of substitutes, via
  `groups`) is removed from the attacker's knowledge, the attacker model is refitted, and the
  drop in RAPID is reported. Each scenario also reports its chance level, from `n_permutations`
  permutations of the true values (default 99) with the refitted attacker held fixed, and the
  drop in the excess over it. Comes with `print()`, `summary()`, `as.data.frame()` and a
  `plot()` dot plot.
* New vignette `vignette("rapid", package = "riskutility")`.

## Deprecated

* `rapid(return_all_records = )` is ignored with a warning. The result always holds all records.
  Use `as.data.frame(x, at_risk_only = TRUE)` for the flagged ones.
* `rapid(num_epsilon_type = )` is honoured with a warning. `"percentage"` becomes a proportion
  in `num_epsilon`, `"absolute"` becomes `num_error_metric = "absolute"`.
* `rapid_synthesizer_cv(return_details = , return_all_records = )` are ignored with a warning.
  The object always holds the fold table `folds` and the fold objects `fits`.

## Fixed

* `model_type = "gbm"` failed with xgboost 3.x (`ALTLIST classes must provide a Set_elt method`)
  for numeric and categorical targets.
* The class probabilities of a gbm attacker with more than two classes were mixed up between
  records and classes under xgboost 3.x. The error was hidden behind the previous one.
* A character target failed in `rapid()` (`non-numeric argument to binary operator`).
* An ordered-factor target failed (`comparison of these types is not implemented`).
* `confint(method = "bootstrap")` returned `[1, 1]` for an object fitted with
  `return_all_records = FALSE`, because it resampled the flagged records only.
* `rapid_synthesizer_cv()` stored `RCS_marginal` as the method of a categorical target while its
  folds used the `rapid()` default `RCS_conditional`. The stored method is now the one the folds
  use.
* `plot(x, which = c(1, 3))` drew only the last type. All requested types are drawn, in the
  order given.
* A class held by every original record gave a gain of `NaN`, counted as not at risk. The gain
  is 0 there, as in the paper.
* `rapid_test()` and `rapid_threshold_select()` on the data fitted the attacker twice and used
  the missing-value strategy `"constant"` for the permutations whatever the call said. They now
  fit once, through `rapid()`, and permute the fitted object's predictions.
* Classes of a categorical target that occur only in the synthetic data became `NA`. The levels
  are now the union, original levels first.
* The references of `?rapid` cited the R-U map paper. They now cite the RAPID paper and
  Meinshausen (2006).

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

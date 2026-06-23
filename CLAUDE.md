# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Test Commands

```bash
# Install dependencies
Rscript -e "devtools::install_deps(dependencies = TRUE)"

# Regenerate documentation (run after modifying roxygen comments)
Rscript -e "devtools::document()"

# Build package
R CMD build .

# Check package (CRAN-like)
_R_CHECK_FORCE_SUGGESTS_=FALSE R CMD check riskutility_*.tar.gz --no-manual --no-vignettes

# Run all tests
Rscript -e "devtools::test()"

# Run a single test file
Rscript -e "testthat::test_file('tests/testthat/test-dcap.R')"

# Load package for interactive development
Rscript -e "devtools::load_all()"
```

## Package Architecture

This R package measures **disclosure risk** and **data utility** for anonymized/synthetic data. The package supports multivariate Risk-Utility (R-U) evaluation following the framework described in "Beyond the Trade-off Curve" (Thees, Müller, Templ 2026).

### S3 Class Structure

All major functions return S3 objects with `print()`, `summary()`, and `plot()` methods:

#### Disclosure Risk Classes

| Class | Function | Purpose |
|-------|----------|---------|
| `dcap` | `dcap()` | Differential Correct Attribution Probability |
| `tcap` | `tcap()` | Targeted CAP (per-record risk) |
| `weap` | `weap()` | Within Equivalence Class Attribution Probability |
| `disco` | `disco()` | Disclosive in Synthetic Correct Original |
| `rapid` | `rapid()` | Risk of Attribute Prediction-Induced Disclosure (ML-based) |
| `dcr` | `dcr()` | Distance to Closest Record |
| `nndr` | `nndr()` | Nearest Neighbor Distance Ratio |
| `ims` | `ims()` | Identical Match Share |
| `recordLinkageRisk` | `recordLinkage()` | Record linkage risk (deterministic/probabilistic/PRAM/predictive/RF/RBRL/Mahalanobis/embedding; independent/bijective/OT matching) |
| `tcloseness` | `tcloseness()` | t-Closeness (EMD-based) |
| `ldiversity` | `ldiversity()` | l-Diversity (distinct/entropy/recursive) |
| `kanonymity` | `kanonymity()` | k-Anonymity assessment |
| `attacker_risk` | `attacker_risk()` | Prosecutor/journalist/marketer attacker models |
| `drisk` | `drisk()` | Distance-based risk (dRisk/dRiskRMD) |
| `population_uniqueness` | `population_uniqueness()` | Population uniqueness (Pitman/Zayatz/SNB) |
| `epsilon_identifiability` | `epsilon_identifiability()` | Epsilon-identifiability |
| `delta_presence` | `delta_presence()` | delta-Presence |
| `hitting_rate` | `hitting_rate()` | Hitting rate |
| `singling_out` | `singling_out()` | Singling out risk |
| `linkability` | `linkability()` | Linkability risk |
| `domias` | `domias()` | DOMIAS membership inference |
| `nnaa` | `nnaa()` | Nearest Neighbor Adversarial Accuracy |
| `disclosure_report` | `disclosure_report()` | Comprehensive risk report |

#### Data Utility Classes

| Class | Function | Purpose |
|-------|----------|---------|
| `propscore` | `propscore()` | Propensity score utility (pMSE) |
| `subgroup_utility` | `subgroup_utility()` | Stratified utility across subgroups |
| `regression_fidelity` | `regression_fidelity()` | Regression coefficient comparison |
| `contingency_fidelity` | `contingency_fidelity()` | Bivariate categorical dependence (TV distance) |
| `rumap` | `rumap()` | Multivariate R-U map |

Each class follows the pattern:
- Main function returns object with class attribute
- `print.classname()` - brief output
- `summary.classname()` - detailed statistics (returns `summary.classname`)
- `plot.classname()` - visualization (often with `which` parameter for multiple plot types)

### Three Risk Metric Families

1. **Attribution-based (CAP family)**: `dcap.R`, `tcap.R`, `weap.R`, `disco.R`
   - Measure if adversaries can infer sensitive attributes from quasi-identifiers
   - Require `key_vars` (quasi-identifiers) and `target_var` (sensitive attribute)
   - Use exact or fuzzy matching between original and synthetic records

2. **ML-based (RAPID)**: `rapid.R`
   - Uses machine learning models to predict sensitive attributes from quasi-identifiers
   - Trains on synthetic data, evaluates on original data
   - Captures complex non-linear relationships that CAP methods may miss
   - Supports: random forest (`rf`, default), linear model (`lm`), CART (`cart`), XGBoost (`gbm`), logistic (`logit`)
   - Requires optional packages: `ranger`, `rpart`, or `xgboost` depending on model choice
   - For numeric targets: `at_risk = TRUE` when prediction error < threshold (accurate prediction = disclosure risk)
   - For categorical targets: `at_risk = TRUE` when ratio/gain/score >= threshold
   - Default `cat_eval_method = "RCS_conditional"` with `cat_tau = 1`
   - Default `return_all_records = TRUE` ensures both at-risk and not-at-risk records are available for plotting
   - ranger defaults to `importance = "impurity"` so variable importance plots work out of the box
   - `store_model = TRUE` needed for `plot(result, which = 4)` (variable importance)
   - 6 plot types: histogram (1), scatter (2), threshold sensitivity (3), QI importance (4), QI attribution (5), QI interactions (6)
   - Attribution plots (5, 6) use sum-to-zero contrasts with reference-level recovery, adaptive quantile binning, and robust CI scaling
   - This is the canonical RAPID implementation; the research version lives in `~/workspace26/RAPID/`

3. **Distance-based (holdout method)**: `dcr.R`, `nndr.R`, `ims.R`
   - Detect memorization by comparing distances to training vs. holdout data
   - Use either explicit `holdout` parameter or `holdout_fraction` for automatic splitting
   - **Important**: See "DCR Delusion" warning below about limitations

### Utility Measures

Current utility functions for comparing original vs synthetic data:

| Function | Purpose | Data Type |
|----------|---------|-----------|
| `propscore()` | Propensity score / pMSE | Mixed |
| `compare_wasserstein()` | Wasserstein distance | Numeric |
| `compare_ks_test()` | Kolmogorov-Smirnov test | Numeric |
| `compare_chisq_gof()` | Chi-squared goodness of fit | Categorical |
| `compare_pca()` | PCA-based comparison | Numeric |
| `compare_embedding()` | t-SNE/UMAP embedding | Mixed |
| `compare_correlation_matrices()` | Correlation comparison | Numeric |
| `compare_feature_importance()` | Feature importance stability | Mixed |
| `compare_model_performance()` | Predictive model evaluation | Mixed |
| `gower()` | Average Gower distance | Mixed |
| `mqs()` | Multivariate quality score | Mixed |
| `ci_overlap()` | Confidence interval overlap | Numeric |

**Additional utility measures** (integrated with rumap):
- `hellinger()` - Hellinger distance for categorical distributions
- `energy_distance()` - Energy distance for multivariate numeric data
- `ci_proximity()` - CI proximity/relative error measure
- `repu()` - Replicated Uniques
- `mmd()` - Maximum Mean Discrepancy (kernel-based, exact + RFF)
- `tstr()` - Train on Synthetic, Test on Real (downstream ML utility)
- `copula_fidelity()` - Empirical copula dependence comparison
- `tail_fidelity()` - Tail preservation (QQ divergence + density ratio + Hill)
- `subgroup_utility()` - Stratified utility assessment across subgroups
- `regression_fidelity()` - Regression coefficient comparison (bias, CI overlap, significance)
- `contingency_fidelity()` - Bivariate contingency table TV distance (categorical complement to copula_fidelity)

### rumap() - Multivariate R-U Framework

The `rumap()` function provides comprehensive multivariate Risk-Utility evaluation:

```r
rumap(original, synthetic,
      risk_measures = c("dcap", "tcap", "disco", "rapid", "ims"),
      utility_measures = c("pmse", "wasserstein", "hellinger", "energy",
                           "ci_proximity", "mmd", "tstr", "copula", "tail",
                           "contingency"),
      key_vars = NULL, target_var = NULL,
      holdout = NULL, holdout_fraction = 0.2,
      normalize = TRUE, ...)
```

**Visualization types** via `plot.rumap(x, which = ...)`:
1. Composite scatterplot (classic R-U map with Pareto front)
2. Heatmap (all measures × methods)
3. Dot plot (individual measures with risk/utility facets)
4. Parallel coordinates (multivariate profiles)
5. Radial/Origami plot (polygonal profiles)
6. PCA biplot (joint PCA of all measures)
7. Blockwise PCA (PC1 utility vs PC1 risk)

**Key features:**
- Min-max normalization to [0,1] with harmonized directions
- Pareto frontier identification
- Internal consistency metrics (Cronbach's α, McDonald's ω)
- Support for comparing multiple synthetic data generators

### Key Dependencies

- **data.table**: Used extensively with NSE patterns (global variables declared in `zzz.R`)
- **ggplot2**: All plot methods
- **simPop**: Integration target for S4 class support
- **reshape2**: Use `reshape2::melt()` explicitly (not imported to avoid data.table conflict)
- **caret**: Model training for propensity scores
- **VIM**: Gower distance calculations
- **rrcov**: Robust PCA (for rumap diagnostics)
- **clue**: Hungarian algorithm for bijective record linkage (GDBRL)
- **robustbase**: MCD covariance estimation for Mahalanobis record linkage
- **torch**: Autoencoder-based embedding record linkage

### Important Files

- `R/zzz.R`: Package startup, global variables for NSE patterns
- `R/comparison-functions.R`: Package-level documentation
- `R/divergence.R` and `R/divergence2.R`: Entropy measures (KL, JS divergence)
- `R/evaluation_stats.R`: Error metrics (MAE, MSE, RMSE, MAPE)
- `R/propscore.R`: Propensity score utility
- `R/rapid.R`: RAPID metric with ML-based attribute inference risk assessment (~1940 lines; includes all S3 methods, internal plot helpers, and model fitting)
- `R/dcap.R`, `R/tcap.R`: CAP-based attribution metrics
- `R/dcr.R`, `R/nndr.R`, `R/ims.R`: Distance-based metrics
- `R/recordLinkage.R`: Record linkage risk (~3385 lines; 8 methods, independent/bijective/OT matching, S3 dispatch, blocking, direction; shared risk back-end `.true_match_risk()`, PRAM `.pram_expected_risk()`)
- `R/embedding_internal.R`: Torch autoencoder engine for embedding-based record linkage (~580 lines; `.ae_model`, `.ae_train`, `.ae_encode`, `.ae_distance`, `.ae_var_importance`)
- `R/sinkhorn_internal.R`: Sinkhorn optimal transport engine (`.sinkhorn()` and `.solve_ot()`)
- `R/mahalanobis_internal.R`: Robust Mahalanobis distance helpers (`.mahal_prepare`, `.mahal_dist`)
- `R/subgroup_utility.R`: Stratified utility assessment (~358 lines)
- `R/regression_fidelity.R`: Regression coefficient comparison (~472 lines; forest plot + CI overlap bar chart)
- `R/contingency_fidelity.R`: Bivariate contingency table fidelity (~423 lines; TV distance heatmap)
- `vignettes/riskutility.Rmd`: Comprehensive JSS-style vignette (main package vignette)
- `vignettes/recordLinkage.Rmd`: Deep-dive vignette on record linkage risk (~2240 lines; JSS paper candidate). **Build-excluded via `.Rbuildignore`** — not shipped to CRAN (sdcMicro name collision, see below); kept in the repo.
- `vignettes/references.bib`: Shared bibliography for all vignettes

### Record Linkage: Independent vs Bijective Matching

`recordLinkage()` supports three matching modes via `matching = c("independent", "bijective", "ot")`:

- **Independent** (default): Each record matched independently (many-to-one). Classical DBRL approach.
- **Bijective**: Global one-to-one assignment via Hungarian algorithm (`clue::solve_LSAP()`). Models GDBRL attacker (Herranz, Nin, Rodriguez & Tassa, 2016). Yields binary risk (0 or 1) and typically higher risk than independent matching.
- **OT** (optimal transport): Entropy-regularized soft global assignment via Sinkhorn-Knopp algorithm. Produces continuous risk in [0,1] that interpolates between independent (smooth) and bijective (hard). Controlled by `ot_epsilon` (regularization) and `ot_max_iter`.

Cost direction per method (bijective/OT cost transform):
- Deterministic/Predictive/Embedding: distance (minimize)
- Probabilistic: likelihood ratio (maximize → `max - LR` transformation)
- PRAM: transition probability (maximize → `max - prob` transformation)
- RF: proximity (maximize)

Key internals: score caching per record during main loop → `.solve_bijective()` or `.solve_ot()` builds cost matrix per block → `clue::solve_LSAP()` (bijective) or `.sinkhorn()` (OT) → override risk. Bijective adds `bijective_assigned` column; OT stores `transport_plans`. `clue` and `torch` are in Suggests.

### Record Linkage: Per-Record Risk Definition, NA Handling, Options

All 8 methods share one risk definition via the internal `.true_match_risk()`: the (weighted) probability mass on the **true** match within the attacker's candidate set (0 if the true record is not in the set), using `risk_weighting` (`uniform`/`softmax`/`kernel`) and the `strategy` candidate-set filter. The rf and embedding methods route their proximity / latent-distance scores through this same helper, so their `risk` is a re-identification probability comparable with the other methods; their former nearest-neighbour similarity is kept in the `nn_similarity` column of `per_record` as a diagnostic.

- **`na_anon`** (`ignore`/`match`/`mismatch`) is honored by **all** methods (Gower, `.fs_log_lr`, `.pram_risk`, RBRL, `.mahal_dist`, embedding deterministic fallback), not just the Gower path. Default `ignore`; NA-free data is byte-identical. PRAM `ignore` multiplies the transition factor by 1 instead of collapsing a missing-QI record to 0 risk. Mahalanobis numeric `ignore` is a zero-contribution approximation (documented).
- **`compute_baseline = TRUE`**: also links X against itself (forward, `truth="row"`) for a no-perturbation reference → `$baseline` (with `risk_reduction`); shown in print/summary and `plot(which=1)`.
- **`expected_risk = TRUE`** (PRAM, forward): perturbation-aware expected risk over the transition distribution `E[r_i] = sum_o P(x_i->o)^2/(S+P(x_i->o))` via `.pram_expected_risk()` (exact ≤ `max_support`, else Monte-Carlo) → `$pram_info$expected_risk`.
- User-supplied `m_probs`/`u_probs` are validated and clamped to (0,1).

### Record Linkage: sdcMicro Name Collision (CRAN Build)

`sdcMicro` (>= 5.8.2) **also exports its own `recordLinkage()`** (a GDBRL/Hungarian implementation) — a lasting name collision with `riskutility::recordLinkage()`. sdcMicro's signature is `recordLinkage(x, y, vars, distance, ...)` with **no `...`, `key`, or `method`**, so if it masks ours on the search path, calls fail with `unused argument (key = ...)` / `(method = ...)`. This broke the CRAN vignette rebuild: CRAN installs sdcMicro 5.8.2 (local dev had 5.8.1, so it never reproduced locally), `recordLinkage.Rmd` attached sdcMicro *after* riskutility, and vignettes rebuild in a shared R session so both vignettes failed.

**CRAN-build resolution (2026-06):** `vignettes/recordLinkage.Rmd` is `.Rbuildignore`d (not shipped; stays in the repo as the JSS candidate) and the live `recordLinkage()` demo in `riskutility.Rmd` is commented out. The `recordLinkage()` **function** (`R/recordLinkage.R`, man page, tests, export) is unchanged and still ships. As a backstop, the committed `recordLinkage.Rmd` loads `sdcMicro` *before* `riskutility`. Do NOT re-enable the recordLinkage vignette in the CRAN build without re-checking this collision (install sdcMicro 5.8.2 into a temp lib to reproduce).

### Embedding Method (Autoencoder)

When `method = "embedding"`, a torch autoencoder learns a joint latent space for mixed data using entity embeddings (Guo & Berkhahn, 2016) for categoricals. Distances are computed in latent space rather than on raw features.

Key internals in `R/embedding_internal.R`:
- `.ae_model()`: `nn_module` with entity embeddings (`min(50, floor(n_levels/2) + 1)` dims per categorical)
- `.ae_train()`: Manual training loop with composite loss (MSE + cross-entropy), early stopping
- `.ae_encode()`: Batch encoding to latent space
- `.ae_distance()`: Euclidean distances in latent space, normalized by 97.5th percentile
- `.ae_var_importance()`: Permutation-based variable importance via reconstruction loss delta
- `.embedding_linkage_block()`: Per-block orchestrator; falls back to deterministic for blocks < `max(30, 5 * latent_dim)` rows

### DCR Delusion Warning

Yao et al. (2025) demonstrate that DCR and related distance-based metrics can fail to detect privacy leakage:

- **False sense of security**: Datasets deemed "private" by DCR can still be vulnerable to Membership Inference Attacks (MIAs)
- **Null distribution matters**: DCR values must be compared against a proper null distribution
- **Not sufficient alone**: DCR should be used alongside other privacy metrics
- **Recommendation**: The `dcr()` function includes statistical tests and warnings, but users should be aware that passing DCR tests does not guarantee privacy protection

Reference: Yao, Z., Krco, N., Ganev, G., & de Montjoye, Y.-A. (2025). arXiv:2505.01524 "The DCR Delusion: Measuring the Privacy Risk of Synthetic Data"

## Git Conventions

- Always commit as `matthias-da`: use `git commit --author="matthias-da <matthias-da@users.noreply.github.com>"`
- Never commit as claude or any other identity

## Code Conventions

- All functions use roxygen2 documentation with `@export` tags
- Examples should use simple synthetic data (avoid external package datasets that may not be available)
- Wrap long-running examples in `\donttest{}`
- Use explicit namespace for functions that conflict (e.g., `reshape2::melt()`)
- Follow existing S3 class pattern: main function + print/summary/plot methods
- `summary()` methods should return a typed summary object (e.g., `summary.rapid`), with its own `print.summary.classname()` method
- `plot()` methods use integer `which` parameter for multiple plot types (e.g., `plot(x, which = 1:6)`)
- Internal helpers prefixed with `.` (e.g., `.fit_rapid_model`, `.plot_qi_importance`)
- Model fitting uses `modifyList(defaults, user_args)` pattern to allow user overrides via `...`
- Use ASCII text in plot labels (no `expression()` with Greek letters) for cross-platform compatibility
- Risk measures: higher values = higher risk (bad)
- Utility measures: higher values = higher utility (good), or transform accordingly

## RAPID Porting Notes

The RAPID implementation in riskutility is the canonical version, ported and improved from the research package at `~/workspace26/RAPID/`. Key differences from the research version:

- **Flat result object** (`result$rapid`, `result$threshold`) vs nested (`result$risk$confidence_rate`)
- **Proper summary pattern** (returns `summary.rapid` object) vs direct printing
- **6 plot types** (vs 3 in research version), with lollipop/forest plots instead of barplots
- **Attribution robustness**: adaptive quantile binning, reference-level recovery, robust CI scaling
- **Default importance**: ranger gets `importance = "impurity"` automatically via `modifyList` pattern
- **Not yet ported**: new error metrics (`mae`, `rmse`, `rmae`, `rrmse`) — the research package has an inconsistency between `rapid.R` and `evaluate_numeric.R`, so this was deferred

## References

- Thees, Müller, Templ (2026). "Beyond the Trade-off Curve: Multivariate and Advanced Risk-Utility Maps for Evaluating Anonymized and Synthetic Data." Journal of Official Statistics.
- Yao, Z., Krco, N., Ganev, G., & de Montjoye, Y.-A. (2025). "The DCR Delusion: Measuring the Privacy Risk of Synthetic Data." arXiv:2505.01524.
- Taub, J., et al. (2018). "Differential Correct Attribution Probability for Synthetic Data." Privacy in Statistical Databases.
- RAPID implementation based on Thees, Müller & Templ (2026) - Risk of Attribute Prediction-Induced Disclosure
- Herranz, J., Nin, J., Rodriguez, P., & Tassa, T. (2016). "Revisiting Distance-Based Record Linkage for Privacy-Preserving Release of Statistical Datasets." Data & Knowledge Engineering, 100, 78-93.
- Guo, C., & Berkhahn, F. (2016). "Entity Embeddings of Categorical Variables." arXiv:1604.06737.
- Karr, A. F., et al. (2006). "A Framework for Evaluating the Utility of Data Altered to Protect Confidentiality." The American Statistician, 60(3), 224-232.
- Snoke, J., et al. (2018). "General and Specific Utility Measures for Synthetic Data." JRSS-A, 181(3), 663-688.
- Related research code: `/Users/matthias/workspace25/PCA_RU/`, `/Users/matthias/workspace26/RAPID/`

## Package Status

- **CRAN**: On CRAN since 2026-06-22 (v0.1.0, first release) — <https://CRAN.R-project.org/package=riskutility>. Title: "Disclosure Risk and Data Utility Metrics for Synthetic and Anonymized Data".
- **Authors**: Matthias Templ (aut, cre); Oscar Thees (ctb)
- **Tests**: 69 test files, ~3900 tests passing
- **R CMD check**: `Status: OK` (0 errors / 0 warnings / 0 notes) with vignettes rebuilt
- **Vignettes**: Main vignette (`riskutility.Rmd`) ships to CRAN; record linkage deep-dive (`recordLinkage.Rmd`) is repo-only (build-excluded — see "Record Linkage: sdcMicro Name Collision")
- **CITATION**: Available in `inst/CITATION`


## Second brain pointers
- Project note: `~/SecondBrain/20-Projects/riskutility.md` — read at session start if it exists; create from `~/SecondBrain/90-Templates/project-note.md` on first write
- Decisions about this project: `grep -rl riskutility ~/SecondBrain/80-Decisions/`
- Related papers: `grep -rl riskutility ~/SecondBrain/30-Papers/`
- Master bibliography: `~/SecondBrain/_refs/library.bib`

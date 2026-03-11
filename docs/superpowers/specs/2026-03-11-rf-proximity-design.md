# RF Proximity Design Spec

**Date**: 2026-03-11
**Status**: Approved (rev 3 — post 6-agent review)
**Scope**: Three RF proximity use cases for the riskutility package

## Overview

Random Forest proximity — the fraction of trees where two observations land in the same terminal node — provides a data-adaptive similarity measure that handles mixed data types, captures non-linear relationships, and respects variable interactions.

### What is RF proximity?

A Random Forest classifies records by routing them through decision trees to terminal (leaf) nodes. Two records that frequently land in the same terminal node are "proximate" — the forest considers them similar based on the features it has learned. Formally, for records *i* and *j*:

    proximity(i, j) = (number of trees where i and j share a terminal node) / (total number of trees)

This produces a similarity score in [0, 1]. Unlike Gower or Euclidean distance, RF proximity is **data-adaptive**: the forest automatically learns which variables matter and captures non-linear relationships and interactions. It handles mixed data (numeric + categorical) natively without requiring explicit distance definitions.

### Three use cases

1. **Record linkage risk**: 5th method in `recordLinkage()` (`method = "rf"`)
2. **Memorization detection**: new `rf_privacy()` function (distance-risk family)
3. **Utility measurement**: new `method = "ranger"` in `propscore()`

All three share a common internal engine (`.rf_proximity()`) based on `ranger`.

### When to use RF proximity vs existing methods

| Scenario | Recommended method | Why |
|----------|-------------------|-----|
| Simple QIs, interpretable risk | `dcr()` with Gower | Transparent distance, fast, well-established |
| Mixed types with complex interactions | `rf_privacy()` | Adapts to variable importance and non-linear patterns |
| Many variables (20+) | `rf_privacy()` or `recordLinkage(method = "rf")` | RF handles high dimensions naturally; Gower gives equal weight to all |
| Speed on small data (n < 5,000) | `dcr()` | Faster, no forest training overhead |
| Structural utility beyond classification | `propscore(method = "ranger")` | Proximity-based structural diagnostics not available with `"rf"` or `"logreg"` |
| ML-savvy adversary model | `recordLinkage(method = "rf")` | Captures adversary who deploys ML-based attacks, not just QI comparison |
| Traditional QI-based adversary | `recordLinkage(method = "deterministic")` | More appropriate for adversaries who literally compare QI values |

Each `rf_privacy()`, `propscore(method = "ranger")`, and `recordLinkage(method = "rf")` help page will include an `@section When to use this method` block with this guidance. All three will cross-reference each other via `@seealso` and share a conceptual `@section RF proximity` block explaining the proximity definition.

### References

- Breiman, L. (2001). Random Forests. *Machine Learning*, 45(1), 5-32. (Original proximity definition)
- Lin, Y. & Jeon, Y. (2006). Random forests and adaptive nearest neighbors. *JASA*, 101(474), 578-590. (RF proximity as adaptive kernel)
- Shi, T. & Horvath, S. (2006). Unsupervised learning with random forest predictors. *JCGS*, 15(1), 118-138. (Unsupervised RF proximity)

These will be included as `@references` on the appropriate help pages.

## Implementation Order

1. `.rf_proximity()` internal engine + `.proximity_from_nodes()` helper
2. `.distance_risk_prepare()` refactor (extract from dcr/nndr)
3. `rf_privacy()` standalone function
4. `propscore(method = "ranger")` extension
5. `recordLinkage(method = "rf")` extension
6. Integration (disclosure_report, rumap) + tests

## Architecture Decision: Distance-Risk Family

**Decision**: Option (B) — shared internal engine, separate public functions.

Keep `dcr()`, `nndr()`, `ims()`, `rf_privacy()` as standalone functions. Extract shared holdout-splitting and validation logic into `.distance_risk_prepare()`.

**Rationale**:
- No breaking API change
- IMS needs no holdout — awkward in a unified API
- RF proximity uses a fundamentally different distance engine (terminal nodes vs Gower/Euclidean)
- Duplication is ~30 lines of boilerplate, easily extracted
- Better discoverability: each function gets its own help page

---

## Component 1: `.rf_proximity()` Internal Engine

**File**: `R/rf_proximity_internal.R` (dedicated file — the engine may grow to 200+ lines; keeps `utils_internal.R` lean)

**Roxygen**: Use `@keywords internal` for all helpers in this file (matching `utils_internal.R` and the majority of the codebase). This generates Rd files but hides them from the default index.

**Dependency guard**: `.rf_proximity()` calls `requireNamespace("ranger", quietly = TRUE)` with an informative error as a safety net. **Additionally**, each public entry point (`rf_privacy()`, `propscore()` ranger branch, `recordLinkage()` rf branch) adds a thin one-line `requireNamespace()` check with a caller-specific error message (e.g., `"Package 'ranger' required for rf_privacy(). Install with install.packages('ranger')"`). This matches the established pattern in `rapid.R` and `recordLinkage.R`. For `disclosure_report()`, use the existing skip-with-message pattern.

```r
.rf_proximity(data1, data2,
              vars = NULL,
              n_trees = 500,
              mtry = NULL,
              importance = TRUE,
              seed = NULL, ...)
```

**Returns**:
```r
list(
  forest         = ranger_object,
  terminal_nodes = matrix(n x n_trees, integer),  # combined data
  n1             = nrow(data1),
  n2             = nrow(data2),
  importance     = named_numeric_or_NULL,
  oob_error      = numeric
)
```

**Implementation**:
1. `requireNamespace("ranger", quietly = TRUE)` or stop with informative message
2. Validate `n_trees >= 10` (proximity from fewer trees is unreliable; minimum ensures meaningful fractions)
3. Check for `.rf_label` column name collision: `if (".rf_label" %in% names(data1) || ".rf_label" %in% names(data2)) stop("Column '.rf_label' already exists in data. Please rename it.")`
4. Combine `rbind(data1[, vars], data2[, vars])`, add `.rf_label = factor(c(rep(0, n1), rep(1, n2)))`
5. Handle high-cardinality unordered factors: if any factor has > 53 levels, pass `respect.unordered.factors = "partition"` to ranger (avoids ranger's fallback to order-based splitting). Emit `message()` noting this.
6. Train ranger using `modifyList` pattern (consistent with `rapid.R` and `recordLinkage.R`):
```r
default_args <- list(
  formula = .rf_label ~ .,
  data = combined,
  num.trees = n_trees,
  probability = TRUE,
  write.forest = TRUE,
  importance = if (importance) "impurity" else "none",
  mtry = mtry,
  seed = seed
)
args <- modifyList(default_args, list(...))
forest <- do.call(ranger::ranger, args)
```
7. Extract terminal nodes via `predict(forest, combined, type = "terminalNodes")$predictions`
8. Handle OOB prediction NAs for small datasets: if any OOB predictions are NA, warn and fall back to in-bag predictions for those records
9. Return forest, terminal nodes, metadata

**Key design choice**: Does NOT compute or store the full proximity matrix. Returns the terminal node matrix (O(n * n_trees) memory). Consumers compute what they need via `.proximity_from_nodes()`.

### `.proximity_from_nodes()` Helper

```r
.proximity_from_nodes(terminal_nodes, idx1, idx2)
# Returns: matrix of dim length(idx2) x length(idx1)
# Entry [i,j] = fraction of trees where idx2[i] and idx1[j] share a terminal node
```

Row-wise computation to avoid materializing the full n x n matrix when only a submatrix is needed.

**Computational cost and optimization**: O(length(idx2) * length(idx1) * n_trees). For `rf_privacy()` with n_syn = n_train = 10,000 and 500 trees this is 50 billion integer comparisons, expected to take 2-5 minutes on a modern laptop.

**Optimization strategy** (v1): Use a tree-by-tree approach rather than record-by-record. For each tree column, group records by terminal node ID using `data.table` or `match()`, then increment a proximity accumulator only for pairs sharing a node. This avoids the full O(n^2 * n_trees) scan when trees have many terminal nodes (each node has few records). If this is still too slow for n > 10,000, a C/Rcpp inner loop can be added as a follow-up.

The `progress` parameter in `rf_privacy()` is wired to this inner loop (one tick per tree or per row batch).

**Terminal node semantics**: `predict(forest, data, type = "terminalNodes")` uses in-bag predictions by default (each observation goes through all trees, including those where it was in-bag). This means training records' terminal nodes are partially overfit. This introduces a small conservative bias toward detecting memorization in the permutation null test. The help page documents this: "The permutation null test is slightly conservative because terminal node assignments for training records reflect partial in-bag overfitting."

### `.proximity_from_nodes_newdata()` Helper

For pushing holdout data through a trained forest:

```r
.proximity_from_nodes_newdata(forest, newdata, terminal_nodes_train, idx_train)
# Predicts terminal nodes for newdata, computes proximity to idx_train subset
# Returns: matrix of dim nrow(newdata) x length(idx_train)
```

Used by `rf_privacy()` to compute holdout proximity without retraining. Holdout records get OOB-style terminal nodes (they were never in-bag) so their proximity estimates are unbiased.

---

## Component 2: `.distance_risk_prepare()` Refactor

**File**: `R/utils_internal.R`

```r
.distance_risk_prepare(X, Y, holdout = NULL, holdout_fraction = 0.5,
                        vars = NULL, na.rm = TRUE, seed = NULL,
                        min_holdout = 1)
```

**Returns**:
```r
list(
  train     = data.frame,   # training portion of X (callers use prep$train)
  synthetic = data.frame,   # validated Y (callers use prep$synthetic)
  holdout   = data.frame,   # user-provided or split from X (callers use prep$holdout)
  vars      = character,    # validated variable names
  was_split = logical       # TRUE if holdout was auto-created from X
)
```

**Handles**:
- Validate X, Y are data frames
- Validate `holdout_fraction` is in (0, 1) exclusive
- Intersect variable names across X, Y, holdout
- Apply `vars` filter
- Remove NAs if `na.rm = TRUE`
- Split X into train/holdout if `holdout` is NULL (using `seed`). When `seed` is non-NULL, calls `set.seed(seed)` — this modifies the global `.Random.seed`. Documented in the function's roxygen.
- Ensure holdout has at least `min_holdout` records after splitting (nndr needs `min_holdout = 2` for nearest-neighbor ratio; dcr and rf_privacy use `min_holdout = 1`)
- Validate holdout has matching columns

**Note on null_test logic**: The null test permutation stays in the callers (`dcr()`, `rf_privacy()`), not in this shared helper. `nndr()` does not have null_test capability, and the permutation logic differs between dcr (re-compute distances) and rf_privacy (re-partition column indices).

**Refactor**: Replace boilerplate in `dcr()` and `nndr()` with call to this helper. Internal variable names in dcr/nndr are mapped to the returned list fields (`prep$train`, `prep$synthetic`, `prep$holdout`). Run existing test suites to verify no regressions.

---

## Component 3: `rf_privacy()`

**File**: `R/rf_privacy.R`
**S3 class**: `"rf_privacy"`
**Family**: `distance-risk`

### S3 dispatch

Follows the standard package pattern:
```r
rf_privacy <- function(X, ...) UseMethod("rf_privacy")
rf_privacy.synth_pair <- function(X, ...) { ... }
rf_privacy.default <- function(X, Y, ...) { ... }
```

All S3 methods (`print.rf_privacy`, `summary.rf_privacy`, `print.summary.rf_privacy`, `plot.rf_privacy`) get `@export` tags.

### Signature

```r
rf_privacy.default(X, Y,
                   holdout = NULL, holdout_fraction = 0.5,
                   vars = NULL, na.rm = FALSE, seed = NULL,
                   progress = FALSE,
                   null_test = TRUE, n_null = 100,
                   n_trees = 500, mtry = NULL, ...)
```

**Parameter order**: matches `dcr.default()` — shared params first (`X, Y, holdout, holdout_fraction, vars`), then `na.rm`, `seed`, `progress`, `null_test`, `n_null`, then RF-specific params (`n_trees`, `mtry`) at the end.

**`na.rm = FALSE` default** (differs from dcr/nndr which use `TRUE`): ranger handles missing values natively via surrogate splits. Stripping NAs is aggressive for real-world datasets with common missingness (e.g., 20% missing in one variable would lose 20% of records, biasing the risk assessment). When `na.rm = FALSE`, ranger handles NAs internally. When `na.rm = TRUE`, records with any NA are removed before training (matching dcr/nndr behavior). The help page documents this difference and recommends `na.rm = FALSE` (the default) for most use cases.

**Seed handling**: `seed` is passed to `.distance_risk_prepare()` for holdout splitting (sets global `.Random.seed`). For the forest, `seed + 1L` is passed to `.rf_proximity()` → ranger's internal seed. This avoids correlation between the holdout split and forest randomization.

### Why a supervised (real-vs-synthetic) forest?

The RF is trained to discriminate training data from synthetic data. This is intentional: the forest learns the boundary between the two distributions. A synthetic record that is memorized (copied from training) will land in terminal nodes dominated by training records, yielding high proximity to training and low proximity to holdout.

Intuitively: if no memorization occurred, synthetic records should land in terminal nodes with roughly equal mixtures of training and holdout observations. If a synthetic record consistently lands in nodes dominated by training records, it may have been copied from the training set.

**Limitation**: if the synthetic data differs systematically from the real data (different marginals, shifted distributions), the forest captures distributional differences, not just memorization. For example, if the synthetic data was generated with a shifted mean, synthetic records will show high proximity to whichever real subset has a similar mean, regardless of whether individual records were copied.

**Mitigation**: the holdout comparison controls for this. Both training and holdout are real data, so distributional differences affect both equally. If memorization is absent, synthetic records should have equal proximity to training and holdout. The `null_test` permutation further calibrates against the null.

**Power inversion**: the privacy assessment is most informative when OOB error is near 0.5 (high-quality synthetic data that the forest cannot easily distinguish from real data). When OOB error is low (poor utility / easy discrimination), the forest has strong discriminative signal but this reflects distributional differences, not necessarily memorization. When OOB error is near 0.5, the forest has weak signal but any remaining asymmetry between training and holdout proximity is more likely to reflect actual memorization. The result reports `oob_error` prominently, and the help page guides interpretation: "When `oob_error > 0.45`, the forest has little discriminative power; interpret memorization results with caution."

The help page includes an `@section Comparison with DCR` explaining that `rf_privacy()` is the RF-proximity analog of `dcr()`, using the same holdout design but replacing Gower/Euclidean distance with terminal-node co-occurrence. It also notes that unsupervised RF proximity (Shi & Horvath 2006) is a potential alternative (parking lot).

### Flow

1. `requireNamespace("ranger")` check with caller-specific message
2. `.distance_risk_prepare(X, Y, holdout, holdout_fraction, vars, na.rm, seed, min_holdout = 1)`
3. `.rf_proximity(train, synthetic, vars, n_trees, mtry, importance = TRUE, seed = seed + 1L)`
4. Per synthetic record: compute **max** and **mean** proximity to all training records
5. Push holdout through forest: `predict(forest, holdout, type = "terminalNodes")`
6. Per synthetic record: compute **max** and **mean** proximity to all holdout records
7. Compute `prox_share`, `prox_ratio`, `max_prox_share`, `max_prox_ratio`, Wilcoxon test
8. If `null_test`: permute train/holdout `n_null` times → null distribution
9. Determine `privacy_pass`

### Metric definitions

Proximity is a **similarity** measure (higher = more similar), opposite to distance.

**Max-based metrics** (primary — analogous to DCR's nearest-neighbor approach):
- **`max_prox_share`**: `mean(max_prox_train > max_prox_holdout) + 0.5 * mean(max_prox_train == max_prox_holdout)` — fraction of synthetic records whose maximum proximity to training exceeds maximum proximity to holdout. Uses mid-rank correction for ties (RF proximity values are rational numbers k/n_trees, so ties are structurally more common than with continuous Gower distances). ~0.5 = no memorization, >0.5 = memorization signal.
- **`max_prox_ratio`**: `mean(max_prox_train) / mean(max_prox_holdout)` — ratio of average max-proximities. ~1 = no memorization, >1 = memorization signal. Guard: if `mean(max_prox_holdout) < 1/n_trees`, return `NA_real_` with warning (degenerate case where holdout records never co-terminate with synthetic records).
- **`max_prox_train`**: per-record max proximity to nearest training record
- **`max_prox_holdout`**: per-record max proximity to nearest holdout record

**Mean-based metrics** (supplementary — captures aggregate distributional similarity):
- **`prox_share`**: `mean(prox_train_mean > prox_holdout_mean) + 0.5 * mean(prox_train_mean == prox_holdout_mean)` — with mid-rank tie correction.
- **`prox_ratio`**: `mean(prox_train_mean) / mean(prox_holdout_mean)` — same zero-guard as above.
- **`prox_train_mean`**: per-record mean proximity to all training records
- **`prox_holdout_mean`**: per-record mean proximity to all holdout records

The help page explains when to use which: "The max-based metrics detect individual memorized records (a single synthetic record very close to a specific training record). The mean-based metrics detect aggregate distributional leakage (systematic over-similarity to training data). Use `max_prox_share` as the primary privacy indicator."

**Relationship between max and mean metrics**: They can disagree. A few highly memorized records can drive `max_prox_ratio` well above 1 while `prox_ratio` stays near 1 (most records are fine). Conversely, a systematic distributional shift can push `prox_share` above 0.5 while `max_prox_share` stays near 0.5 (no specific record is copied, but synthetic data is globally too similar to training). Both scenarios are documented in the help page.

### `privacy_pass` threshold

When `null_test = TRUE` (the default): `privacy_pass` is derived from the permutation null distribution. `privacy_pass = TRUE` when the observed `max_prox_share` falls within the 95th percentile of the null distribution AND `max_prox_ratio` does not exceed the 95th percentile of the null ratio distribution. The 0.55 heuristic is not used.

When `null_test = FALSE`: `privacy_pass = max_prox_share <= 0.55 & wilcox_p > 0.05`. The 0.55 threshold is a screening heuristic carried over from DCR. Proximity-based share may have a different null distribution than distance-based share, so this is a rough approximation.

Permutation p-value uses `(sum(null >= observed) + 1) / (n_null + 1)` following Phipson & Smyth (2010) to avoid p-value = 0.

### Wilcoxon test

Uses `wilcox.test(max_prox_train, max_prox_holdout, paired = TRUE, alternative = "greater")`. Direction is `"greater"` (memorization means higher proximity to training, opposite of DCR which uses `"less"` for distance).

**Caveat**: the Wilcoxon signed-rank test assumes paired differences are independently drawn from a symmetric distribution. This independence assumption is violated because all proximity values share the same trained forest. The p-value is **anti-conservative** (too many small p-values). The help page documents this: "The Wilcoxon p-value is a heuristic indicator. The permutation null test (`null_test = TRUE`, the default) provides a principled alternative that accounts for the shared forest structure." When both are available, prioritize the null test.

### Permutation null test

Permutation scheme: shuffle which real records are labeled "train" vs "holdout" and recompute `max_prox_share` (and `prox_share`, `max_prox_ratio`). The forest is NOT retrained — the terminal node structure is fixed. This tests the conditional null: given the forest, are synthetic records' proximities symmetric between any random partition of real records?

**Efficiency**: Pre-compute the full synthetic-vs-all-real proximity vectors once (both max and mean per synthetic record to each real record). Store as matrix `n_syn x n_real`. Then each permutation only re-partitions columns into train/holdout and takes column-wise max/mean — O(n_null * n_syn * n_real) additions, trivial compared to the initial proximity computation.

**`n_null = 100`**: adequate for 0.05 significance screening. Users wanting precise p-values should increase `n_null` (documented in help page).

### Result fields

- `max_prox_share`: primary metric — fraction with higher max proximity to training (with tie correction)
- `max_prox_ratio`: ratio of mean max-proximities (train/holdout)
- `max_prox_train`: per-record max proximity to nearest training record
- `max_prox_holdout`: per-record max proximity to nearest holdout record
- `prox_share`: supplementary — fraction with higher mean proximity to training (with tie correction)
- `prox_ratio`: ratio of mean mean-proximities (train/holdout)
- `prox_train_mean`: per-record mean proximity to training
- `prox_holdout_mean`: per-record mean proximity to holdout
- `privacy_pass`: logical (null-test-derived when available, heuristic otherwise)
- `wilcox_test`: Wilcoxon test object
- `null_distribution`: null stats including null shares, null ratios, p-values (if `null_test = TRUE`)
- `oob_error`: OOB classification error from the forest
- `var_importance`: named numeric vector from the RF (uses `var_importance` to match `recordLinkage`)
- `n_synthetic`, `n_train`, `n_holdout`, `vars`

### S3 methods

- `print.rf_privacy`: one-line summary showing `max_prox_share`, pass/fail, and `oob_error`
- `summary.rf_privacy`: detailed statistics (both max and mean metrics, null distribution, Wilcoxon, per-record outliers)
- `print.summary.rf_privacy`: formatted summary output

**Mock print output**:
```
RF Privacy Assessment (rf_privacy)
  Max proximity share:  0.52 (training not preferred)
  Max proximity ratio:  1.03
  OOB error:            0.47
  Null test:            PASS (observed within 95% null interval, p = 0.38)
  Privacy:              PASS
```

- `plot.rf_privacy` with `which`:
  1. Paired density (max-proximity-to-training vs max-proximity-to-holdout)
  2. Per-record max proximity difference histogram (with flagged outliers)
  3. Null distribution with observed value (if `null_test = TRUE`)

### Help page structure

- `@details`: conceptual explanation of RF proximity for SDC (self-contained, no assumed RF knowledge)
- `@section When to use this method`: decision guide (from Overview table)
- `@section Interpretation`: metric interpretation guide for max and mean metrics
- `@section Comparison with DCR`: relationship to `dcr()`, advantages and limitations
- `@section Limitations`: supervised forest conflation, power inversion, Wilcoxon anti-conservatism
- `@section Computational considerations`: expected runtime for typical sizes (~5 seconds for n=1,000; 2-5 minutes for n=10,000)
- `@seealso`: `dcr`, `nndr`, `ims`, `propscore`, `recordLinkage`
- `@references`: Breiman (2001), Lin & Jeon (2006)
- `@examples`: two examples — memorized data (copy some training records) and random synthetic data

---

## Component 4: `propscore(method = "ranger")`

**File**: `R/propscore.R` (extend existing)

### Method naming

The existing `method = "rf"` uses `randomForest::randomForest()`. The new method uses `ranger::ranger()`. Both are random forests, but "ranger" is chosen because "rf" is taken. The distinction is:

- `method = "rf"`: standard RF propensity scores via randomForest (existing, backward-compatible)
- `method = "ranger"`: RF propensity scores via ranger + optional proximity-based structural utility analysis (new). Use when you want to understand not just whether records can be classified as original or synthetic, but also how structurally similar the two datasets are in terms of which records cluster together.
- `method = "logreg"`: logistic regression propensity scores (existing)

The help page will prominently state: "This is a random forest method (via the ranger package) that extends `method = 'rf'` with proximity-based structural diagnostics."

**Add `match.arg()`**: The existing `propscore()` lacks `match.arg(method)`. As part of this change, add `method <- match.arg(method)` with valid values `c("rf", "ranger", "logreg")` to prevent silent fallthrough.

### Parameter additions

```r
propscore(X, Y, form = NULL,
          method = c("rf", "ranger", "logreg"),   # match.arg, default "rf"
          proximity = c("summary", "full", "none"),  # new, ranger only
          importance = TRUE,                          # new, ranger only
          adjust_size = TRUE, cluster = NULL, na = "impute", ...)
```

When `proximity` or `importance` are explicitly passed with `method != "ranger"`, emit a `message()` noting they are only used with `method = "ranger"`. Silent ignoring only when they are at their defaults.

The `synth_pair` method passes `proximity` and `importance` through `...`. Test coverage includes `synth_pair` dispatch with these new params.

### New branch

(alongside existing `"logreg"` and `"rf"` branches):
1. `requireNamespace("ranger")` check with `"Package 'ranger' required for propscore(method = 'ranger')"`
2. Call `.rf_proximity()` with the combined labeled data (using `modifyList` pattern)
3. Propensity scores from `forest$predictions[, 2]` (probability mode, OOB predictions). Handle NA OOB predictions for small datasets.
4. pMSE via existing formula: `1/(2*n) * sum((p - cr)^2)`
5. If `proximity == "summary"`: compute structural metrics on-the-fly from terminal node matrix
6. If `proximity == "full"`: additionally store full n x n matrix
7. If `proximity == "none"`: skip proximity computation

### Additional result fields

- `oob_error`: OOB classification error (~0.5 = good utility). When `oob_error > 0.45`, emit `message()` noting the forest has little discriminative power.
- `var_importance`: named numeric vector (when `importance = TRUE`). Uses `var_importance` to match `recordLinkage` convention.
- When `proximity != "none"`:
  - `within_orig_prox`: mean proximity among original records
  - `within_synth_prox`: mean proximity among synthetic records
  - `cross_prox`: mean proximity between original and synthetic
  - `structure_ratio`: `cross_prox / ((within_orig_prox + within_synth_prox) / 2)` — closer to 1 = better utility. Uses arithmetic mean of the two within-class proximities as denominator. This is a descriptive index, not a test statistic. Range is typically [0, 1], where values near 1 indicate that original and synthetic records are equally intermixed in the forest's learned structure. Values well below 1 suggest the forest can separate the distributions. Typical values for good synthetic data: 0.8-1.0; below 0.5 indicates substantial distributional divergence. Note: values > 1 are theoretically possible when OOB error is ~0.5 (forest cannot discriminate). Document the expected range and interpretation in the help page.
- When `proximity == "full"`:
  - `proximity_matrix`: full n x n matrix

### Scaling

`warning()` (not just `message()`) when n > 10,000 with `proximity = "full"`, reporting estimated memory in MB and suggesting `proximity = "summary"`. At n = 10,000 combined, the full matrix is ~800 MB; at n = 20,000 it is ~3 GB.

### Plot

Existing: `which = 1` (density), `which = 2` (density ratio). New:
- `which = 3`: proximity structure boxplot (within-original, within-synthetic, cross-class distributions). Only available when `proximity != "none"`.
- `which = 4`: variable importance bar chart (when `importance = TRUE`).

---

## Component 5: `recordLinkage(method = "rf")`

**File**: `R/recordLinkage.R` (extend existing)

### Parameter additions

```r
recordLinkage(...,
              method = c("deterministic", "probabilistic", "pram", "predictive", "rf"),
              n_trees = 500,       # only for rf
              rf_global = FALSE,   # only for rf
              ...)
```

`n_trees` and `rf_global` silently ignored for other methods (these are new RF-specific parameters with no pre-existing user expectation; silent ignoring is obvious).

### Parameter handling for existing parameters

Parameters meaningful for other methods but not for RF:

- **`key`**: maps to `vars` in `.rf_proximity()`. The RF trains on `key` variables as features.
- **`type`** (ordinal detection): ignored for `method = "rf"`. RF handles ordinal/nominal natively.
- **`weights`**: ignored for `method = "rf"`. RF determines variable importance internally. A `message()` is emitted if `weights` is explicitly provided with `method = "rf"`, noting that RF uses data-driven importance instead.
- **`strategy`**: The existing strategies are `"nearest"`, `"threshold"`, `"topk"`, `"topk_threshold"`, `"nearest_threshold"`. For RF, `"nearest"` means highest proximity. The `"threshold"` and `"topk"` variants can also work (proximity as score). `"nearest_threshold"` and `"topk_threshold"` require a score threshold — documented that proximity thresholds in [0, 1] differ from Gower distance thresholds. `message()` if `strategy != "nearest"` with `method = "rf"`, noting that thresholds are proximity-based (not distance-based).
- **`risk_weighting`**, **`kappa`**, **`bandwidth`**, **`kernel`**: ignored for `method = "rf"`.
- **`direction`**: supported. Controls which dataset is query vs search. `"forward"` (default): synthetic records are matched to original. `"reverse"`: original records matched to synthetic. The RF always trains on `rbind(original, synthetic)` regardless of direction; direction only determines which records are queries and which are candidates in the proximity-based matching step.
- **`truth`**: both `"row"` and `"id"` modes supported, same as other methods.

The help page includes a `@section RF method` block that explicitly lists which parameters are used vs. ignored, following the pattern of `@section PRAM method` in the existing help page.

### Dependency guard

`requireNamespace("ranger")` check with `"Package 'ranger' required for recordLinkage(method = 'rf')"` at the top of the rf branch.

### Variable importance

Returned by default in `result$var_importance` (matching the existing `recordLinkage` field name). Shows which QIs drive linkage risk.

### Matching logic

**Independent matching**: each query record matched to candidate with highest proximity.

**Bijective matching**: cost matrix = `1 - cross_proximity` → `clue::solve_LSAP()`. Binary risk (1 if matched to true record, 0 otherwise). Note: proximity values from a supervised forest cluster near 0 and 1, so the cost matrix may have low variance. This can make the Hungarian algorithm's solution sensitive to forest hyperparameters. Document this in the help page.

### Blocking modes

**`rf_global = FALSE` (default)**: train separate RF per block.
- Each forest learns within-block patterns only
- More precise, avoids cross-block contamination
- **Small block handling**: blocks with < 2 records in either class (original or synthetic) fall back to `rf_global = TRUE` for that block with a `message()`. Blocks with 2+ records per class but < 50 total records are trained normally but flagged.
- **Message aggregation**: instead of one message per small block, emit a single summary: `"17 of 42 blocks have < 50 records. Consider rf_global = TRUE or coarser blocking (e.g., block on fewer variables)."`
- Per-block variable importance aggregated (weighted by block size) into `result$var_importance`

**`rf_global = TRUE`**: train one global RF.
- Extract terminal nodes for all records
- Restrict proximity to within-block pairs
- Robust to small blocks (the forest trains on all data, within-block proximity just restricts matching)
- Single variable importance from the global forest
- Note: global proximity reflects cross-block patterns, so a record in block A may get its proximity score influenced by block B patterns

Help page documents the trade-off: "Per-block forests learn block-specific patterns and are more precise, but require sufficient block size (~50+ records per block, i.e., ~25 per class). Use `rf_global = TRUE` when most blocks are small, or consider coarser blocking (e.g., block on region alone instead of region x sex x age-group). The 50-record threshold is a heuristic; datasets with many variables may need larger blocks."

**No blocking**: `message()` when n > 10,000, suggesting blocking and reporting estimated computation time.

### Helper function

The RF branch logic (especially per-block training) is factored into `.rf_linkage_block()` to keep the main `recordLinkage()` function readable (it is already ~2230 lines).

### Result

Same structure as other methods:
- `per_record`: risk, cand_n, true_in_set, bijective_assigned (for bijective)
- `var_importance`: named numeric (matching existing field name)
- Score interpretation: proximity 0–1, higher = riskier match

RF importance is displayed via existing `which = 3` ("Per-variable importance") in `plot.recordLinkageRisk`. No new plot type needed — the existing infrastructure handles it.

---

## Scaling Strategy

Informational messages for large datasets, with appropriate severity:
- n > 10,000 with `proximity = "full"` in propscore → `warning()` reporting estimated memory, suggesting `"summary"`
- n > 10,000 without blocking in recordLinkage → `message()` suggesting blocking with estimated compute time
- n > 10,000 in rf_privacy → `message()` noting expected runtime (2-5 minutes)
- Uses `message()` for informational guidance, `warning()` for memory-intensive operations
- No approximate algorithms in v1

**Expected runtimes** (n_trees = 500, modern laptop, documented in help pages):
- n = 1,000: ~5 seconds
- n = 5,000: ~30 seconds
- n = 10,000: 2-5 minutes
- n = 50,000: 30+ minutes (suggest blocking or subsampling)

---

## Dependencies

- `ranger`: already in Suggests (used by RAPID). No new hard dependencies.
- `clue`: already in Suggests (used by bijective matching). Needed for `recordLinkage(method = "rf", matching = "bijective")`.

---

## Test Plan

| File | Tests |
|------|-------|
| `test-rf-proximity-internal.R` | `.rf_proximity()`: correct dimensions, proximity in [0,1], symmetric, identical data → high proximity, `n_trees >= 10` validation, `.rf_label` column collision, high-cardinality factor handling, `modifyList` override. `.proximity_from_nodes()`: submatrix extraction, tie handling. `.proximity_from_nodes_newdata()`: holdout push-through. |
| `test-distance-risk-prepare.R` | `.distance_risk_prepare()`: holdout splitting, variable intersection, NA removal, `holdout_fraction` in (0,1) validation, `min_holdout` enforcement. Regression: `dcr()` and `nndr()` existing tests still pass after refactor. |
| `test-rf-privacy.R` | Memorized data → high `max_prox_share`, random data → ~0.5, `privacy_pass` logic (both null-test-derived and heuristic fallback), null_test permutation, tie correction, `prox_ratio` zero-guard, `na.rm = FALSE` with missing data, seed separation, OOB error reporting, S3 methods (print/summary/plot), `synth_pair` dispatch. |
| `test-propscore-ranger.R` | Identical distributions → OOB ~0.5, different distributions → OOB < 0.5, `proximity = "summary"/"full"/"none"`, `structure_ratio`, `var_importance`, `warning()` for large n + full, `match.arg` validation, `message()` for ranger params with non-ranger method, `synth_pair` dispatch with `proximity`/`importance`. New plot types (`which = 3` boxplot, `which = 4` importance). |
| `test-recordLinkage-rf.R` | Independent + bijective matching, with/without blocking, `rf_global` TRUE/FALSE, `var_importance`, small block fallback (<2 per class), aggregated small block message, large n message, direction forward/reverse, truth = "row"/"id", `weights` message, `strategy` handling, existing `which = 3` plot with RF importance. |

---

## Integration

- `disclosure_report()`: add `rf_privacy` as optional metric. Use skip-with-message pattern when ranger is not installed.
- `rumap()`: add `"rf_privacy"` to valid risk measures
- Both already discover metrics dynamically — minimal wiring

---

## Vignette

Add a section to the main `riskutility.Rmd` vignette covering:
- Brief RF proximity explanation (2-3 paragraphs)
- `rf_privacy()` worked example: same dataset evaluated with both `dcr()` and `rf_privacy()`, showing where results agree and disagree
- `propscore(method = "ranger")` worked example showing structural utility metrics
- Comparison table: RF proximity vs Gower distance vs Euclidean (handles mixed types, captures interactions, requires holdout, computation cost, established in literature)

---

## Reproducibility Note

Ranger's terminal node assignments can vary across package versions (due to internal tie-breaking or split-point selection changes). The help page notes: "Results depend on the ranger version and the random seed. Pin the ranger version for reproducible risk assessments."

---

## Parking Lot (future work)

- Approximate proximity for very large n (subsampled proximity, sparse terminal node matching)
- Unsupervised RF proximity (Shi & Horvath 2006, Breiman-Cutler) as alternative to supervised real-vs-synthetic forest — addresses limitation of conflating distributional similarity with memorization
- C/Rcpp inner loop for `.proximity_from_nodes()` if tree-by-tree optimization is insufficient
- Comparison vignette: RF proximity vs Gower vs Euclidean for record linkage (deeper than vignette section)
- ranger hyperparameter tuning beyond n_trees/mtry (expose min.node.size?)
- Geometric mean alternative for `structure_ratio` denominator when within-class proximities differ substantially
- Log-transform of bijective cost matrix (`-log(proximity)`) for better LSAP stability
- `rf_global = "auto"` option that switches to global when >25% of blocks fall below size threshold

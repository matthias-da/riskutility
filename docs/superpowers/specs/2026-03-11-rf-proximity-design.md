# RF Proximity Design Spec

**Date**: 2026-03-11
**Status**: Approved (rev 2 — post spec review)
**Scope**: Three RF proximity use cases for the riskutility package

## Overview

Random Forest proximity — the fraction of trees where two observations land in the same terminal node — provides a data-adaptive similarity measure that handles mixed data types, captures non-linear relationships, and respects variable interactions. This spec covers three use cases:

1. **Record linkage risk**: 5th method in `recordLinkage()` (`method = "rf"`)
2. **Memorization detection**: new `rf_privacy()` function (distance-risk family)
3. **Utility measurement**: new `method = "ranger"` in `propscore()`

All three share a common internal engine (`.rf_proximity()`) based on `ranger`.

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

**Dependency guard**: `.rf_proximity()` calls `requireNamespace("ranger", quietly = TRUE)` with an informative error. This is the single guard point — public callers (`rf_privacy`, `propscore`, `recordLinkage`) do not need their own checks.

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
2. Combine `rbind(data1[, vars], data2[, vars])`, add `.rf_label = factor(c(rep(0, n1), rep(1, n2)))`
3. Train `ranger::ranger(.rf_label ~ ., data = combined, num.trees = n_trees, probability = TRUE, write.forest = TRUE, importance = if(importance) "impurity" else "none", mtry = mtry, seed = seed, ...)`
4. Extract terminal nodes via `predict(forest, combined, type = "terminalNodes")$predictions`
5. Return forest, terminal nodes, metadata

**Key design choice**: Does NOT compute or store the full proximity matrix. Returns the terminal node matrix (O(n * n_trees) memory). Consumers compute what they need via `.proximity_from_nodes()`.

### `.proximity_from_nodes()` Helper

```r
.proximity_from_nodes(terminal_nodes, idx1, idx2)
# Returns: matrix of dim length(idx2) x length(idx1)
# Entry [i,j] = fraction of trees where idx2[i] and idx1[j] share a terminal node
```

Row-wise computation to avoid materializing the full n x n matrix when only a submatrix is needed.

**Computational cost**: O(length(idx2) * length(idx1) * n_trees). For `rf_privacy()` with n_syn = n_train = 10,000 and 500 trees this is 50 billion comparisons. In practice, the comparison is vectorized per-row: for each synthetic record, `colSums(nodes[idx1, ] == nodes_i)` across all trees. This is O(n_train * n_trees) per synthetic record, feasible but not cheap. Approximate subsampling is deferred to the parking lot.

### `.proximity_from_nodes_newdata()` Helper

For pushing holdout data through a trained forest:

```r
.proximity_from_nodes_newdata(forest, newdata, terminal_nodes_train, idx_train)
# Predicts terminal nodes for newdata, computes proximity to idx_train subset
# Returns: matrix of dim nrow(newdata) x length(idx_train)
```

Used by `rf_privacy()` to compute holdout proximity without retraining.

---

## Component 2: `.distance_risk_prepare()` Refactor

**File**: `R/utils_internal.R`

```r
.distance_risk_prepare(X, Y, holdout = NULL, holdout_fraction = 0.5,
                        vars = NULL, na.rm = TRUE, seed = NULL)
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
- Intersect variable names across X, Y, holdout
- Apply `vars` filter
- Remove NAs if `na.rm = TRUE`
- Split X into train/holdout if `holdout` is NULL (using `seed`)
- Validate holdout has matching columns

**Refactor**: Replace boilerplate in `dcr()` and `nndr()` with call to this helper. Internal variable names in dcr/nndr are mapped to the returned list fields (`prep$train`, `prep$synthetic`, `prep$holdout`). Run existing test suites to verify no regressions.

---

## Component 3: `rf_privacy()`

**File**: `R/rf_privacy.R`
**S3 class**: `"rf_privacy"`
**Family**: `distance-risk`

```r
rf_privacy(X, Y,
           holdout = NULL, holdout_fraction = 0.5,
           vars = NULL, na.rm = TRUE, seed = NULL,
           n_trees = 500, mtry = NULL,
           null_test = TRUE, n_null = 100,
           progress = FALSE, ...)
```

### Why a supervised (real-vs-synthetic) forest?

The RF is trained to discriminate training data from synthetic data. This is intentional: the forest learns the boundary between the two distributions. A synthetic record that is memorized (copied from training) will land in terminal nodes dominated by training records, yielding high proximity to training and low proximity to holdout.

**Limitation**: if the synthetic data differs systematically from the real data (different marginals, shifted distributions), the forest captures distributional differences, not just memorization. The proximity metric then conflates "distributional similarity" with "memorization." This is the same limitation as DCR with Gower distance — a synthetic record can be close to training simply because it's realistic.

**Mitigation**: the holdout comparison controls for this. Both training and holdout are real data, so distributional differences affect both equally. If memorization is absent, synthetic records should have equal proximity to training and holdout. The `null_test` permutation further calibrates against the null.

The help page will document this limitation and note that unsupervised RF proximity (Shi & Horvath 2006) is a potential alternative for future work.

### Flow

1. `.distance_risk_prepare(X, Y, holdout, holdout_fraction, vars, na.rm, seed)`
2. `.rf_proximity(train, synthetic, vars, n_trees, mtry, importance = TRUE, seed)`
3. Per synthetic record: mean proximity to all training records (from terminal node matrix)
4. Push holdout through forest: `predict(forest, holdout, type = "terminalNodes")`
5. Per synthetic record: mean proximity to all holdout records
6. Compute `prox_share`, `prox_ratio`, Wilcoxon test
7. If `null_test`: permute train/holdout `n_null` times → null distribution
8. `privacy_pass = prox_share <= 0.55 & wilcox_p > 0.05`

### Metric definitions

Proximity is a **similarity** measure (higher = more similar), opposite to distance. Naming reflects this:

- **`prox_share`**: `mean(prox_train > prox_holdout)` — fraction of synthetic records with higher mean proximity to training than to holdout. ~0.5 = no memorization, >0.5 = memorization signal.
- **`prox_ratio`**: `mean(prox_train) / mean(prox_holdout)` — ratio of average proximities. ~1 = no memorization, >1 = memorization signal.
- **`prox_train`**: per-record vector of mean proximity to training records
- **`prox_holdout`**: per-record vector of mean proximity to holdout records

### `privacy_pass` threshold

The 0.55 threshold is carried over from DCR for consistency. Proximity-based share and distance-based share may have different null distributions. The `null_test = TRUE` permutation provides a principled, data-adaptive alternative: if the observed `prox_share` falls outside the 95th percentile of the null distribution, the test fails regardless of the 0.55 threshold. The help page will recommend interpreting the null test result when available, using 0.55 as a quick heuristic.

### Result fields

- `prox_share`: fraction with higher proximity to training than holdout
- `prox_ratio`: mean proximity ratio train/holdout
- `prox_train`: per-record mean proximity to training
- `prox_holdout`: per-record mean proximity to holdout
- `privacy_pass`: logical
- `wilcox_test`: Wilcoxon test object
- `null_distribution`: null stats (if `null_test = TRUE`)
- `variable_importance`: from the RF
- `n_synthetic`, `n_train`, `n_holdout`, `vars`

### S3 methods

- `print.rf_privacy`: one-line summary (prox_share, prox_ratio, pass/fail)
- `summary.rf_privacy`: detailed statistics, null distribution
- `print.summary.rf_privacy`: formatted summary output
- `plot.rf_privacy` with `which`:
  1. Paired density (proximity-to-training vs proximity-to-holdout)
  2. Per-record proximity difference histogram
  3. Null distribution with observed value

---

## Component 4: `propscore(method = "ranger")`

**File**: `R/propscore.R` (extend existing)

### Method naming

The existing `method = "rf"` uses `randomForest::randomForest()`. The new method uses `ranger::ranger()`. Both are random forests, but "ranger" is chosen because "rf" is taken. The distinction is:

- `method = "rf"`: standard RF propensity scores via randomForest (existing, backward-compatible)
- `method = "ranger"`: RF propensity scores via ranger + optional proximity-based structural utility analysis (new)
- `method = "logreg"`: logistic regression propensity scores (existing)

The help page will document when to prefer `"ranger"` over `"rf"`: when proximity-based structural analysis is desired, or for large datasets where ranger is faster.

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

### New branch

(alongside existing `"logreg"` and `"rf"` branches):
1. Call `.rf_proximity()` with the combined labeled data
2. Propensity scores from `forest$predictions[, 2]` (probability mode, OOB predictions)
3. pMSE via existing formula: `1/(2*n) * sum((p - cr)^2)`
4. If `proximity == "summary"`: compute structural metrics on-the-fly from terminal node matrix
5. If `proximity == "full"`: additionally store full n x n matrix
6. If `proximity == "none"`: skip proximity computation

### Additional result fields

- `oob_error`: OOB classification error (~0.5 = good utility)
- `variable_importance`: named numeric vector (when `importance = TRUE`)
- When `proximity != "none"`:
  - `within_orig_prox`: mean proximity among original records
  - `within_synth_prox`: mean proximity among synthetic records
  - `cross_prox`: mean proximity between original and synthetic
  - `structure_ratio`: `cross_prox / ((within_orig_prox + within_synth_prox) / 2)` — closer to 1 = better utility. Uses arithmetic mean of the two within-class proximities as denominator.
- When `proximity == "full"`:
  - `proximity_matrix`: full n x n matrix

### Scaling

`message()` when n > 10,000 with `proximity = "full"`, reporting estimated memory in MB and suggesting `proximity = "summary"`.

### Plot

New `which` type — boxplot of within-original, within-synthetic, cross-class proximity distributions. Only available when `proximity != "none"`. Also: new `which` type for variable importance bar chart (when `importance = TRUE`).

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

`n_trees` and `rf_global` silently ignored for other methods.

### Parameter handling for existing parameters

Parameters meaningful for other methods but not for RF:

- **`key`**: maps to `vars` in `.rf_proximity()`. The RF trains on `key` variables as features.
- **`type`** (ordinal detection): ignored for `method = "rf"`. RF handles ordinal/nominal natively.
- **`weights`**: ignored for `method = "rf"`. RF determines variable importance internally. A `message()` is emitted if `weights` is explicitly provided with `method = "rf"`, noting that RF uses data-driven importance instead.
- **`strategy`**: for RF, `"nearest"` means highest proximity. Other strategies (`"weighted"`, `"probabilistic"`) are not applicable — `message()` if explicitly set to non-default with `method = "rf"`.
- **`risk_weighting`**, **`kappa`**, **`bandwidth`**, **`kernel`**: ignored for `method = "rf"`.
- **`direction`**: supported. Controls which dataset is query vs search. `"forward"` (default): synthetic records are matched to original. `"reverse"`: original records matched to synthetic. The RF always trains on `rbind(original, synthetic)` regardless of direction; direction only determines which records are queries and which are candidates in the proximity-based matching step.
- **`truth`**: both `"row"` and `"id"` modes supported, same as other methods.

### Variable importance

Returned by default in `result$variable_importance`. Shows which QIs drive linkage risk.

### Matching logic

**Independent matching**: each query record matched to candidate with highest proximity.

**Bijective matching**: cost matrix = `1 - cross_proximity` → `clue::solve_LSAP()`. Binary risk (1 if matched to true record, 0 otherwise).

### Blocking modes

**`rf_global = FALSE` (default)**: train separate RF per block.
- Each forest learns within-block patterns only
- More precise, avoids cross-block contamination
- `message()` if any block has < 50 records, suggesting `rf_global = TRUE` (small blocks yield unreliable forests — ranger with `probability = TRUE` needs sufficient observations per class)
- Per-block variable importance aggregated (weighted by block size) into `result$variable_importance`

**`rf_global = TRUE`**: train one global RF.
- Extract terminal nodes for all records
- Restrict proximity to within-block pairs
- Robust to small blocks
- Single variable importance from the global forest

Help page documents the trade-off: per-block is more precise but needs sufficient block size (~50+ records per block); global is robust to small blocks but proximity reflects cross-block patterns.

**No blocking**: `message()` when n > 10,000, suggesting blocking.

### Helper function

The RF branch logic (especially per-block training) is factored into `.rf_linkage_block()` to keep the main `recordLinkage()` function readable (it is already ~2230 lines).

### Result

Same structure as other methods:
- `per_record`: risk, cand_n, true_in_set, bijective_assigned (for bijective)
- `variable_importance`: named numeric (new for rf)
- Score interpretation: proximity 0–1, higher = riskier match

**New plot `which` type**: variable importance bar chart.

---

## Scaling Strategy

Informational `message()` for large datasets, no hard stops:
- n > 10,000 with `proximity = "full"` in propscore → suggest `"summary"`, report estimated memory
- n > 10,000 without blocking in recordLinkage → suggest blocking
- Uses `message()` (suppressible with `suppressMessages()`), not `warning()`
- No approximate algorithms in v1

---

## Dependencies

- `ranger`: already in Suggests (used by RAPID). No new hard dependencies.
- `clue`: already in Suggests (used by bijective matching). Needed for `recordLinkage(method = "rf", matching = "bijective")`.

---

## Test Plan

| File | Tests |
|------|-------|
| `test-rf-proximity-internal.R` | `.rf_proximity()`: correct dimensions, proximity in [0,1], symmetric, identical data → high proximity. `.proximity_from_nodes()`: submatrix extraction. `.proximity_from_nodes_newdata()`: holdout push-through. |
| `test-distance-risk-prepare.R` | `.distance_risk_prepare()`: holdout splitting, variable intersection, NA removal. Regression: `dcr()` and `nndr()` existing tests still pass after refactor. |
| `test-rf-privacy.R` | Memorized data → high prox_share, random data → ~0.5, privacy_pass logic, null_test, S3 methods (print/summary/plot). |
| `test-propscore-ranger.R` | Identical distributions → OOB ~0.5, different distributions → OOB < 0.5, `proximity = "summary"/"full"/"none"`, structure_ratio, variable importance, message for large n, `match.arg` validation. New plot types (proximity boxplot, importance bar chart). |
| `test-recordLinkage-rf.R` | Independent + bijective matching, with/without blocking, rf_global TRUE/FALSE, variable importance, small block message, large n message, direction forward/reverse, truth = "row"/"id". New plot type (importance bar chart). |

---

## Integration

- `disclosure_report()`: add `rf_privacy` as optional metric
- `rumap()`: add `"rf_privacy"` to valid risk measures
- Both already discover metrics dynamically — minimal wiring

---

## Parking Lot (future work)

- Approximate proximity for very large n (subsampled proximity, sparse terminal node matching)
- Unsupervised RF proximity (Shi & Horvath 2006) as alternative to supervised real-vs-synthetic forest — addresses limitation of conflating distributional similarity with memorization
- Comparison vignette: RF proximity vs Gower vs Euclidean for record linkage
- Literature references: Breiman (2001) proximity; Lin & Jeon (2006) RF as adaptive kernel; Shi & Horvath (2006) unsupervised RF
- ranger hyperparameter tuning beyond n_trees/mtry (expose min.node.size?)
- Geometric mean alternative for `structure_ratio` denominator when within-class proximities differ substantially

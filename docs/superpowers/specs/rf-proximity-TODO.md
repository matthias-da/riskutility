# RF Proximity Approach — Brainstorming TODO

## 1. Unify distance-risk family into a common framework?

**Question**: Should `rf_privacy()`, `dcr()`, `nndr()`, `ims()` be layered into a single function with four methods?

**Current state**: `dcr()`, `nndr()`, `ims()` are standalone functions with nearly identical signatures:
- Same parameters: `X, Y, holdout, holdout_fraction, vars, method, na.rm, seed, ...`
- Same holdout-splitting boilerplate (~30 lines duplicated)
- Same distance computation (Gower/Euclidean), then diverge in what they compute from distances
- All return S3 objects with `privacy_pass`, `print/summary/plot` methods

**Decision**: **(B) Shared internal engine** — keep separate public functions, extract shared logic.

**Rationale**:
- No breaking API change for `dcr()`, `nndr()`, `ims()`
- IMS needs no holdout — awkward fit in a unified API
- RF proximity uses a different distance engine (terminal nodes, not Gower/Euclidean) — a unified function would have two separate code paths behind one `method` parameter
- The actual duplication is ~30 lines of holdout splitting + variable intersection + NA removal
- Better discoverability: each function gets its own help page

**Internal helper**:
```r
.distance_risk_prepare(X, Y, holdout, holdout_fraction,
                        vars, na.rm, seed)
# Returns list: $train, $synthetic, $holdout, $vars
```

- Called by `dcr()`, `nndr()`, `rf_privacy()` (all need holdout splitting)
- `ims()` calls a simpler variant (no holdout)
- Lives in `R/utils_internal.R`
- Refactor existing dcr/nndr to use it, then rf_privacy() gets it for free

---

## 2. RF proximity as 5th method in `recordLinkage()`

**Decision**: YES — add `method = "rf"` to `recordLinkage()`.

**Design sketch**:
- Train RF (ranger) on `rbind(original, synthetic)` with label column
- Extract cross-proximity matrix (synthetic x original block)
- Highest proximity = best match (independent matching)
- Cost = `1 - proximity` for bijective matching via `clue::solve_LSAP()`
- Supports `matching = c("independent", "bijective")` like other methods
- `key_vars` selects which variables the RF trains on
- `weights` parameter: not applicable (RF handles variable importance internally)

**Variable importance**: YES — returned by default in `result$variable_importance`. Shows which QIs drive linkage risk (data-driven, unlike user-specified weights in other methods). New `plot()` `which` type: variable importance bar chart.

**Blocking**:
- **Default**: train separate RF per block **(b)**. Conceptually cleaner — each forest learns within-block patterns only. Avoids cross-block contamination.
- **Optional**: `rf_global = TRUE` trains one global RF, restricts proximity to within-block pairs **(a)**. Useful when blocks are small (not enough data per block for a good forest).
- Help page documents the trade-off: per-block is more precise but needs sufficient block size; global is robust to small blocks but proximity reflects cross-block patterns.

**Scaling**: `message()` when n > 10,000 without blocking, suggesting blocking to reduce computation.

---

## 3. RF proximity as memorization/privacy metric — `rf_privacy()`

**Decision**: Standalone function in the distance-risk family.

**Design sketch**:
- Train RF on `rbind(training, synthetic)` with label
- For each synthetic record, compute proximity to all training records and all holdout records
- Compare: are synthetic records systematically closer (higher proximity) to training than holdout?
- Mirrors DCR logic but uses RF proximity instead of Gower/Euclidean distance
- Same parameters: `X, Y, holdout, holdout_fraction, vars, na.rm, seed, ...`
- Holdout records predicted through trained forest via `predict(..., type = "terminalNodes")` — not used in training
- Returns S3 object with `privacy_pass`, `print/summary/plot`

**Advantages over DCR**:
- Handles mixed data natively (no need for Gower)
- Captures non-linear relationships and variable interactions
- Data-adaptive distance (not predetermined metric)
- Potentially more sensitive to subtle memorization

**Output metrics** (parallel to DCR):
- `prox_share`: fraction of synthetic with higher mean proximity to training than holdout (~0.5 = good, >0.5 = memorization)
- `prox_ratio`: mean proximity to training / mean proximity to holdout (~1 = good, >1 = memorization)
- `prox_train` / `prox_holdout`: per-record mean proximity vectors

**Statistical tests**:
- Wilcoxon signed-rank test on paired proximity differences
- Permutation null distribution test (`null_test = TRUE`)

**`privacy_pass`**: `prox_share <= 0.55 AND wilcoxon p > 0.05` (mirrors DCR)

---

## 4. RF proximity as utility measure — extend `propscore()`

**Decision**: YES — add as new method inside `propscore()`.

**Current `propscore()` methods**:
- `method = "rf"` — uses `randomForest::randomForest()`, extracts `$votes[,2]` as propensity scores
- `method = "logreg"` — uses `glm(..., family = binomial())`, extracts `predict(..., type = "response")`

**Design**: Add `method = "ranger"` (since `"rf"` is taken by randomForest):
- Uses `ranger::ranger()` with `probability = TRUE` for propensity scores
- Additionally extracts RF proximity matrix via terminal node co-occurrence
- pMSE computed from predicted probabilities (same formula as existing methods)
- Proximity matrix stored in result for structural utility analysis

**What the proximity adds beyond classification accuracy**:
- Within-class proximity structure: do original-original pairs and synthetic-synthetic pairs have similar proximity distributions?
- Cross-class proximity: are original-synthetic pairs as proximate as within-class pairs? (if yes → high utility)
- This goes beyond "can you tell them apart?" to "do they have the same internal structure?"

**Additional result fields** (beyond existing propscore fields):
- `proximity_matrix`: full proximity matrix (optional, controlled by parameter)
- `within_orig_prox`: mean proximity among original records
- `within_synth_prox`: mean proximity among synthetic records
- `cross_prox`: mean proximity between original and synthetic
- `structure_ratio`: `cross_prox / mean(within_orig_prox, within_synth_prox)` — closer to 1 = better utility
- `oob_error`: OOB classification error (simpler utility summary; ~0.5 = good)
- `variable_importance`: which variables are most distinguishable

**`proximity` parameter**: `proximity = c("summary", "full", "none")`
- `"summary"` (default): compute structural metrics on-the-fly from terminal node matrix without storing full n×n matrix. Returns `structure_ratio`, `within_orig_prox`, `within_synth_prox`, `cross_prox`, `oob_error`.
- `"full"`: additionally store the full proximity matrix in result (for inspection, heatmap plots, custom analysis). O(n^2) memory.
- `"none"`: ranger used purely for propensity scores, no proximity overhead.

**Variable importance**: `importance = TRUE` by default (cheap, informative). Stored in `result$variable_importance`.

**Scaling**: `message()` when n > 10,000 with `proximity = "full"`, reporting estimated memory and suggesting `"summary"`.

**Plot**: add new `which` type for proximity structure comparison (within-class vs cross-class proximity distributions).

---

## 5. Shared `.rf_proximity()` internal engine

**Decision**: YES — single internal helper used by all three use cases.

**Interface sketch**:
```r
.rf_proximity(data1, data2,
              vars = NULL,
              n_trees = 500,
              mtry = NULL,
              seed = NULL, ...)
# Returns list:
#   $forest   — trained ranger object
#   $proximity — full proximity matrix (n1+n2) x (n1+n2)
#   $cross     — cross-block submatrix n2 x n1
#   $idx1, $idx2 — row indices for data1/data2 in combined data
```

**Implementation notes**:
- Uses `ranger::ranger()` with `write.forest = TRUE`
- Proximity via terminal node co-occurrence: predict all data, count shared leaves / n_trees
- For holdout (use case b): separate `predict(forest, holdout, type = "terminalNodes")` call
- `ranger` already in Suggests (used by RAPID)
- Memory: full proximity matrix is O(n^2). For n > ~10k, need sparse/chunked approach

---

## 6. Computational scaling strategy

**Decision**: Informational `message()` for large datasets, no hard stops.

- **n > 10,000 with `proximity = "full"`** in propscore: `message()` reporting estimated memory, suggesting `proximity = "summary"`
- **n > 10,000 without blocking** in recordLinkage: `message()` suggesting blocking
- No `force` parameter, no hard errors — user is informed, user decides
- Uses `message()` (not `warning()`) — informational, suppressible with `suppressMessages()`
- No approximate/chunked algorithms in v1. Can add later if demand arises.

---

## Parking Lot

- Should `rf_privacy()` be integrated into `disclosure_report()` and `rumap()`?
- ranger hyperparameter tuning: expose `num.trees`, `mtry`, `min.node.size`?
- Comparison vignette: RF proximity vs Gower vs Euclidean for record linkage
- Literature: Breiman (2001) proximity; Lin & Jeon (2006) RF as adaptive kernel; Shi & Horvath (2006) unsupervised RF

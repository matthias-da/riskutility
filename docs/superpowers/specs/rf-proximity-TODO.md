# RF Proximity Approach — Brainstorming TODO

**Status**: All design decisions made. Spec finalized (rev 3). Ready for implementation planning.

## Decisions Made

| # | Question | Decision |
|---|----------|----------|
| 1 | Unify distance-risk family? | **(B)** Shared `.distance_risk_prepare()`, separate public functions |
| 2 | Record linkage RF | 5th method in `recordLinkage()`, var_importance by default, per-block RF default + global option |
| 3 | rf_privacy() metrics | Max + mean proximity, null-test-derived privacy_pass, Wilcoxon as heuristic |
| 4 | Utility measure | Extend `propscore()` with `method = "ranger"`, `proximity = c("summary", "full", "none")` |
| 5 | Shared engine | `.rf_proximity()` internal helper with `modifyList` pattern |
| 6 | Scaling | `message()` for info, `warning()` for memory-intensive ops, no hard stops |

## Review Findings Addressed (rev 3)

All items from the 6-agent review have been incorporated:

### Must-fix (addressed)
- **A**: Max proximity added as primary metric (nearest-neighbor analog of DCR)
- **B**: `var_importance` field name (matches existing `recordLinkage`)
- **C**: `privacy_pass` derived from null distribution when `null_test = TRUE`
- **D**: Wilcoxon documented as heuristic, not formal test
- **E**: "When to use" decision guide table added
- **F**: Optimization strategy (tree-by-tree), runtime estimates, progress wiring
- **G**: Parameter order matches `dcr.default()`
- **H**: `na.rm = FALSE` default for RF (ranger handles NAs natively)
- **I**: Per-caller `requireNamespace()` guards + deep safety net
- **J**: `prox_ratio` zero-guard (`< 1/n_trees → NA`)
- **K**: In-bag vs OOB terminal nodes documented (conservative bias)
- **L**: Blocks < 2 per class fall back to global RF
- **M**: Aggregated small-block message (single summary, not per-block)
- **N**: propscore plot `which = 3` (proximity boxplot), `which = 4` (importance)
- **O**: `min_holdout` parameter in `.distance_risk_prepare()`
- **P**: `modifyList` pattern for ranger call
- **Q**: S3 generic + `synth_pair` method in spec
- **R**: High-cardinality factors → `respect.unordered.factors = "partition"`
- **S**: Self-contained RF proximity explanation in help pages
- **T**: Vignette section prescribed

### Minor items (addressed)
- Mid-rank tie correction for `prox_share` and `max_prox_share`
- Wilcoxon `alternative = "greater"` (not `"less"`)
- Permutation p-value: `(sum+1)/(n_null+1)`
- Pre-compute full syn-vs-real proximity for efficient null test
- `n_trees >= 10` validation
- `holdout_fraction` in (0,1) validation
- `@keywords internal` for new file
- `.rf_label` column name collision check
- Seed separation (`seed` for holdout, `seed + 1` for forest)
- `prox_*` prefix rationale documented
- `strategy` actual values listed (nearest, threshold, topk, etc.)
- Ranger version reproducibility note
- `structure_ratio` range and interpretation documented
- OOB prediction NA handling for small datasets
- `null_test` logic stays in callers, not shared helper
- recordLinkage RF uses existing `which = 3` for importance plot
- `synth_pair` dispatch tested with new params
- `set.seed()` side effect in `.distance_risk_prepare()` documented
- Mock print output for rf_privacy
- Help page structure prescribed (@details, @section, @seealso, @references, @examples)
- `warning()` (not message) for `proximity = "full"` with large n

# Embedding-Based Record Linkage via Autoencoder

## Goal

Add `method = "embedding"` to `recordLinkage()` — an autoencoder-based distance method that captures nonlinear relationships between quasi-identifiers via learned latent representations.

## Attacker Model

The adversary has access to the population's data structure and trains an autoencoder to learn a compact, informative representation. Records are linked via Euclidean distance in the latent space. This captures nonlinear dependencies between quasi-identifiers that linear methods (Gower, Mahalanobis) miss.

The autoencoder trains only on the original (query) data, modelling the attacker's knowledge of population structure without access to the specific anonymization. This is a deliberate design choice: the attacker knows the data-generating distribution but not which specific perturbation was applied. Training on query-only data means embeddings are optimized for the original distribution; search-side records that closely resemble originals will project nearby in the latent space, while heavily perturbed records will project further away. The trade-off is potential domain shift when anonymization changes the marginal distributions substantially — in such cases, the Mahalanobis or deterministic methods may be more appropriate.

## Dependencies

- `torch` in Suggests (not Imports); `luz` is **not** used — a manual training loop is simpler and avoids the complexity of `luz`'s `set_hparams()` pipeline for our custom composite loss
- `torch` must be added to `DESCRIPTION` Suggests field
- Gated at runtime: `requireNamespace("torch", quietly = TRUE)`
- Tests use `skip_if_not_installed("torch")`

## Architecture

### Autoencoder Model

```
Input: [entity_embed(cat1), ..., entity_embed(catK), scaled_num1, ..., scaled_numP]
  -> Encoder: Linear(input_dim, hidden_dim) + ReLU + Linear(hidden_dim, latent_dim) + ReLU
  -> Bottleneck: latent_dim features (the embedding)
  -> Decoder: Linear(latent_dim, hidden_dim) + ReLU + Linear(hidden_dim, output_dim)
  -> Output activations: sigmoid for numeric outputs, raw logits for categorical outputs
Output: reconstruct numerics (MSE on sigmoid outputs) + categoricals (cross-entropy on logits)
```

Where `hidden_dim = max(latent_dim + 2, floor(input_dim * 2 / 3))`.

- **Entity embeddings** for categorical variables: `min(50, floor(n_levels / 2) + 1)` dimensions per variable (Guo & Berkhahn 2016). R torch uses 1-based indexing, so category indices run from 1 to `n_levels`, with index `n_levels + 1` reserved as the `<UNK>` token for unseen categories (using `padding_idx = n_levels + 1L` in `nn_embedding`). The `<UNK>` embedding is initialized to the mean of all learned embeddings after training.
- **Numerics/ordinals**: min-max scaled to [0,1]
- **Latent dimension**: `max(2, floor(input_dim / 3))`, user-overridable via `emb_latent_dim`
- **Loss**: MSE (numeric reconstruction) + cross-entropy (categorical reconstruction), weighted by proportion of each variable type
- **Edge cases**: All-numeric data uses pure MSE loss (no entity embeddings, no cross-entropy). All-categorical data uses pure cross-entropy loss (no numeric scaling, no MSE).

### Training

- Manual training loop (no `luz` dependency): forward pass, composite loss, backward, optimizer step
- Optimizer: Adam, lr = 0.001
- Default epochs: 50 (`emb_epochs`)
- Batch size: 32 (or nrow if smaller)
- Early stopping: 20% validation holdout, patience = 5
- CPU by default (sufficient for typical SDC datasets < 50k rows)
- Reproducibility: `torch::torch_manual_seed()` called with the user's R seed (via `sample.int(.Machine$integer.max, 1L)`) before training, ensuring reproducible results when `set.seed()` is used

### Distance Computation

1. Train autoencoder on query data only
2. Encode both query and search data through the trained encoder
3. Compute pairwise Euclidean distances in latent space
4. Normalize to [0,1]: use the 97.5th percentile of within-query pairwise distances as threshold, with a fallback to `max(within-query distances)` when the block has fewer than 40 records (small-sample robustness)
5. Clamp: `pmin(d_normalized, 1)`

This normalization parallels the Mahalanobis chi-squared threshold approach.

### Variable Importance

Permutation-based: for each variable, shuffle it in the query data, re-encode, measure mean embedding shift compared to original embeddings. Normalized to sum to 1.

Computational cost is O(p * n * latent_dim) where p = number of variables and n = number of query records. For datasets with many variables (p > 30), this can be slow. Variable importance is computed by default but can be skipped by the orchestrator for large-p blocks in future optimizations.

## New Parameters

All ignored unless `method = "embedding"`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `emb_latent_dim` | integer or NULL | NULL (auto) | Latent dimension. Auto: `max(2, floor(input_dim / 3))` |
| `emb_epochs` | integer | 50L | Training epochs |
| `emb_global` | logical | FALSE | If TRUE, train one autoencoder globally; if FALSE, per block |

## Integration into recordLinkage()

### Pattern

Follows the deterministic method pattern (distance-based scores fed to strategy/risk_weighting):

- **Validation block** (pre-loop): Check `requireNamespace("torch")`, warn about ignored params
- **Per-block processing**: `.embedding_linkage_block()` orchestrates train/encode/distance
- **Score cache**: `score_cache[[i]] <- list(cand = indices, scores = distances, maximize = FALSE)` — distances are minimized (lower = closer match = higher risk), same direction as deterministic and Mahalanobis methods
- **Post-loop**: All matching modes (independent, bijective, OT) work automatically via score_cache. All strategies (`minmax`, `topk`, `rank`, `softmax`) and `risk_weighting` modes work as with deterministic (distance-based).
- **Direction**: `direction = "reverse"` is handled at the `recordLinkage()` level (swapping query/search before the method runs), not inside the embedding internals — same as all other methods.

### Per-Block vs Global

- `emb_global = FALSE` (default): Train separate autoencoder per block
- `emb_global = TRUE`: Train once on all query data, restrict distance lookups within blocks
- Small-block fallback: blocks with fewer than `max(30, 5 * latent_dim)` rows fall back to deterministic method (need sufficient data for meaningful autoencoder training)

### Unseen Categories

Categorical levels in the search data not seen during training get the `<UNK>` embedding (index `n_levels + 1`), initialized after training to the mean of all learned level embeddings for that variable. This uses R torch's 1-based indexing (no index 0).

## Internal Functions

File: `R/embedding_internal.R` (~300 lines)

| Function | Purpose |
|----------|---------|
| `.ae_model()` | Returns `torch::nn_module` defining the autoencoder (entity embeddings + 2-layer encoder + 2-layer decoder) |
| `.ae_preprocess()` | Min-max scale numerics, build category-to-index mappings (1-based). Returns list with tensors, mappings, metadata |
| `.ae_train()` | Manual training loop with Adam optimizer, early stopping, composite loss. Returns fitted model + preprocessing metadata |
| `.ae_encode()` | Pass data through encoder, return embedding matrix (n x latent_dim) |
| `.ae_distance()` | Euclidean distances from one embedding to candidates, normalized by threshold |
| `.ae_var_importance()` | Permutation importance: shuffle each var, re-encode, measure embedding shift |
| `.embedding_linkage_block()` | Orchestrator: preprocess -> train -> encode -> distance matrix. Returns `list(distances, var_importance, threshold)` |

## Result Object

- `risk`: Normalized distance-based risk in [0,1] (lower distance = higher risk)
- `d_true`, `d_min`, `d_rank`: Distances in latent space (normalized)
- `var_importance`: Permutation-based importance
- `settings$emb_latent_dim`, `settings$emb_epochs`, `settings$emb_global`: Stored for reproducibility
- All matching modes (independent, bijective, OT) work via score_cache
- All strategies and risk_weighting modes work (distance-based, like deterministic)

## Print/Summary Display

```
Method:      embedding
Matching:    independent
Embedding:   autoencoder (latent_dim = 5, epochs = 50)
```

## Roxygen Documentation

- `@section Embedding method:` describing the autoencoder approach, attacker model, and when to use it
- `@param emb_latent_dim`, `@param emb_epochs`, `@param emb_global`
- Update `@param method` to include `"embedding"`
- Reference: Guo, C. & Berkhahn, F. (2016). Entity Embeddings of Categorical Variables. arXiv:1604.06737.

## Test Plan

File: `tests/testthat/test-recordLinkage-embedding.R`

- Basic functionality with numeric data
- Basic functionality with all-categorical data
- Mixed data (numeric + categorical)
- Near-copy detection (truth = "row")
- Bijective matching
- OT matching
- Blocking (per-block training)
- Global training with blocking (`emb_global = TRUE`, `block` set)
- Direction = "reverse"
- Custom `emb_latent_dim` and `emb_epochs`
- Small-block fallback to deterministic
- Unseen categorical levels in search data
- Print and summary output
- Variable importance
- Reproducibility (same seed = same result)
- All tests gated with `skip_if_not_installed("torch")`

# Embedding-Based Record Linkage via Autoencoder

## Goal

Add `method = "embedding"` to `recordLinkage()` — an autoencoder-based distance method that captures nonlinear relationships between quasi-identifiers via learned latent representations.

## Attacker Model

The adversary has access to the population's data structure and trains an autoencoder to learn a compact, informative representation. Records are linked via Euclidean distance in the latent space. This captures nonlinear dependencies between quasi-identifiers that linear methods (Gower, Mahalanobis) miss.

The autoencoder trains only on the original (query) data, modelling the attacker's knowledge of population structure without access to the specific anonymization.

## Dependencies

- `torch` and `luz` in Suggests (not Imports)
- Gated at runtime: `requireNamespace("torch", quietly = TRUE)`
- Tests use `skip_if_not_installed("torch")`

## Architecture

### Autoencoder Model

```
Input: [entity_embed(cat1), ..., entity_embed(catK), scaled_num1, ..., scaled_numP]
  -> Encoder: Linear(input_dim, latent_dim) + ReLU
  -> Bottleneck: latent_dim features (the embedding)
  -> Decoder: Linear(latent_dim, input_dim) + output activations
Output: reconstruct numerics (MSE) + categoricals (cross-entropy)
```

- **Entity embeddings** for categorical variables: `min(50, n_levels / 2)` dimensions per variable (Guo & Berkhahn 2016)
- **Numerics/ordinals**: min-max scaled to [0,1]
- **Latent dimension**: `max(2, floor(input_dim / 3))`, user-overridable via `emb_latent_dim`
- **Loss**: MSE (numeric reconstruction) + cross-entropy (categorical reconstruction), weighted by proportion of each variable type

### Training

- Optimizer: Adam, lr = 0.001
- Default epochs: 50 (`emb_epochs`)
- Batch size: 32 (or nrow if smaller)
- Early stopping: 20% validation holdout, patience = 5
- CPU by default (sufficient for typical SDC datasets < 50k rows)

### Distance Computation

1. Train autoencoder on query data only
2. Encode both query and search data through the trained encoder
3. Compute pairwise Euclidean distances in latent space
4. Normalize to [0,1] using the 97.5th percentile of within-query pairwise distances as threshold
5. Clamp: `pmin(d_normalized, 1)`

This normalization parallels the Mahalanobis chi-squared threshold approach.

### Variable Importance

Permutation-based: for each variable, shuffle it in the query data, re-encode, measure mean embedding shift compared to original embeddings. Normalized to sum to 1.

## New Parameters

All ignored unless `method = "embedding"`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `emb_latent_dim` | integer or NULL | NULL (auto) | Latent dimension. Auto: `max(2, floor(input_dim / 3))` |
| `emb_epochs` | integer | 50L | Training epochs |
| `emb_global` | logical | FALSE | If TRUE, train one autoencoder globally; if FALSE, per block |

## Integration into recordLinkage()

### Pattern

Follows the RF method pattern exactly:

- **Validation block** (pre-loop): Check `requireNamespace("torch")`, warn about ignored params
- **Per-block processing**: `.embedding_linkage_block()` orchestrates train/encode/distance
- **Score cache**: `score_cache[[i]] <- list(cand = indices, scores = distances, maximize = FALSE)`
- **Post-loop**: bijective/OT matching work automatically via score_cache

### Per-Block vs Global

- `emb_global = FALSE` (default): Train separate autoencoder per block
- `emb_global = TRUE`: Train once on all query data, restrict distance lookups within blocks
- Small-block fallback: blocks with fewer than `2 * latent_dim` rows fall back to deterministic method

### Unseen Categories

Categorical levels in the search data not seen during training get a dedicated `<UNK>` embedding index (index 0), initialized to the mean of all learned embeddings.

## Internal Functions

File: `R/embedding_internal.R` (~250 lines)

| Function | Purpose |
|----------|---------|
| `.ae_model()` | Returns `torch::nn_module` defining the autoencoder (entity embeddings + encoder + decoder) |
| `.ae_preprocess()` | Min-max scale numerics, build category-to-index mappings. Returns list with tensors, mappings, metadata |
| `.ae_train()` | Train via `luz::fit()` with early stopping. Returns fitted model + preprocessing metadata |
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
- Mixed data (numeric + categorical)
- Near-copy detection (truth = "row")
- Bijective matching
- OT matching
- Blocking (per-block training)
- Global training (`emb_global = TRUE`)
- Direction = "reverse"
- Custom `emb_latent_dim` and `emb_epochs`
- Small-block fallback to deterministic
- Unseen categorical levels in search data
- Print and summary output
- Variable importance
- All tests gated with `skip_if_not_installed("torch")`

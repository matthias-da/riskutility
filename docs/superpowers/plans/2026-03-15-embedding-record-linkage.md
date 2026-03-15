# Embedding-Based Record Linkage Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `method = "embedding"` to `recordLinkage()` — an autoencoder-based distance method using torch entity embeddings for mixed data.

**Architecture:** A new internal file `R/embedding_internal.R` contains the autoencoder model, preprocessing, training, encoding, and distance functions. The main `recordLinkage()` function gets a new branch paralleling the RF method pattern but with `maximize = FALSE` (distance-based, like deterministic/Mahalanobis). Score cache integration means bijective/OT matching work automatically.

**Tech Stack:** R torch (nn_module, nn_embedding, nn_linear, optim_adam), no luz dependency. Manual training loop with early stopping.

**Spec:** `docs/superpowers/specs/2026-03-15-embedding-record-linkage-design.md`

---

## File Structure

| File | Action | Responsibility |
|------|--------|---------------|
| `R/embedding_internal.R` | Create | Autoencoder model, preprocessing, training, encoding, distance, variable importance |
| `R/recordLinkage.R` | Modify | Add "embedding" method: signature, validation, per-block processing, var importance, settings, print/summary |
| `DESCRIPTION` | Modify | Add `torch` to Suggests |
| `tests/testthat/test-recordLinkage-embedding.R` | Create | All embedding method tests |

---

## Chunk 1: Internal Autoencoder Engine

### Task 1: Preprocessing — `.ae_preprocess()`

**Files:**
- Create: `R/embedding_internal.R`
- Test: `tests/testthat/test-recordLinkage-embedding.R`

This function takes a data frame and key variables, identifies numeric vs categorical columns, builds category-to-index mappings (1-based for torch), min-max scales numerics to [0,1], and returns everything needed to construct tensors.

- [ ] **Step 1: Write `.ae_preprocess()` implementation**

Create `R/embedding_internal.R` with:

```r
#' Preprocess data for autoencoder embedding
#'
#' Min-max scales numeric columns to [0,1], builds 1-based category-to-index
#' mappings for categorical columns, and computes metadata for model construction.
#'
#' @param data data.frame to preprocess
#' @param key character vector of column names to use
#' @param type named character vector of variable types (from recordLinkage)
#' @param cat_maps optional named list of existing category mappings (for
#'   encoding new data with a trained model's vocabulary)
#' @param num_ranges optional named list of existing numeric ranges (for
#'   encoding new data with training-time scaling)
#' @return list with:
#'   \describe{
#'     \item{num_mat}{numeric matrix of scaled numerics (n x p_num), or NULL}
#'     \item{cat_idx}{list of integer vectors, one per categorical variable,
#'       each of length n with 1-based indices (unseen levels get n_levels + 1)}
#'     \item{cat_maps}{named list: variable -> named integer vector mapping
#'       level -> index (1-based)}
#'     \item{num_ranges}{named list: variable -> c(min, max)}
#'     \item{num_keys}{character vector of numeric key names}
#'     \item{cat_keys}{character vector of categorical key names}
#'     \item{n_levels}{named integer vector: n_levels per categorical var
#'       (training levels only, excluding UNK)}
#'     \item{emb_dims}{named integer vector: embedding dimension per categorical var}
#'     \item{input_dim}{integer, total input dimension after entity embeddings}
#'   }
#' @keywords internal
.ae_preprocess <- function(data, key, type, cat_maps = NULL,
                           num_ranges = NULL) {
  n <- nrow(data)

  # Classify variables
  num_keys <- character(0)
  cat_keys <- character(0)
  for (v in key) {
    if (!is.null(type) && v %in% names(type)) {
      if (type[v] %in% c("numeric", "ordinal")) {
        num_keys <- c(num_keys, v)
      } else {
        cat_keys <- c(cat_keys, v)
      }
    } else if (is.numeric(data[[v]]) || is.integer(data[[v]])) {
      num_keys <- c(num_keys, v)
    } else {
      cat_keys <- c(cat_keys, v)
    }
  }

  # Numeric: min-max scale to [0,1]
  num_mat <- NULL
  if (length(num_keys) > 0) {
    if (is.null(num_ranges)) {
      num_ranges <- setNames(lapply(num_keys, function(v) {
        vals <- as.numeric(data[[v]])
        c(min = min(vals, na.rm = TRUE), max = max(vals, na.rm = TRUE))
      }), num_keys)
    }
    num_mat <- matrix(0, nrow = n, ncol = length(num_keys))
    colnames(num_mat) <- num_keys
    for (j in seq_along(num_keys)) {
      v <- num_keys[j]
      rng <- num_ranges[[v]]
      vals <- as.numeric(data[[v]])
      span <- rng[2] - rng[1]
      if (span > 0) {
        num_mat[, j] <- (vals - rng[1]) / span
      } else {
        num_mat[, j] <- 0.5
      }
      # Clamp to [0,1] for out-of-range values in new data
      num_mat[, j] <- pmin(pmax(num_mat[, j], 0), 1)
    }
  } else {
    num_ranges <- list()
  }

  # Categorical: 1-based index mapping
  if (is.null(cat_maps)) {
    cat_maps <- setNames(lapply(cat_keys, function(v) {
      lvls <- sort(unique(as.character(data[[v]])))
      setNames(seq_along(lvls), lvls)
    }), cat_keys)
  }

  n_levels <- setNames(vapply(cat_keys, function(v) {
    length(cat_maps[[v]])
  }, integer(1)), cat_keys)

  cat_idx <- setNames(lapply(cat_keys, function(v) {
    vals <- as.character(data[[v]])
    mapped <- cat_maps[[v]][vals]
    # Unseen levels -> n_levels + 1 (UNK index)
    mapped[is.na(mapped)] <- n_levels[v] + 1L
    as.integer(mapped)
  }), cat_keys)

  # Embedding dimensions: min(50, floor(n_levels / 2) + 1)
  emb_dims <- setNames(vapply(cat_keys, function(v) {
    as.integer(min(50L, floor(n_levels[v] / 2) + 1L))
  }, integer(1)), cat_keys)

  # Total input dimension = sum(emb_dims) + n_numeric
  input_dim <- sum(emb_dims) + length(num_keys)

  list(
    num_mat = num_mat,
    cat_idx = cat_idx,
    cat_maps = cat_maps,
    num_ranges = num_ranges,
    num_keys = num_keys,
    cat_keys = cat_keys,
    n_levels = n_levels,
    emb_dims = emb_dims,
    input_dim = as.integer(input_dim)
  )
}
```

- [ ] **Step 2: Write basic preprocessing test**

Add to `tests/testthat/test-recordLinkage-embedding.R`:

```r
# tests/testthat/test-recordLinkage-embedding.R

test_that(".ae_preprocess handles numeric data", {
  skip_if_not_installed("torch")
  df <- data.frame(a = c(0, 5, 10), b = c(100, 200, 300))
  res <- .ae_preprocess(df, key = c("a", "b"), type = NULL)

  expect_equal(res$num_keys, c("a", "b"))
  expect_equal(length(res$cat_keys), 0)
  expect_equal(nrow(res$num_mat), 3)
  expect_equal(ncol(res$num_mat), 2)
  # Min-max scaled: a -> 0, 0.5, 1; b -> 0, 0.5, 1
  expect_equal(res$num_mat[, "a"], c(0, 0.5, 1))
  expect_equal(res$num_mat[, "b"], c(0, 0.5, 1))
  expect_equal(res$input_dim, 2L)
})

test_that(".ae_preprocess handles categorical data", {
  skip_if_not_installed("torch")
  df <- data.frame(x = c("A", "B", "C", "A"), stringsAsFactors = FALSE)
  res <- .ae_preprocess(df, key = "x", type = NULL)

  expect_equal(res$cat_keys, "x")
  expect_equal(length(res$num_keys), 0)
  expect_null(res$num_mat)
  expect_equal(res$n_levels[["x"]], 3L)
  # 1-based indices: A=1, B=2, C=3
  expect_equal(res$cat_idx[["x"]], c(1L, 2L, 3L, 1L))
  # emb_dim = min(50, floor(3/2)+1) = 2
  expect_equal(res$emb_dims[["x"]], 2L)
})

test_that(".ae_preprocess handles mixed data", {
  skip_if_not_installed("torch")
  df <- data.frame(a = c(1, 2, 3), b = c("X", "Y", "X"),
                   stringsAsFactors = FALSE)
  res <- .ae_preprocess(df, key = c("a", "b"), type = NULL)

  expect_equal(res$num_keys, "a")
  expect_equal(res$cat_keys, "b")
  # input_dim = 1 (numeric) + emb_dim(b)
  expect_equal(res$input_dim, 1L + res$emb_dims[["b"]])
})

test_that(".ae_preprocess maps unseen levels to UNK", {
  skip_if_not_installed("torch")
  df_train <- data.frame(x = c("A", "B"), stringsAsFactors = FALSE)
  prep_train <- .ae_preprocess(df_train, key = "x", type = NULL)

  df_new <- data.frame(x = c("A", "C"), stringsAsFactors = FALSE)
  prep_new <- .ae_preprocess(df_new, key = "x", type = NULL,
                             cat_maps = prep_train$cat_maps)

  # A -> 1, C (unseen) -> n_levels + 1 = 3
  expect_equal(prep_new$cat_idx[["x"]], c(1L, 3L))
})
```

- [ ] **Step 3: Run tests to verify preprocessing**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-recordLinkage-embedding.R')"`
Expected: 4 tests pass (or skip if torch not installed)

- [ ] **Step 4: Commit**

```bash
git add R/embedding_internal.R tests/testthat/test-recordLinkage-embedding.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add .ae_preprocess() for embedding record linkage"
```

---

### Task 2: Autoencoder Model — `.ae_model()`

**Files:**
- Modify: `R/embedding_internal.R`
- Modify: `tests/testthat/test-recordLinkage-embedding.R`

Defines the `torch::nn_module` with entity embeddings for categoricals and a 2-layer encoder/decoder.

- [ ] **Step 1: Write `.ae_model()` implementation**

Append to `R/embedding_internal.R`:

```r
#' Define autoencoder nn_module for mixed data
#'
#' Creates a torch nn_module with entity embeddings for categorical variables
#' and a 2-layer encoder (input -> hidden -> latent) / decoder (latent -> hidden -> output).
#'
#' @param prep list from `.ae_preprocess()` with input_dim, cat_keys, n_levels,
#'   emb_dims, num_keys
#' @param latent_dim integer, bottleneck dimension
#' @return instantiated torch nn_module
#' @keywords internal
.ae_model <- function(prep, latent_dim) {
  n_num <- length(prep$num_keys)
  cat_keys <- prep$cat_keys
  n_levels <- prep$n_levels
  emb_dims <- prep$emb_dims
  input_dim <- prep$input_dim
  hidden_dim <- as.integer(max(latent_dim + 2L, floor(input_dim * 2 / 3)))

  # Output dimension: n_numeric + sum(n_levels) for categorical logits
  output_dim <- n_num + sum(n_levels)

  ae <- torch::nn_module(
    "RecordLinkageAE",

    initialize = function() {
      # Entity embeddings: one nn_embedding per categorical variable
      # num_embeddings = n_levels + 1 (extra for UNK), with padding_idx
      self$embeddings <- torch::nn_module_list()
      for (v in cat_keys) {
        emb <- torch::nn_embedding(
          num_embeddings = n_levels[v] + 1L,
          embedding_dim = emb_dims[v],
          padding_idx = n_levels[v] + 1L  # UNK index (1-based)
        )
        self$embeddings$append(emb)
      }

      # Encoder: 2 layers
      self$enc1 <- torch::nn_linear(input_dim, hidden_dim)
      self$enc2 <- torch::nn_linear(hidden_dim, latent_dim)

      # Decoder: 2 layers
      self$dec1 <- torch::nn_linear(latent_dim, hidden_dim)
      self$dec2 <- torch::nn_linear(hidden_dim, output_dim)

      # Store metadata for forward pass
      self$n_num <- n_num
      self$cat_keys <- cat_keys
      self$n_levels <- n_levels
    },

    encode = function(num_tensor, cat_tensors) {
      parts <- list()

      # Entity embeddings for each categorical variable
      for (k in seq_along(self$cat_keys)) {
        parts[[k]] <- self$embeddings[[k]](cat_tensors[[k]])
      }

      # Append numeric features
      if (self$n_num > 0 && !is.null(num_tensor)) {
        parts[[length(parts) + 1]] <- num_tensor
      }

      x <- torch::torch_cat(parts, dim = 2)
      x <- torch::nnf_relu(self$enc1(x))
      torch::nnf_relu(self$enc2(x))
    },

    decode = function(z) {
      x <- torch::nnf_relu(self$dec1(z))
      self$dec2(x)
    },

    forward = function(num_tensor, cat_tensors) {
      z <- self$encode(num_tensor, cat_tensors)
      out <- self$decode(z)

      # Split output into numeric (sigmoid) and categorical (logits) parts
      list(
        latent = z,
        num_recon = if (self$n_num > 0) {
          torch::torch_sigmoid(out[, 1:self$n_num, drop = FALSE])
        } else {
          NULL
        },
        cat_logits = if (length(self$cat_keys) > 0) {
          offset <- self$n_num
          logits <- list()
          for (k in seq_along(self$cat_keys)) {
            v <- self$cat_keys[k]
            nl <- self$n_levels[v]
            logits[[k]] <- out[, (offset + 1):(offset + nl), drop = FALSE]
            offset <- offset + nl
          }
          logits
        } else {
          NULL
        }
      )
    }
  )

  ae()
}
```

- [ ] **Step 2: Write model construction test**

Append to test file:

```r
test_that(".ae_model constructs valid autoencoder", {
  skip_if_not_installed("torch")
  df <- data.frame(a = c(1, 2, 3, 4, 5),
                   b = c("X", "Y", "X", "Y", "X"),
                   stringsAsFactors = FALSE)
  prep <- .ae_preprocess(df, key = c("a", "b"), type = NULL)
  model <- .ae_model(prep, latent_dim = 2L)

  expect_true(inherits(model, "nn_module"))

  # Forward pass should work
  num_t <- torch::torch_tensor(prep$num_mat, dtype = torch::torch_float())
  cat_t <- lapply(prep$cat_idx, function(idx) {
    torch::torch_tensor(idx, dtype = torch::torch_long())
  })
  out <- model(num_t, cat_t)
  expect_equal(length(out$latent$shape), 2)
  expect_equal(out$latent$shape[1], 5)  # n rows
  expect_equal(out$latent$shape[2], 2)  # latent_dim
})

test_that(".ae_model works with all-numeric data", {
  skip_if_not_installed("torch")
  df <- data.frame(a = 1:5, b = 6:10)
  prep <- .ae_preprocess(df, key = c("a", "b"), type = NULL)
  model <- .ae_model(prep, latent_dim = 2L)

  num_t <- torch::torch_tensor(prep$num_mat, dtype = torch::torch_float())
  out <- model(num_t, list())
  expect_equal(out$latent$shape[2], 2)
  expect_true(!is.null(out$num_recon))
  expect_null(out$cat_logits)
})

test_that(".ae_model works with all-categorical data", {
  skip_if_not_installed("torch")
  df <- data.frame(x = c("A", "B", "C", "A", "B"),
                   y = c("P", "Q", "P", "Q", "P"),
                   stringsAsFactors = FALSE)
  prep <- .ae_preprocess(df, key = c("x", "y"), type = NULL)
  model <- .ae_model(prep, latent_dim = 2L)

  cat_t <- lapply(prep$cat_idx, function(idx) {
    torch::torch_tensor(idx, dtype = torch::torch_long())
  })
  out <- model(NULL, cat_t)
  expect_equal(out$latent$shape[2], 2)
  expect_null(out$num_recon)
  expect_true(!is.null(out$cat_logits))
})
```

- [ ] **Step 3: Run tests**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-recordLinkage-embedding.R')"`
Expected: 7 tests pass

- [ ] **Step 4: Commit**

```bash
git add R/embedding_internal.R tests/testthat/test-recordLinkage-embedding.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add .ae_model() autoencoder nn_module for mixed data"
```

---

### Task 3: Training Loop — `.ae_train()`

**Files:**
- Modify: `R/embedding_internal.R`
- Modify: `tests/testthat/test-recordLinkage-embedding.R`

Manual training loop with Adam optimizer, composite MSE + cross-entropy loss, early stopping on validation holdout.

- [ ] **Step 1: Write `.ae_train()` implementation**

Append to `R/embedding_internal.R`:

```r
#' Train autoencoder on data
#'
#' Manual training loop with Adam optimizer, composite loss (MSE for numerics +
#' cross-entropy for categoricals), and early stopping on a 20% validation
#' holdout.
#'
#' @param data data.frame training data
#' @param key character vector of column names
#' @param type named character vector of variable types (or NULL)
#' @param latent_dim integer, bottleneck dimension
#' @param epochs integer, maximum training epochs
#' @param batch_size integer, mini-batch size
#' @param lr numeric, learning rate for Adam
#' @param patience integer, early stopping patience
#' @param val_fraction numeric, fraction of data for validation
#' @return list with:
#'   \describe{
#'     \item{model}{trained torch nn_module}
#'     \item{prep}{preprocessing metadata from .ae_preprocess()}
#'     \item{latent_dim}{integer, latent dimension used}
#'   }
#' @keywords internal
.ae_train <- function(data, key, type = NULL, latent_dim = NULL,
                      epochs = 50L, batch_size = 32L, lr = 0.001,
                      patience = 5L, val_fraction = 0.2) {

  prep <- .ae_preprocess(data, key, type)

  # Auto latent dimension
  if (is.null(latent_dim)) {
    latent_dim <- as.integer(max(2L, floor(prep$input_dim / 3)))
  }

  # Seed torch from R's RNG for reproducibility
  torch::torch_manual_seed(sample.int(.Machine$integer.max, 1L))

  model <- .ae_model(prep, latent_dim)

  optimizer <- torch::optim_adam(model$parameters, lr = lr)

  n <- nrow(data)
  batch_size <- min(batch_size, n)

  # Validation split
  n_val <- max(1L, as.integer(floor(n * val_fraction)))
  n_train <- n - n_val
  perm <- sample.int(n)
  train_idx <- perm[seq_len(n_train)]
  val_idx <- perm[(n_train + 1L):n]

  # Helper: build tensors for a subset of rows
  .make_tensors <- function(rows) {
    num_t <- if (!is.null(prep$num_mat)) {
      torch::torch_tensor(prep$num_mat[rows, , drop = FALSE],
                          dtype = torch::torch_float())
    } else {
      NULL
    }
    cat_t <- lapply(prep$cat_idx, function(idx) {
      torch::torch_tensor(idx[rows], dtype = torch::torch_long())
    })
    list(num = num_t, cat = cat_t)
  }

  # Helper: compute composite loss
  .compute_loss <- function(out, num_target, cat_targets) {
    losses <- list()
    n_terms <- 0L

    if (!is.null(out$num_recon) && !is.null(num_target)) {
      losses$mse <- torch::nnf_mse_loss(out$num_recon, num_target)
      n_terms <- n_terms + length(prep$num_keys)
    }

    if (!is.null(out$cat_logits) && length(out$cat_logits) > 0) {
      ce_sum <- torch::torch_tensor(0, dtype = torch::torch_float())
      for (k in seq_along(prep$cat_keys)) {
        # Target indices must be 1-based, matching training levels
        # Cross-entropy expects (batch, n_classes) logits vs (batch,) targets
        ce_sum <- ce_sum + torch::nnf_cross_entropy(
          out$cat_logits[[k]], cat_targets[[k]]
        )
      }
      losses$ce <- ce_sum / length(prep$cat_keys)
      n_terms <- n_terms + length(prep$cat_keys)
    }

    if (length(losses) == 0) return(torch::torch_tensor(0))

    # Weight by proportion of each type
    total <- torch::torch_tensor(0, dtype = torch::torch_float())
    n_vars <- length(prep$num_keys) + length(prep$cat_keys)
    if (!is.null(losses$mse)) {
      total <- total + losses$mse * (length(prep$num_keys) / n_vars)
    }
    if (!is.null(losses$ce)) {
      total <- total + losses$ce * (length(prep$cat_keys) / n_vars)
    }
    total
  }

  # Prepare validation tensors (fixed)
  val_tensors <- .make_tensors(val_idx)
  val_cat_targets <- lapply(prep$cat_idx, function(idx) {
    torch::torch_tensor(idx[val_idx], dtype = torch::torch_long())
  })

  best_val_loss <- Inf
  patience_counter <- 0L

  for (epoch in seq_len(epochs)) {
    model$train()

    # Shuffle training data
    train_perm <- sample(train_idx)
    n_batches <- ceiling(n_train / batch_size)

    for (b in seq_len(n_batches)) {
      start <- (b - 1L) * batch_size + 1L
      end <- min(b * batch_size, n_train)
      batch_rows <- train_perm[start:end]

      batch_t <- .make_tensors(batch_rows)
      batch_cat_targets <- lapply(prep$cat_idx, function(idx) {
        torch::torch_tensor(idx[batch_rows], dtype = torch::torch_long())
      })

      optimizer$zero_grad()
      out <- model(batch_t$num, batch_t$cat)
      loss <- .compute_loss(out, batch_t$num, batch_cat_targets)
      loss$backward()
      optimizer$step()
    }

    # Validation loss
    model$eval()
    torch::with_no_grad({
      val_out <- model(val_tensors$num, val_tensors$cat)
      val_loss <- .compute_loss(val_out, val_tensors$num, val_cat_targets)
      val_loss_val <- val_loss$item()
    })

    # Early stopping
    if (val_loss_val < best_val_loss) {
      best_val_loss <- val_loss_val
      patience_counter <- 0L
    } else {
      patience_counter <- patience_counter + 1L
      if (patience_counter >= patience) break
    }
  }

  # Initialize UNK embeddings to mean of learned embeddings
  torch::with_no_grad({
    for (k in seq_along(prep$cat_keys)) {
      v <- prep$cat_keys[k]
      nl <- prep$n_levels[v]
      # Learned embeddings are indices 1:nl; UNK is nl+1
      learned <- model$embeddings[[k]]$weight[1:nl, ]
      mean_emb <- learned$mean(dim = 1)
      model$embeddings[[k]]$weight[nl + 1L, ] <- mean_emb
    }
  })

  list(
    model = model,
    prep = prep,
    latent_dim = latent_dim
  )
}
```

- [ ] **Step 2: Write training test**

Append to test file:

```r
test_that(".ae_train produces a trained model", {
  skip_if_not_installed("torch")
  set.seed(42)
  df <- data.frame(a = rnorm(50), b = sample(letters[1:3], 50, TRUE),
                   stringsAsFactors = FALSE)
  res <- .ae_train(df, key = c("a", "b"), epochs = 5L)

  expect_true(inherits(res$model, "nn_module"))
  expect_true(is.list(res$prep))
  expect_true(is.integer(res$latent_dim))
  expect_true(res$latent_dim >= 2L)
})

test_that(".ae_train respects custom latent_dim", {
  skip_if_not_installed("torch")
  set.seed(42)
  df <- data.frame(a = rnorm(30), b = rnorm(30))
  res <- .ae_train(df, key = c("a", "b"), latent_dim = 4L, epochs = 3L)

  expect_equal(res$latent_dim, 4L)
})

test_that(".ae_train is reproducible with set.seed", {
  skip_if_not_installed("torch")
  df <- data.frame(a = rnorm(30), b = rnorm(30))

  set.seed(123)
  res1 <- .ae_train(df, key = c("a", "b"), latent_dim = 2L, epochs = 5L)

  set.seed(123)
  res2 <- .ae_train(df, key = c("a", "b"), latent_dim = 2L, epochs = 5L)

  # Same seed -> same model weights
  w1 <- as.numeric(res1$model$enc1$weight$cpu())
  w2 <- as.numeric(res2$model$enc1$weight$cpu())
  expect_equal(w1, w2, tolerance = 1e-6)
})
```

- [ ] **Step 3: Run tests**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-recordLinkage-embedding.R')"`
Expected: 10 tests pass

- [ ] **Step 4: Commit**

```bash
git add R/embedding_internal.R tests/testthat/test-recordLinkage-embedding.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add .ae_train() manual training loop with early stopping"
```

---

### Task 4: Encoding and Distance — `.ae_encode()`, `.ae_distance()`

**Files:**
- Modify: `R/embedding_internal.R`
- Modify: `tests/testthat/test-recordLinkage-embedding.R`

- [ ] **Step 1: Write `.ae_encode()` and `.ae_distance()` implementations**

Append to `R/embedding_internal.R`:

```r
#' Encode data through trained autoencoder
#'
#' Passes data through the encoder portion of the trained autoencoder,
#' returning an n x latent_dim embedding matrix.
#'
#' @param model trained torch nn_module from .ae_train()
#' @param data data.frame to encode
#' @param prep preprocessing metadata from .ae_preprocess() (training-time)
#' @param key character vector of column names
#' @param type named character vector of variable types (or NULL)
#' @return numeric matrix (n x latent_dim)
#' @keywords internal
.ae_encode <- function(model, data, prep, key, type = NULL) {
  # Preprocess with training-time mappings
  new_prep <- .ae_preprocess(data, key, type,
                             cat_maps = prep$cat_maps,
                             num_ranges = prep$num_ranges)

  num_t <- if (!is.null(new_prep$num_mat)) {
    torch::torch_tensor(new_prep$num_mat, dtype = torch::torch_float())
  } else {
    NULL
  }
  cat_t <- lapply(new_prep$cat_idx, function(idx) {
    torch::torch_tensor(idx, dtype = torch::torch_long())
  })

  model$eval()
  emb <- torch::with_no_grad({
    model$encode(num_t, cat_t)
  })

  as.matrix(emb$cpu())
}


#' Compute normalized Euclidean distances in latent space
#'
#' Computes pairwise Euclidean distances between query and search embeddings,
#' normalized to [0,1] using the 97.5th percentile of within-query distances
#' (with fallback to max for small datasets).
#'
#' @param emb_query numeric matrix (n_q x latent_dim)
#' @param emb_search numeric matrix (n_s x latent_dim)
#' @return list with:
#'   \describe{
#'     \item{dist_mat}{numeric matrix (n_q x n_s) of normalized distances}
#'     \item{threshold}{numeric, normalization threshold used}
#'   }
#' @keywords internal
.ae_distance <- function(emb_query, emb_search) {
  n_q <- nrow(emb_query)

  n_s <- nrow(emb_search)

  # Cross-distance: ||q_i - s_j||^2 = ||q_i||^2 + ||s_j||^2 - 2*q_i.s_j
  # Memory-efficient: only compute the n_q x n_s cross block
  q_sq <- rowSums(emb_query^2)
  s_sq <- rowSums(emb_search^2)
  cross_sq <- outer(q_sq, s_sq, "+") - 2 * tcrossprod(emb_query, emb_search)
  cross_sq[cross_sq < 0] <- 0  # numerical cleanup
  cross_dist <- sqrt(cross_sq)

  # Within-query distances for threshold
  if (n_q > 1) {
    q_sq_mat <- outer(q_sq, q_sq, "+") - 2 * tcrossprod(emb_query)
    q_sq_mat[q_sq_mat < 0] <- 0
    q_dist <- sqrt(q_sq_mat)
    q_upper <- q_dist[upper.tri(q_dist)]
    if (length(q_upper) > 0 && any(q_upper > 0)) {
      if (n_q < 40) {
        threshold <- max(q_upper)
      } else {
        threshold <- as.numeric(stats::quantile(q_upper, 0.975))
      }
    } else {
      threshold <- 1
    }
  } else {
    threshold <- if (max(cross_dist) > 0) max(cross_dist) else 1
  }

  if (threshold <= 0) threshold <- 1

  # Normalize and clamp
  cross_dist <- cross_dist / threshold
  cross_dist <- pmin(cross_dist, 1)

  list(dist_mat = cross_dist, threshold = threshold)
}
```

- [ ] **Step 2: Write encoding and distance tests**

Append to test file:

```r
test_that(".ae_encode returns correct dimensions", {
  skip_if_not_installed("torch")
  set.seed(42)
  df <- data.frame(a = rnorm(30), b = rnorm(30))
  trained <- .ae_train(df, key = c("a", "b"), latent_dim = 3L, epochs = 3L)

  emb <- .ae_encode(trained$model, df, trained$prep,
                    key = c("a", "b"))
  expect_equal(nrow(emb), 30)
  expect_equal(ncol(emb), 3)
  expect_true(all(is.finite(emb)))
})

test_that(".ae_encode handles unseen categorical levels", {
  skip_if_not_installed("torch")
  set.seed(42)
  df_train <- data.frame(a = 1:20,
                         b = rep(c("X", "Y"), 10),
                         stringsAsFactors = FALSE)
  trained <- .ae_train(df_train, key = c("a", "b"),
                       latent_dim = 2L, epochs = 3L)

  df_new <- data.frame(a = 1:5,
                       b = c("X", "Y", "Z", "X", "W"),  # Z, W unseen
                       stringsAsFactors = FALSE)
  emb <- .ae_encode(trained$model, df_new, trained$prep,
                    key = c("a", "b"))
  expect_equal(nrow(emb), 5)
  expect_true(all(is.finite(emb)))
})

test_that(".ae_distance returns normalized distances", {
  skip_if_not_installed("torch")
  set.seed(42)
  emb_q <- matrix(rnorm(20), 10, 2)
  emb_s <- matrix(rnorm(14), 7, 2)
  res <- .ae_distance(emb_q, emb_s)

  expect_equal(nrow(res$dist_mat), 10)
  expect_equal(ncol(res$dist_mat), 7)
  expect_true(all(res$dist_mat >= 0))
  expect_true(all(res$dist_mat <= 1))
  expect_true(res$threshold > 0)
})
```

- [ ] **Step 3: Run tests**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-recordLinkage-embedding.R')"`
Expected: 13 tests pass

- [ ] **Step 4: Commit**

```bash
git add R/embedding_internal.R tests/testthat/test-recordLinkage-embedding.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add .ae_encode() and .ae_distance() for latent space distances"
```

---

### Task 5: Variable Importance and Block Orchestrator

**Files:**
- Modify: `R/embedding_internal.R`
- Modify: `tests/testthat/test-recordLinkage-embedding.R`

- [ ] **Step 1: Write `.ae_var_importance()` and `.embedding_linkage_block()`**

Append to `R/embedding_internal.R`:

```r
#' Permutation-based variable importance in latent space
#'
#' For each variable, shuffles its values in the data, re-encodes, and measures
#' the mean Euclidean shift in the latent space compared to the original
#' embeddings. Normalized to sum to 1.
#'
#' @param model trained torch nn_module
#' @param data data.frame
#' @param prep preprocessing metadata
#' @param key character vector of column names
#' @param type variable types (or NULL)
#' @param emb_original numeric matrix of original embeddings (n x latent_dim)
#' @return named numeric vector of importance scores summing to 1
#' @keywords internal
.ae_var_importance <- function(model, data, prep, key, type = NULL,
                               emb_original = NULL) {
  if (is.null(emb_original)) {
    emb_original <- .ae_encode(model, data, prep, key, type)
  }
  n <- nrow(data)
  importance <- setNames(numeric(length(key)), key)

  for (v in key) {
    shuffled <- data
    shuffled[[v]] <- data[[v]][sample.int(n)]
    emb_shuffled <- .ae_encode(model, shuffled, prep, key, type)
    # Mean Euclidean shift
    shifts <- sqrt(rowSums((emb_original - emb_shuffled)^2))
    importance[v] <- mean(shifts)
  }

  # Normalize to sum to 1 (avoid division by zero)
  total <- sum(importance)
  if (total > 0) importance <- importance / total

  importance
}


#' Orchestrate embedding-based linkage for one block
#'
#' Preprocesses, trains autoencoder, encodes query and search data,
#' computes normalized distances, and variable importance.
#'
#' @param query data.frame of query records
#' @param search data.frame of search records
#' @param key character vector of column names
#' @param type variable types (or NULL)
#' @param latent_dim integer or NULL (auto)
#' @param epochs integer, max training epochs
#' @return list with:
#'   \describe{
#'     \item{dist_mat}{numeric matrix (n_q x n_s) of normalized distances}
#'     \item{var_importance}{named numeric vector}
#'     \item{threshold}{numeric, distance normalization threshold}
#'     \item{latent_dim}{integer, latent dimension used}
#'   }
#' @keywords internal
.embedding_linkage_block <- function(query, search, key, type = NULL,
                                     latent_dim = NULL, epochs = 50L) {
  # Train on query data only (attacker model)
  trained <- .ae_train(query, key, type,
                       latent_dim = latent_dim, epochs = epochs)

  # Encode both datasets
  emb_query <- .ae_encode(trained$model, query, trained$prep, key, type)
  emb_search <- .ae_encode(trained$model, search, trained$prep, key, type)

  # Distances
  dist_res <- .ae_distance(emb_query, emb_search)

  # Variable importance
  var_imp <- .ae_var_importance(trained$model, query, trained$prep, key, type,
                                emb_original = emb_query)

  list(
    dist_mat = dist_res$dist_mat,
    var_importance = var_imp,
    threshold = dist_res$threshold,
    latent_dim = trained$latent_dim
  )
}
```

- [ ] **Step 2: Write variable importance and orchestrator tests**

Append to test file:

```r
test_that(".ae_var_importance returns valid importance scores", {
  skip_if_not_installed("torch")
  set.seed(42)
  df <- data.frame(a = rnorm(30), b = sample(letters[1:3], 30, TRUE),
                   stringsAsFactors = FALSE)
  trained <- .ae_train(df, key = c("a", "b"), latent_dim = 2L, epochs = 5L)
  emb <- .ae_encode(trained$model, df, trained$prep, key = c("a", "b"))
  vi <- .ae_var_importance(trained$model, df, trained$prep,
                           key = c("a", "b"), emb_original = emb)

  expect_equal(names(vi), c("a", "b"))
  expect_true(all(vi >= 0))
  expect_equal(sum(vi), 1, tolerance = 1e-10)
})

test_that(".embedding_linkage_block returns expected structure", {
  skip_if_not_installed("torch")
  set.seed(42)
  X <- data.frame(a = rnorm(25), b = rnorm(25))
  Y <- data.frame(a = rnorm(25), b = rnorm(25))
  res <- .embedding_linkage_block(X, Y, key = c("a", "b"), epochs = 3L)

  expect_equal(nrow(res$dist_mat), 25)
  expect_equal(ncol(res$dist_mat), 25)
  expect_true(all(res$dist_mat >= 0))
  expect_true(all(res$dist_mat <= 1))
  expect_true(is.numeric(res$var_importance))
  expect_true(res$threshold > 0)
  expect_true(is.integer(res$latent_dim))
})
```

- [ ] **Step 3: Run tests**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-recordLinkage-embedding.R')"`
Expected: 15 tests pass

- [ ] **Step 4: Commit**

```bash
git add R/embedding_internal.R tests/testthat/test-recordLinkage-embedding.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add .ae_var_importance() and .embedding_linkage_block() orchestrator"
```

---

## Chunk 2: Integration into recordLinkage()

### Task 6: Wire Embedding Method into recordLinkage()

**Files:**
- Modify: `R/recordLinkage.R` (6 locations)
- Modify: `DESCRIPTION`
- Modify: `tests/testthat/test-recordLinkage-embedding.R`

This is the largest task. It adds the "embedding" method to the `recordLinkage()` function signature, validation, main processing loop, variable importance, output settings, and print/summary/plot display.

**Integration points** (line numbers from current `R/recordLinkage.R`):

1. **Function signature** (line 507-509): Add `"embedding"` to method choices
2. **Validation block** (after line 589): Add torch dependency check + strategy warning
3. **Main loop** (after line 1317): Add `} else if (method == "embedding") {` block
4. **Variable importance** (after line 1398): Add embedding case
5. **Settings** (line 1465-1469): Add emb_* settings
6. **Print method** (after line 2308): Add embedding-specific display
7. **Summary method** (line 2476-2482): Add embedding importance label
8. **Plot method** (line 2606-2612): Add embedding importance label in vi_label switch

- [ ] **Step 1: Add `torch` to DESCRIPTION Suggests**

In `DESCRIPTION`, add `torch` after `testthat`:

```
Suggests:
    ...
    testthat (>= 3.0.0),
    torch,
    knitr,
    ...
```

- [ ] **Step 2: Add "embedding" to method choices and new params to signature**

In `R/recordLinkage.R` line 507-509, change:

```r
method = c("deterministic", "probabilistic",
           "pram", "predictive", "rf",
           "rbrl", "mahalanobis"),
```

to:

```r
method = c("deterministic", "probabilistic",
           "pram", "predictive", "rf",
           "rbrl", "mahalanobis", "embedding"),
```

After `ot_max_iter = 100L,` (line 541), add three new parameters:

```r
emb_latent_dim = NULL,
emb_epochs = 50L,
emb_global = FALSE,
```

- [ ] **Step 3: Add validation block for embedding method**

After the RF validation block (line 589), add:

```r
if (method == "embedding") {
    if (!requireNamespace("torch", quietly = TRUE)) {
        stop("Package 'torch' required for recordLinkage(method = 'embedding'). ",
             "Install with install.packages('torch')", call. = FALSE)
    }
    if (!missing(weights)) {
        message("method = 'embedding': weights ignored. ",
                "Embedding uses learned variable representations instead.")
    }
    if (!missing(strategy) && strategy != "nearest") {
        message("method = 'embedding': strategy '", strategy,
                "' uses distance-based scores [0,1], ",
                "not threshold-based matching.")
    }
}
```

- [ ] **Step 4: Add embedding processing block in main loop**

After the RF block closing `}` (line 1317), add:

```r
} else if (method == "embedding") {

    # Embedding method: autoencoder-based record linkage ----
    emb_var_importance <- NULL
    emb_actual_latent_dim <- NULL

    .emb_process_block <- function(q_idx, s_idx, dist_mat) {
        for (r in seq_along(q_idx)) {
            qi <- q_idx[r]
            dist_vec <- dist_mat[r, ]
            cand_n[qi] <<- length(s_idx)
            best_col <- which.min(dist_vec)
            risk[qi] <<- 1 - dist_vec[best_col]  # lower distance = higher risk

            # Truth evaluation
            tpos <- true_idx[qi]
            if (!is.na(tpos) && tpos > 0L) {
                t_col <- match(tpos, s_idx)
                if (!is.na(t_col)) {
                    true_in_set[qi] <<- TRUE
                    d_true[qi] <<- dist_vec[t_col]
                    d_min[qi] <<- min(dist_vec)
                    d_rank[qi] <<- as.integer(sum(dist_vec <= dist_vec[t_col]))
                }
            }

            if (!is.null(score_cache)) {
                score_cache[[qi]] <<- list(cand = s_idx, scores = dist_vec,
                                           maximize = FALSE)
            }

            if (isTRUE(return_matches)) {
                matches[[qi]] <<- s_idx[best_col]
            }
        }
    }

    if (is.null(block) || length(split_search) == 1) {
        # No blocking: single autoencoder
        block_res <- .embedding_linkage_block(
            query_data, search_data, key, type,
            latent_dim = emb_latent_dim, epochs = emb_epochs
        )
        s_idx <- seq_len(nrow(search_data))
        .emb_process_block(seq_len(n_query), s_idx, block_res$dist_mat)
        emb_var_importance <- block_res$var_importance
        emb_actual_latent_dim <- block_res$latent_dim

    } else {
        # Blocked matching
        all_blocks <- names(split_search)
        small_blocks <- 0L
        total_blocks <- length(all_blocks)
        fallback_blocks <- character(0)
        importance_list <- list()
        block_sizes <- integer(0)

        # Determine latent_dim for fallback threshold
        if (is.null(emb_latent_dim)) {
            # Estimate from full data
            prep_est <- .ae_preprocess(query_data, key, type)
            auto_latent <- as.integer(max(2L, floor(prep_est$input_dim / 3)))
        } else {
            auto_latent <- emb_latent_dim
        }
        min_block_size <- as.integer(max(30L, 5L * auto_latent))

        if (emb_global) {
            # Global embedding: train once on all query data
            trained <- .ae_train(query_data, key, type,
                                 latent_dim = emb_latent_dim,
                                 epochs = emb_epochs)
            emb_actual_latent_dim <- trained$latent_dim

            emb_all_query <- .ae_encode(trained$model, query_data,
                                        trained$prep, key, type)
            emb_all_search <- .ae_encode(trained$model, search_data,
                                         trained$prep, key, type)

            # Variable importance from global model
            emb_var_importance <- .ae_var_importance(
                trained$model, query_data, trained$prep, key, type,
                emb_original = emb_all_query
            )

            for (blk in all_blocks) {
                s_idx <- split_search[[blk]]
                q_idx <- which(blk_query == blk)
                if (length(q_idx) == 0 || length(s_idx) == 0) next

                # Within-block distances from global embeddings
                dist_res <- .ae_distance(emb_all_query[q_idx, , drop = FALSE],
                                         emb_all_search[s_idx, , drop = FALSE])
                .emb_process_block(q_idx, s_idx, dist_res$dist_mat)
            }

        } else {
            # Per-block embedding
            for (blk in all_blocks) {
                s_idx <- split_search[[blk]]
                q_idx <- which(blk_query == blk)
                if (length(q_idx) == 0 || length(s_idx) == 0) next

                if (length(q_idx) < min_block_size) {
                    small_blocks <- small_blocks + 1L
                    fallback_blocks <- c(fallback_blocks, blk)
                    next
                }

                block_res <- .embedding_linkage_block(
                    query_data[q_idx, , drop = FALSE],
                    search_data[s_idx, , drop = FALSE],
                    key, type,
                    latent_dim = emb_latent_dim, epochs = emb_epochs
                )
                .emb_process_block(q_idx, s_idx, block_res$dist_mat)
                importance_list[[blk]] <- block_res$var_importance
                block_sizes <- c(block_sizes, length(q_idx))
                if (is.null(emb_actual_latent_dim)) {
                    emb_actual_latent_dim <- block_res$latent_dim
                }
            }

            # Fallback blocks: use deterministic (Gower distance)
            if (length(fallback_blocks) > 0) {
                for (blk in fallback_blocks) {
                    s_idx <- split_search[[blk]]
                    q_idx <- which(blk_query == blk)
                    if (length(q_idx) == 0 || length(s_idx) == 0) next

                    anon_block <- search_data[s_idx, , drop = FALSE]
                    for (r in seq_along(q_idx)) {
                        qi <- q_idx[r]
                        di <- .dist_to_candidates(
                            query_data[qi, , drop = FALSE],
                            anon_block, key, type, weights,
                            wsum, rng, na_anon
                        )
                        cand_n[qi] <- length(s_idx)
                        best <- which.min(di)
                        risk[qi] <- 1 - di[best]

                        tpos <- true_idx[qi]
                        if (!is.na(tpos) && tpos > 0L) {
                            t_col <- match(tpos, s_idx)
                            if (!is.na(t_col)) {
                                true_in_set[qi] <- TRUE
                                d_true[qi] <- di[t_col]
                                d_min[qi] <- min(di)
                                d_rank[qi] <- as.integer(
                                    sum(di <= di[t_col]))
                            }
                        }
                        if (!is.null(score_cache)) {
                            score_cache[[qi]] <- list(
                                cand = s_idx, scores = di,
                                maximize = FALSE)
                        }
                        if (isTRUE(return_matches)) {
                            matches[[qi]] <- s_idx[best]
                        }
                    }
                }
            }

            # Aggregate variable importance
            if (length(importance_list) > 0) {
                imp_mat <- do.call(rbind, importance_list)
                emb_var_importance <- colSums(imp_mat * block_sizes) /
                    sum(block_sizes)
            }

            if (small_blocks > 0) {
                message(small_blocks, " of ", total_blocks,
                        " blocks have < ", min_block_size,
                        " query records; using deterministic fallback.")
            }
        }
    }

}
```

**Important:** The `} else if (method == "embedding") {` block goes right after the RF block's final closing `}` at line 1317 and before the bijective matching override section at line 1319. The existing code at line 1317 has `}` which closes `} else if (method == "rf") {`. The new code must be inserted so the structure becomes:

```
    } else if (method == "rf") {
        ...
    } else if (method == "embedding") {
        ...
    }

    # bijective matching override ----
```

- [ ] **Step 5: Add variable importance case**

After the mahalanobis case (line 1398), before the `else` block (line 1399), insert:

```r
} else if (method == "embedding") {
    var_importance <- if (!is.null(emb_var_importance)) {
        emb_var_importance
    } else {
        setNames(rep(NA_real_, length(key)), key)
    }
```

- [ ] **Step 6: Add embedding settings to output**

After line 1467 (`robust = if (method == "mahalanobis") robust else NULL,`), add:

```r
emb_latent_dim = if (method == "embedding") {
    emb_actual_latent_dim
} else NULL,
emb_epochs = if (method == "embedding") emb_epochs else NULL,
emb_global = if (method == "embedding") emb_global else NULL,
```

- [ ] **Step 7: Add embedding display to print method**

After the predictive-specific print block (line 2308), add:

```r
if (meth == "embedding") {
    ldim <- if (!is.null(s$emb_latent_dim)) s$emb_latent_dim else "?"
    epc <- if (!is.null(s$emb_epochs)) s$emb_epochs else "?"
    cat("Embedding:   autoencoder (latent_dim = ", ldim,
        ", epochs = ", epc, ")\n", sep = "")
}
```

- [ ] **Step 8: Add embedding label to summary variable importance switch**

In the summary method, in the `vi_label <- switch(...)` block (line 2476-2482), add before the default fallback:

```r
embedding   = "Variable Importance (permutation-based embedding shift):",
```

So the switch becomes:

```r
vi_label <- switch(x$method,
    deterministic = "Variable Importance (mean weighted distance):",
    probabilistic = "Variable Importance (log-LR on agreement):",
    predictive    = "Variable Importance (model coefficients):",
    pram          = "Variable Importance (perturbation strength):",
    mahalanobis   = "Variable Importance (precision matrix proportion):",
    embedding     = "Variable Importance (permutation-based embedding shift):",
    "Variable Importance:")
```

- [ ] **Step 8b: Update plot method vi_label switch**

In the `plot.recordLinkageRisk` method, in the `vi_label <- switch(...)` block (line 2606-2612), add before the default fallback:

```r
embedding     = "Permutation Embedding Shift",
```

So the switch becomes:

```r
vi_label <- switch(x$method,
    deterministic = "Mean Weighted Distance",
    predictive    = "Absolute Model Coefficient",
    pram          = "Perturbation Strength (1 - diag)",
    rf            = "RF Impurity Importance",
    mahalanobis   = "Precision Matrix Proportion",
    embedding     = "Permutation Embedding Shift",
    "Importance")
```

- [ ] **Step 9: Update roxygen documentation**

Update `@param method` (line 241-246) to include `"embedding"`:

```r
#' @param method character. Linkage method: \code{"deterministic"} (default),
#'   \code{"probabilistic"} (Fellegi-Sunter), \code{"pram"} (transition matrix),
#'   \code{"predictive"} (propensity-score-based), \code{"rf"}
#'   (random forest proximity-based; requires \pkg{ranger}),
#'   \code{"rbrl"} (rank-based record linkage),
#'   \code{"mahalanobis"} (Mahalanobis distance with robust covariance), or
#'   \code{"embedding"} (autoencoder latent-space distance; requires \pkg{torch}).
```

Add new `@param` entries (after existing `rf_global` or `ot_max_iter` params):

```r
#' @param emb_latent_dim integer or NULL. For \code{method = "embedding"}:
#'   bottleneck dimension. If \code{NULL} (default), auto-computed as
#'   \code{max(2, floor(input_dim / 3))} where \code{input_dim} includes
#'   entity embedding dimensions for categoricals.
#' @param emb_epochs integer. For \code{method = "embedding"}: maximum
#'   training epochs (default 50). Early stopping with patience 5 may
#'   terminate training earlier.
#' @param emb_global logical. For \code{method = "embedding"} with blocking:
#'   if \code{TRUE}, train a single autoencoder on all query data and
#'   restrict distance computation within blocks. If \code{FALSE} (default),
#'   train a separate autoencoder per block (with deterministic fallback for
#'   blocks smaller than \code{max(30, 5 * latent_dim)}).
```

Add `@section Embedding method:` (after existing method sections):

```r
#' @section Embedding method:
#' The \code{method = "embedding"} approach trains an autoencoder with entity
#' embeddings (Guo & Berkhahn, 2016) on the original (query) data, projects
#' both query and search records into a latent space, and measures
#' re-identification risk via Euclidean distance. This captures nonlinear
#' dependencies between quasi-identifiers that Gower and Mahalanobis distances
#' may miss.
#'
#' The autoencoder trains only on query data, modeling an attacker who knows
#' the population structure but not the specific anonymization. Distances are
#' normalized to [0,1] using the 97.5th percentile of within-query pairwise
#' distances. Variable importance is permutation-based: each variable is
#' shuffled and the mean embedding shift is measured.
#'
#' Requires the \pkg{torch} package. Entity embeddings handle mixed data
#' naturally; all-numeric and all-categorical datasets are supported.
```

- [ ] **Step 10: Write integration tests**

Append to test file:

```r
# --- Integration tests via recordLinkage() ---

test_that("recordLinkage(method = 'embedding') basic numeric", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       emb_epochs = 5L)

  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(nrow(res$per_record), 30)
  expect_true(all(res$per_record$risk >= 0))
  expect_true(all(res$per_record$risk <= 1))
  expect_equal(res$settings$matching, "independent")
})

test_that("recordLinkage(method = 'embedding') all-categorical", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(
    sex = sample(c("M", "F"), 30, TRUE),
    edu = sample(c("low", "mid", "high"), 30, TRUE),
    stringsAsFactors = FALSE
  )
  Y <- data.frame(
    sex = sample(c("M", "F"), 30, TRUE),
    edu = sample(c("low", "mid", "high"), 30, TRUE),
    stringsAsFactors = FALSE
  )
  res <- recordLinkage(X, Y, key = c("sex", "edu"), method = "embedding",
                       emb_epochs = 5L)

  expect_s3_class(res, "recordLinkageRisk")
  expect_true(all(is.finite(res$per_record$risk)))
})

test_that("recordLinkage(method = 'embedding') mixed data", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(a = rnorm(30),
                  b = sample(c("X", "Y", "Z"), 30, TRUE),
                  stringsAsFactors = FALSE)
  Y <- data.frame(a = rnorm(30),
                  b = sample(c("X", "Y", "Z"), 30, TRUE),
                  stringsAsFactors = FALSE)
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       emb_epochs = 5L)

  expect_s3_class(res, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'embedding') near-copy detection", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- X + rnorm(40, 0, 0.01)  # near-copies
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       truth = "row", emb_epochs = 10L)

  expect_true(res$overall$mean_risk > 0.1)
})

test_that("recordLinkage(method = 'embedding') bijective matching", {
  skip_if_not_installed("torch")
  skip_if_not_installed("clue")
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       matching = "bijective", emb_epochs = 5L)

  expect_s3_class(res, "recordLinkageRisk")
  expect_true(all(res$per_record$risk %in% c(0, 1)))
})

test_that("recordLinkage(method = 'embedding') OT matching", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       matching = "ot", emb_epochs = 5L)

  expect_s3_class(res, "recordLinkageRisk")
  expect_true(any(res$per_record$risk > 0 & res$per_record$risk < 1))
})

test_that("recordLinkage(method = 'embedding') with blocking", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(a = rnorm(80), b = rep(c("A", "B"), 40))
  Y <- data.frame(a = rnorm(80), b = rep(c("A", "B"), 40))
  res <- recordLinkage(X, Y, key = c("a", "b"), block = "b",
                       method = "embedding", emb_epochs = 3L)

  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(nrow(res$per_record), 80)
})

test_that("recordLinkage(method = 'embedding') emb_global = TRUE with blocking", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(a = rnorm(80), b = rep(c("A", "B"), 40))
  Y <- data.frame(a = rnorm(80), b = rep(c("A", "B"), 40))
  res <- recordLinkage(X, Y, key = c("a", "b"), block = "b",
                       method = "embedding", emb_epochs = 3L,
                       emb_global = TRUE)

  expect_s3_class(res, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'embedding') direction = 'reverse'", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       direction = "reverse", emb_epochs = 5L)

  expect_equal(res$direction, "reverse")
})

test_that("recordLinkage(method = 'embedding') custom latent_dim and epochs", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       emb_latent_dim = 4L, emb_epochs = 3L)

  expect_equal(res$settings$emb_latent_dim, 4L)
  expect_equal(res$settings$emb_epochs, 3L)
})

test_that("recordLinkage(method = 'embedding') small-block fallback", {
  skip_if_not_installed("torch")
  set.seed(1)
  # Block C has only 5 records -> should fall back to deterministic
  X <- data.frame(a = rnorm(85),
                  b = c(rep("A", 40), rep("B", 40), rep("C", 5)))
  Y <- data.frame(a = rnorm(85),
                  b = c(rep("A", 40), rep("B", 40), rep("C", 5)))
  expect_message(
    res <- recordLinkage(X, Y, key = c("a", "b"), block = "b",
                         method = "embedding", emb_epochs = 3L),
    "deterministic fallback"
  )
  expect_s3_class(res, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'embedding') unseen categories", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(a = rnorm(30),
                  b = sample(c("X", "Y"), 30, TRUE),
                  stringsAsFactors = FALSE)
  Y <- data.frame(a = rnorm(30),
                  b = sample(c("X", "Y", "Z"), 30, TRUE),  # Z unseen
                  stringsAsFactors = FALSE)
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       emb_epochs = 5L)

  expect_s3_class(res, "recordLinkageRisk")
  expect_true(all(is.finite(res$per_record$risk)))
})

test_that("recordLinkage(method = 'embedding') print output", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       emb_epochs = 3L)

  out <- capture.output(print(res))
  expect_true(any(grepl("embedding", out)))
  expect_true(any(grepl("autoencoder", out)))
  expect_true(any(grepl("latent_dim", out)))
})

test_that("recordLinkage(method = 'embedding') summary output", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       emb_epochs = 3L)

  out <- capture.output(print(summary(res)))
  expect_true(any(grepl("permutation", out, ignore.case = TRUE)))
})

test_that("recordLinkage(method = 'embedding') variable importance", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  res <- recordLinkage(X, Y, key = c("a", "b", "c"), method = "embedding",
                       emb_epochs = 5L)

  expect_true(!is.null(res$var_importance))
  expect_equal(length(res$var_importance), 3)
  expect_true(all(res$var_importance >= 0))
  expect_equal(sum(res$var_importance), 1, tolerance = 0.01)
})

test_that("recordLinkage(method = 'embedding') weights message", {
  skip_if_not_installed("torch")
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  expect_message(
    recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                  weights = c(a = 1, b = 2), emb_epochs = 3L),
    "importance"
  )
})

test_that("recordLinkage(method = 'embedding') reproducibility", {
  skip_if_not_installed("torch")
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))

  set.seed(42)
  res1 <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                        emb_epochs = 5L)

  set.seed(42)
  res2 <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                        emb_epochs = 5L)

  expect_equal(res1$per_record$risk, res2$per_record$risk)
})
```

- [ ] **Step 11: Run tests**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-recordLinkage-embedding.R')"`
Expected: All ~33 tests pass (or skip if torch not installed)

- [ ] **Step 12: Run devtools::document() to update NAMESPACE and .Rd**

Run: `Rscript -e "devtools::document()"`
Expected: Updated man/recordLinkage.Rd with new params and sections

- [ ] **Step 13: Run R CMD check**

Run: `Rscript -e "devtools::check(args = c('--no-manual', '--no-vignettes'))"`
Expected: 0 errors, 0 warnings (notes OK)

- [ ] **Step 14: Commit**

```bash
git add DESCRIPTION R/recordLinkage.R R/embedding_internal.R \
  tests/testthat/test-recordLinkage-embedding.R man/ NAMESPACE
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Wire method='embedding' into recordLinkage() with full integration"
```

---

## Chunk 3: Verification and Cleanup

### Task 7: Full Test Suite and R CMD Check

**Files:**
- All modified files from Tasks 1-6

- [ ] **Step 1: Run complete test suite**

Run: `Rscript -e "devtools::test()"`
Expected: All existing tests still pass, all new embedding tests pass

- [ ] **Step 2: Run R CMD check**

Run: `_R_CHECK_FORCE_SUGGESTS_=FALSE R CMD build . && _R_CHECK_FORCE_SUGGESTS_=FALSE R CMD check riskutility_*.tar.gz --no-manual --no-vignettes`
Expected: 0 errors, 0 notes (2 expected vignette warnings)

- [ ] **Step 3: Fix any issues found**

Address any failing tests or check warnings.

- [ ] **Step 4: Final commit if fixes were needed**

```bash
git add -A
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Fix issues found during full test suite verification"
```

---

## Summary

| Task | Component | New/Modified | Lines |
|------|-----------|-------------|-------|
| 1 | `.ae_preprocess()` | Create `R/embedding_internal.R` | ~100 |
| 2 | `.ae_model()` | Append to `R/embedding_internal.R` | ~90 |
| 3 | `.ae_train()` | Append to `R/embedding_internal.R` | ~120 |
| 4 | `.ae_encode()`, `.ae_distance()` | Append to `R/embedding_internal.R` | ~70 |
| 5 | `.ae_var_importance()`, `.embedding_linkage_block()` | Append to `R/embedding_internal.R` | ~60 |
| 6 | Integration into `recordLinkage()` | Modify `R/recordLinkage.R` + `DESCRIPTION` | ~200 |
| 7 | Verification | All files | — |

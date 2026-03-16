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

    if (!is.null(out$num_recon) && !is.null(num_target)) {
      losses$mse <- torch::nnf_mse_loss(out$num_recon, num_target)
    }

    if (!is.null(out$cat_logits) && length(out$cat_logits) > 0) {
      ce_sum <- torch::torch_tensor(0, dtype = torch::torch_float())
      for (k in seq_along(prep$cat_keys)) {
        ce_sum <- ce_sum + torch::nnf_cross_entropy(
          out$cat_logits[[k]], cat_targets[[k]]
        )
      }
      losses$ce <- ce_sum / length(prep$cat_keys)
    }

    if (length(losses) == 0) return(torch::torch_tensor(0))

    # Weight by proportion of each type
    n_vars <- length(prep$num_keys) + length(prep$cat_keys)
    total <- torch::torch_tensor(0, dtype = torch::torch_float())
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
      mean_emb <- if (nl == 1L) learned else learned$mean(dim = 1)
      model$embeddings[[k]]$weight[nl + 1L, ] <- mean_emb
    }
  })

  list(
    model = model,
    prep = prep,
    latent_dim = latent_dim
  )
}

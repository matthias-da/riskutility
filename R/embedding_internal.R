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

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

test_that(".ae_preprocess handles constant numeric column", {
  skip_if_not_installed("torch")
  df <- data.frame(a = c(5, 5, 5), b = c(1, 2, 3))
  res <- .ae_preprocess(df, key = c("a", "b"), type = NULL)

  # Constant column -> 0.5
  expect_equal(res$num_mat[, "a"], c(0.5, 0.5, 0.5))
  # Normal column still scales correctly
  expect_equal(res$num_mat[, "b"], c(0, 0.5, 1))
})

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

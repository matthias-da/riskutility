# tests/testthat/test-recordLinkage-embedding.R

test_that(".ae_preprocess handles numeric data", {
  skip_if_no_torch()
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
  skip_if_no_torch()
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
  skip_if_no_torch()
  df <- data.frame(a = c(1, 2, 3), b = c("X", "Y", "X"),
                   stringsAsFactors = FALSE)
  res <- .ae_preprocess(df, key = c("a", "b"), type = NULL)

  expect_equal(res$num_keys, "a")
  expect_equal(res$cat_keys, "b")
  # input_dim = 1 (numeric) + emb_dim(b)
  expect_equal(res$input_dim, 1L + res$emb_dims[["b"]])
})

test_that(".ae_preprocess maps unseen levels to UNK", {
  skip_if_no_torch()
  df_train <- data.frame(x = c("A", "B"), stringsAsFactors = FALSE)
  prep_train <- .ae_preprocess(df_train, key = "x", type = NULL)

  df_new <- data.frame(x = c("A", "C"), stringsAsFactors = FALSE)
  prep_new <- .ae_preprocess(df_new, key = "x", type = NULL,
                             cat_maps = prep_train$cat_maps)

  # A -> 1, C (unseen) -> n_levels + 1 = 3
  expect_equal(prep_new$cat_idx[["x"]], c(1L, 3L))
})

test_that(".ae_preprocess handles constant numeric column", {
  skip_if_no_torch()
  df <- data.frame(a = c(5, 5, 5), b = c(1, 2, 3))
  res <- .ae_preprocess(df, key = c("a", "b"), type = NULL)

  # Constant column -> 0.5
  expect_equal(res$num_mat[, "a"], c(0.5, 0.5, 0.5))
  # Normal column still scales correctly
  expect_equal(res$num_mat[, "b"], c(0, 0.5, 1))
})

test_that(".ae_model constructs valid autoencoder", {
  skip_if_no_torch()
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
  skip_if_no_torch()
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
  skip_if_no_torch()
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
  skip_if_no_torch()
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
  skip_if_no_torch()
  set.seed(42)
  df <- data.frame(a = rnorm(30), b = rnorm(30))
  res <- .ae_train(df, key = c("a", "b"), latent_dim = 4L, epochs = 3L)

  expect_equal(res$latent_dim, 4L)
})

test_that(".ae_train is reproducible with set.seed", {
  skip_if_no_torch()
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
  skip_if_no_torch()
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
  skip_if_no_torch()
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

test_that(".ae_var_importance returns valid importance scores", {
  skip_if_no_torch()
  set.seed(42)
  df <- data.frame(a = rnorm(30), b = sample(letters[1:3], 30, TRUE),
                   stringsAsFactors = FALSE)
  trained <- .ae_train(df, key = c("a", "b"), latent_dim = 2L, epochs = 5L)
  emb <- .ae_encode(trained$model, df, trained$prep, key = c("a", "b"))
  vi <- .ae_var_importance(trained$model, df, trained$prep,
                           key = c("a", "b"), emb_original = emb)

  expect_equal(names(vi), c("a", "b"))
  expect_true(all(vi >= 0))
  # Sum to 1 unless all shifts are zero (degenerate model)
  if (sum(vi) > 0) expect_equal(sum(vi), 1, tolerance = 1e-10)
})

test_that(".embedding_linkage_block returns expected structure", {
  skip_if_no_torch()
  set.seed(42)
  X <- data.frame(a = rnorm(25), b = rnorm(25))
  Y <- data.frame(a = rnorm(25), b = rnorm(25))
  res <- .embedding_linkage_block(X, Y, key = c("a", "b"), epochs = 3L)

  expect_equal(nrow(res$dist_mat), 25)
  expect_equal(ncol(res$dist_mat), 25)
  expect_true(all(res$dist_mat >= 0))
  expect_true(all(res$dist_mat <= 1))
  expect_true(is.numeric(res$var_importance))
  expect_equal(names(res$var_importance), c("a", "b"))
  expect_true(res$threshold > 0)
  expect_true(is.integer(res$latent_dim))
})

# --- Integration tests via recordLinkage() ---

test_that("recordLinkage(method = 'embedding') basic numeric", {
  skip_if_no_torch()
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
  skip_if_no_torch()
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
  skip_if_no_torch()
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
  skip_if_no_torch()
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- X + rnorm(40, 0, 0.01)  # near-copies
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       truth = "row", emb_epochs = 10L)

  expect_true(res$overall$mean_risk > 0.1)
})

test_that("recordLinkage(method = 'embedding') bijective matching", {
  skip_if_no_torch()
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
  skip_if_no_torch()
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       matching = "ot", emb_epochs = 5L)

  expect_s3_class(res, "recordLinkageRisk")
  expect_true(any(res$per_record$risk > 0 & res$per_record$risk < 1))
})

test_that("recordLinkage(method = 'embedding') with blocking", {
  skip_if_no_torch()
  set.seed(1)
  X <- data.frame(a = rnorm(80), b = rep(c("A", "B"), 40))
  Y <- data.frame(a = rnorm(80), b = rep(c("A", "B"), 40))
  res <- recordLinkage(X, Y, key = c("a", "b"), block = "b",
                       method = "embedding", emb_epochs = 3L)

  expect_s3_class(res, "recordLinkageRisk")
  expect_equal(nrow(res$per_record), 80)
})

test_that("recordLinkage(method = 'embedding') emb_global = TRUE with blocking", {
  skip_if_no_torch()
  set.seed(1)
  X <- data.frame(a = rnorm(80), b = rep(c("A", "B"), 40))
  Y <- data.frame(a = rnorm(80), b = rep(c("A", "B"), 40))
  res <- recordLinkage(X, Y, key = c("a", "b"), block = "b",
                       method = "embedding", emb_epochs = 3L,
                       emb_global = TRUE)

  expect_s3_class(res, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'embedding') direction = 'reverse'", {
  skip_if_no_torch()
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       direction = "reverse", emb_epochs = 5L)

  expect_equal(res$direction, "reverse")
})

test_that("recordLinkage(method = 'embedding') custom latent_dim and epochs", {
  skip_if_no_torch()
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       emb_latent_dim = 4L, emb_epochs = 3L)

  expect_equal(res$settings$emb_latent_dim, 4L)
  expect_equal(res$settings$emb_epochs, 3L)
})

test_that("recordLinkage(method = 'embedding') small-block fallback", {
  skip_if_no_torch()
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
  skip_if_no_torch()
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
  skip_if_no_torch()
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
  skip_if_no_torch()
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                       emb_epochs = 3L)

  out <- capture.output(print(summary(res)))
  expect_true(any(grepl("permutation", out, ignore.case = TRUE)))
})

test_that("recordLinkage(method = 'embedding') variable importance", {
  skip_if_no_torch()
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  res <- recordLinkage(X, Y, key = c("a", "b", "c"), method = "embedding",
                       emb_epochs = 5L)

  expect_true(!is.null(res$var_importance))
  expect_equal(length(res$var_importance), 3)
  expect_true(all(res$var_importance >= 0))
})

test_that("recordLinkage(method = 'embedding') weights message", {
  skip_if_no_torch()
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- data.frame(a = rnorm(20), b = rnorm(20))
  expect_message(
    recordLinkage(X, Y, key = c("a", "b"), method = "embedding",
                  weights = c(a = 1, b = 2), emb_epochs = 3L),
    "ignored"
  )
})

test_that("recordLinkage(method = 'embedding') reproducibility", {
  skip_if_no_torch()
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

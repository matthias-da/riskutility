#' Compare Feature Importance between Real and Synthetic Data
#'
#' Compares feature importance rankings from models trained on real (X) and synthetic/anonymized (Y) datasets.
#'
#' This function supports three types of feature importance measures:
#' \itemize{
#'   \item \strong{Model-based:} Uses \code{caret::varImp}.
#'   \item \strong{Permutation-based:} Permutes each feature \code{nsim} times
#'     and measures the resulting drop in model performance (computed
#'     internally, no additional packages required). A prediction wrapper is
#'     required.
#'   \item \strong{SHAP-based:} Computes mean absolute Shapley values using the \code{kernelshap} package.
#' }
#'
#' @param X Original dataset (data.frame or data.table). Must include the target variable.
#' @param Y Synthetic/anonymized dataset (data.frame or data.table). Must include the target variable.
#' @param formula A formula specifying the target and predictor variables (e.g., \code{target ~ .}).
#' @param method Modeling method (e.g., "rf", "rpart", "glm", etc.).
#' @param importance_type Type of importance measure: "model", "permutation", or "shap".
#' @param pred_wrapper Prediction function wrapper required for permutation and SHAP importance.
#'   If \code{NULL} and \code{importance_type = "shap"}, a default wrapper is used which returns the probability
#'   of the first class for classification tasks, or raw predictions for regression.
#' @param metric Performance metric for permutation importance. For classification tasks, common choices include
#'   \itemize{
#'     \item \code{"Accuracy"}
#'     \item \code{"Kappa"}
#'     \item \code{"ROC AUC"} (if probabilities are provided)
#'   }
#'   For regression tasks, you can use:
#'   \itemize{
#'     \item \code{"RMSE"} (Root Mean Squared Error)
#'     \item \code{"MAE"} (Mean Absolute Error)
#'     \item \code{"R-squared"}
#'     \item \code{"MAPE"} (Mean Absolute Percentage Error)
#'   }
#'   Metric names are matched case-insensitively. Error metrics (RMSE, MAE,
#'   MAPE) are sign-flipped so that a larger importance always indicates a
#'   more influential feature. For \code{"ROC AUC"} the positive class is the
#'   second factor level and \code{pred_wrapper} must return probabilities.
#'   Additionally, a custom metric function can be supplied that takes the observed values and predictions as inputs and
#'   returns a single numeric performance value; it is assumed to be
#'   larger-is-better.
#' @param nsim Number of simulations used by permutation importance.
#'
#' @param ... additional arguments passed to methods
#'
#' @return A list containing:
#' \item{comparison}{A data.table with columns for feature names, importance values for datasets X and Y, and the absolute difference.}
#' \item{model_X}{The model trained on X.}
#' \item{model_Y}{The model trained on Y.}
#'
#' @details
#' \itemize{
#'   \item If the target variable in \code{X} is a factor, the target variable in \code{Y} is coerced to a factor with the same levels.
#'   \item For SHAP importance, the function uses the \code{kernelshap} package. The returned importance is the mean absolute SHAP value per feature.
#'   \item For classification tasks in the SHAP branch, if no \code{pred_wrapper} is provided, predictions are taken as the probability of the first class.
#' }
#'
#' @examples
#' # Model-based importance (Decision Tree) - basic example
#' set.seed(123)
#' X <- data.frame(
#'   age = rnorm(100, 40, 10),
#'   income = 50000 + rnorm(100, 0, 10000),
#'   outcome = factor(sample(c("Yes", "No"), 100, replace = TRUE))
#' )
#' Y <- data.frame(
#'   age = rnorm(100, 42, 12),
#'   income = 48000 + rnorm(100, 0, 12000),
#'   outcome = factor(sample(c("Yes", "No"), 100, replace = TRUE))
#' )
#' res_tree <- compare_feature_importance(X, Y, outcome ~ .,
#'                                        method = "rpart",
#'                                        importance_type = "model")
#'
#' \donttest{
#' # Permutation importance
#' if (requireNamespace("caret", quietly = TRUE)) {
#'   pred_fn <- function(object, newdata) predict(object, newdata, type = "raw")
#'   res_perm <- compare_feature_importance(
#'     X, Y, outcome ~ ., method = "rf",
#'     importance_type = "permutation",
#'     pred_wrapper = pred_fn, metric = "Accuracy", nsim = 5
#'   )
#' }
#' }
#'
#' @importFrom data.table data.table
#' @family comparison
#' @author Matthias Templ
#' @export
compare_feature_importance <- function(X, ...) {
  UseMethod("compare_feature_importance")
}

#' @rdname compare_feature_importance
#' @export
compare_feature_importance.synth_pair <- function(X, ...) {
  compare_feature_importance.default(X = X$original, Y = X$synthetic, ...)
}

#' @rdname compare_feature_importance
#' @export
compare_feature_importance.default <- function(X, Y, formula, method,
                                       importance_type = c("model", "permutation", "shap"),
                                       pred_wrapper = NULL, metric = "Accuracy", nsim = 5, ...) {
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required for compare_feature_importance(). Please install it.")
  }
  importance_type <- match.arg(importance_type)
  if (importance_type == "permutation") {
    if (is.null(pred_wrapper)) stop("pred_wrapper must be provided for permutation importance.")
    .permutation_metric(metric)  # validate the metric before fitting models
  }
  target_var <- all.vars(formula)[1]

  # Ensure consistent target variable type between X and Y.
  if (is.factor(X[[target_var]]) && !is.factor(Y[[target_var]])) {
    Y[[target_var]] <- factor(Y[[target_var]], levels = levels(X[[target_var]]))
  }
  if (is.factor(Y[[target_var]]) && !is.factor(X[[target_var]])) {
    X[[target_var]] <- factor(X[[target_var]], levels = levels(Y[[target_var]]))
  }

  ctrl <- caret::trainControl(method = "cv", number = 5)
  model_X <- caret::train(formula, data = X, method = method, trControl = ctrl)
  model_Y <- caret::train(formula, data = Y, method = method, trControl = ctrl)

  features <- setdiff(names(X), target_var)

  if (importance_type == "model") {
    vi_X <- caret::varImp(model_X)$importance
    vi_Y <- caret::varImp(model_Y)$importance

    # Extract importance values safely, handling missing features
    extract_importance <- function(vi, features) {
      col_name <- if ("Overall" %in% colnames(vi)) "Overall" else colnames(vi)[1]
      # Get available features in the importance matrix
      available <- intersect(features, rownames(vi))
      imp <- setNames(rep(NA_real_, length(features)), features)
      if (length(available) > 0) {
        imp[available] <- vi[available, col_name]
      }
      return(imp)
    }

    importance_X <- extract_importance(vi_X, features)
    importance_Y <- extract_importance(vi_Y, features)
  } else if (importance_type == "permutation") {
    importance_X <- .permutation_importance(model_X, X, target_var, features,
                                            metric, pred_wrapper, nsim)
    importance_Y <- .permutation_importance(model_Y, Y, target_var, features,
                                            metric, pred_wrapper, nsim)
  } else if (importance_type == "shap") {
    if (!requireNamespace("kernelshap", quietly = TRUE)) {
      stop("Package 'kernelshap' is required for SHAP importance. Please install it.")
    }
    # If no prediction wrapper is provided, define a default one.
    # For classification tasks, return the probability of the first class.
    if (is.null(pred_wrapper)) {
      pred_wrapper <- function(object, newdata) {
        if (is.factor(X[[target_var]])) {
          pred <- predict(object, newdata, type = "prob")
          return(pred[, 1])
        } else {
          predict(object, newdata)
        }
      }
    }
    pred_fun <- function(object, X, ...) pred_wrapper(object, X)
    feat_X <- as.data.frame(X)[, features, drop = FALSE]
    feat_Y <- as.data.frame(Y)[, features, drop = FALSE]
    # Background samples for the SHAP reference distribution (capped for speed).
    bg_X <- feat_X[sample.int(nrow(feat_X), min(nrow(feat_X), 100L)), , drop = FALSE]
    bg_Y <- feat_Y[sample.int(nrow(feat_Y), min(nrow(feat_Y), 100L)), , drop = FALSE]
    # Compute SHAP values using kernelshap (model-agnostic, available on CRAN).
    shap_X <- kernelshap::kernelshap(model_X, X = feat_X, bg_X = bg_X,
                                     pred_fun = pred_fun, verbose = FALSE)
    shap_Y <- kernelshap::kernelshap(model_Y, X = feat_Y, bg_X = bg_Y,
                                     pred_fun = pred_fun, verbose = FALSE)
    importance_X <- colMeans(abs(shap_X$S))
    importance_Y <- colMeans(abs(shap_Y$S))
  }

  # Replace NaN with NA for consistent handling
  importance_X[is.nan(importance_X)] <- NA_real_
  importance_Y[is.nan(importance_Y)] <- NA_real_

  comparison <- data.table::data.table(
    feature = features,
    importance_X = importance_X[features],
    importance_Y = importance_Y[features],
    difference = abs(importance_X[features] - importance_Y[features])
  )

  # Compute rank correlation (Spearman) to assess importance ranking similarity
  # Handle case where all values might be NA
  valid_pairs <- !is.na(importance_X[features]) & !is.na(importance_Y[features])
  if (sum(valid_pairs) >= 2) {
    rank_corr <- cor(importance_X[features], importance_Y[features],
                     method = "spearman", use = "complete.obs")
    pearson_corr <- cor(importance_X[features], importance_Y[features],
                        method = "pearson", use = "complete.obs")
  } else {
    rank_corr <- NA_real_
    pearson_corr <- NA_real_
    warning("Not enough valid feature importance pairs to compute correlation.")
  }

  result <- list(
    comparison = comparison,
    model_X = model_X,
    model_Y = model_Y,
    rank_correlation = rank_corr,
    pearson_correlation = pearson_corr,
    mean_abs_difference = mean(comparison$difference, na.rm = TRUE)
  )

  class(result) <- "compare_feature_importance"
  return(result)
}

#' Print method for compare_feature_importance objects
#'
#' @param x an object of class "compare_feature_importance"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.compare_feature_importance <- function(x, ...) {
  cat("Feature Importance Comparison\n")
  cat("=============================\n\n")
  cat("Rank Correlation (Spearman):", round(x$rank_correlation, 4), "\n")
  cat("Pearson Correlation:", round(x$pearson_correlation, 4), "\n")
  cat("Mean Absolute Difference:", round(x$mean_abs_difference, 4), "\n\n")
  cat("Feature-wise Comparison:\n")
  print(x$comparison)
  invisible(x)
}


#' Summary method for compare_feature_importance objects
#'
#' @param object an object of class "compare_feature_importance"
#' @param ... additional arguments (ignored)
#' @return An object of class "summary.compare_feature_importance"
#' @export
summary.compare_feature_importance <- function(object, ...) {
  comp <- object$comparison
  summ <- list(
    n_features = nrow(comp),
    rank_correlation = object$rank_correlation,
    pearson_correlation = object$pearson_correlation,
    mean_abs_difference = object$mean_abs_difference,
    max_abs_difference = max(comp$difference, na.rm = TRUE),
    top_discrepant = comp[order(-comp$difference), , drop = FALSE],
    comparison = comp
  )
  class(summ) <- "summary.compare_feature_importance"
  summ
}


#' Print method for summary.compare_feature_importance objects
#'
#' @param x an object of class "summary.compare_feature_importance"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.compare_feature_importance <- function(x, ...) {
  cat("Summary: Feature Importance Comparison\n")
  cat("=======================================\n\n")

  cat("Features compared:", x$n_features, "\n\n")

  cat("Correlation:\n")
  cat("  Spearman (rank): ", sprintf("%.4f", x$rank_correlation), "\n")
  cat("  Pearson:         ", sprintf("%.4f", x$pearson_correlation), "\n\n")

  cat("Importance Differences:\n")
  cat("  Mean absolute: ", sprintf("%.4f", x$mean_abs_difference), "\n")
  cat("  Max absolute:  ", sprintf("%.4f", x$max_abs_difference), "\n\n")

  cat("Feature-wise Comparison (sorted by discrepancy):\n")
  print(x$top_discrepant, row.names = FALSE)

  invisible(x)
}


#' Plot method for compare_feature_importance objects
#'
#' @param x an object of class "compare_feature_importance"
#' @param y not used
#' @param which integer, which plot to produce (default 1). Currently only 1 is available.
#' @param ... additional arguments passed to plotting functions
#' @importFrom graphics barplot par legend
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.compare_feature_importance <- function(x, y = NULL, which = 1, ...) {
  comp <- x$comparison

  # Create side-by-side bar plot
  mat <- rbind(comp$importance_X, comp$importance_Y)
  colnames(mat) <- comp$feature

  op <- par(mar = c(8, 4, 4, 2))
  on.exit(par(op))

  barplot(mat,
          beside = TRUE,
          col = c("steelblue", "coral"),
          main = paste("Feature Importance Comparison\nRank Corr =",
                       round(x$rank_correlation, 3)),
          ylab = "Importance",
          las = 2, ...)
  legend("topright", legend = c("Original (X)", "Synthetic (Y)"),
         fill = c("steelblue", "coral"))

  invisible(x)
}


# ------------------------------------------------------------------------------
# Internal permutation importance (replaces vip::vi_permute; vip was archived
# from CRAN on 2026-07-08). Resolves a caret-style metric name, or a custom
# metric function(obs, pred), to a scoring function plus its direction.
.permutation_metric <- function(metric) {
  if (is.function(metric)) {
    return(list(fun = metric, smaller_is_better = FALSE))
  }
  key <- gsub("[^a-z]", "", tolower(metric))
  fun <- switch(key,
    accuracy = function(obs, pred) mean(pred == obs),
    kappa = function(obs, pred) {
      obs <- factor(obs)
      tab <- table(factor(pred, levels = levels(obs)), obs)
      po <- sum(diag(tab)) / sum(tab)
      pe <- sum(rowSums(tab) * colSums(tab)) / sum(tab)^2
      if (pe >= 1) return(NA_real_)
      (po - pe) / (1 - pe)
    },
    rocauc = ,
    auc = function(obs, pred) {
      obs <- as.integer(factor(obs)) - 1L
      if (length(unique(obs)) != 2L) {
        stop("metric 'ROC AUC' requires a binary target.")
      }
      # Mann-Whitney formulation; positive class = second factor level
      r <- rank(as.numeric(pred))
      n1 <- sum(obs == 1L)
      n0 <- sum(obs == 0L)
      (sum(r[obs == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
    },
    rmse = function(obs, pred) sqrt(mean((as.numeric(pred) - as.numeric(obs))^2)),
    mae = function(obs, pred) mean(abs(as.numeric(pred) - as.numeric(obs))),
    mape = function(obs, pred) {
      obs <- as.numeric(obs)
      mean(abs((obs - as.numeric(pred)) / obs)) * 100
    },
    rsquared = ,
    rsq = function(obs, pred) {
      obs <- as.numeric(obs)
      pred <- as.numeric(pred)
      1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
    },
    stop("Unknown metric '", metric, "'. Supported metric names: 'Accuracy', ",
         "'Kappa', 'ROC AUC', 'RMSE', 'MAE', 'MAPE', 'R-squared'; ",
         "alternatively supply a custom metric function(obs, pred).")
  )
  list(fun = fun, smaller_is_better = key %in% c("rmse", "mae", "mape"))
}

# For each feature: permute the column nsim times, score predictions through
# pred_wrapper, and report the mean performance degradation relative to the
# unpermuted baseline. Error metrics are sign-flipped so that a larger
# importance always means a more influential feature.
.permutation_importance <- function(model, data, target_var, features,
                                    metric, pred_wrapper, nsim) {
  m <- .permutation_metric(metric)
  data <- as.data.frame(data)
  obs <- data[[target_var]]
  baseline <- m$fun(obs, pred_wrapper(model, data))
  vapply(features, function(f) {
    permuted <- vapply(seq_len(nsim), function(i) {
      perm <- data
      perm[[f]] <- sample(perm[[f]])
      m$fun(obs, pred_wrapper(model, perm))
    }, numeric(1))
    if (m$smaller_is_better) mean(permuted) - baseline else baseline - mean(permuted)
  }, numeric(1))
}

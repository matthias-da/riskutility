#' Compare Feature Importance between Real and Synthetic Data
#'
#' Compares feature importance rankings from models trained on real (X) and synthetic/anonymized (Y) datasets.
#'
#' This function supports three types of feature importance measures:
#' \itemize{
#'   \item \strong{Model-based:} Uses \code{caret::varImp}.
#'   \item \strong{Permutation-based:} Uses \code{vip::vi_permute}. A prediction wrapper is required.
#'   \item \strong{SHAP-based:} Computes mean absolute Shapley values using the \code{fastshap} package.
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
#'   Additionally, a custom metric function can be supplied that takes the observed values and predictions as inputs and
#'   returns a single numeric performance value.
#' @param nsim Number of simulations or permutations (also used by \code{fastshap::explain}).
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
#'   \item For SHAP importance, the function uses the \code{fastshap} package. The returned importance is the mean absolute SHAP value per feature.
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
#' # Permutation importance (requires vip package)
#' if (requireNamespace("vip", quietly = TRUE)) {
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
    if (is.null(pred_wrapper)) stop("pred_wrapper must be provided for permutation importance.")
    if (!requireNamespace("vip", quietly = TRUE)) {
      stop("Package 'vip' is required for this function. Please install it.")
    }
    imp_X <- vip::vi_permute(model_X, train = X, target = target_var, metric = metric,
                             pred_wrapper = pred_wrapper, nsim = nsim)
    importance_X <- setNames(imp_X$Importance, imp_X$Variable)
    imp_Y <- vip::vi_permute(model_Y, train = Y, target = target_var, metric = metric,
                             pred_wrapper = pred_wrapper, nsim = nsim)
    importance_Y <- setNames(imp_Y$Importance, imp_Y$Variable)
  } else if (importance_type == "shap") {
    if (!requireNamespace("fastshap", quietly = TRUE)) {
      stop("Package 'fastshap' is required for SHAP importance. Please install it.")
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
    # Compute SHAP values using fastshap.
    shap_X <- fastshap::explain(model_X,
                                X = as.data.frame(X)[, features, drop = FALSE],
                                pred_wrapper = pred_wrapper,
                                nsim = nsim
    )
    importance_X <- colMeans(abs(shap_X))
    shap_Y <- fastshap::explain(model_Y,
                                X = as.data.frame(Y)[, features, drop = FALSE],
                                pred_wrapper = pred_wrapper,
                                nsim = nsim
    )
    importance_Y <- colMeans(abs(shap_Y))
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

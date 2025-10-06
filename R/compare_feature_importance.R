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
#' library(data.table)
#' set.seed(123)
#' X <- data.table(
#'   age = rnorm(500, 40, 10),
#'   income = 50000 + rnorm(500, 0, 10000),
#'   default = factor(sample(c("Yes", "No"), 500, replace = TRUE))
#' )
#'
#' Y <- data.table(
#'   age = rnorm(500, 42, 12),
#'   income = 3000 + 1200 * rnorm(500, 40, 12) + rnorm(500, 0, 10000),
#'   default = sample(c("Yes", "No"), 500, TRUE)
#' )
#'
#' # Model-based importance (Decision Tree)
#' res_tree <- compare_feature_importance(X, Y, default ~ ., method = "rpart", importance_type = "model")
#' print(res_tree$comparison)
#'
#' # Permutation importance (Random Forest)
#' pred_wrapper_classif <- function(object, newdata) predict(object, newdata, type = "raw")
#' res_perm <- compare_feature_importance(
#'   X, Y, default ~ ., method = "rf", importance_type = "permutation",
#'   pred_wrapper = pred_wrapper_classif, metric = "Accuracy", nsim = 5
#' )
#' print(res_perm$comparison)
#'
#' # Permutation importance (Random Forest)
#' pred_wrapper_classif <- function(object, newdata) predict(object, newdata, type = "raw")
#' res_perm <- compare_feature_importance(
#'   X, Y, income ~ ., method = "rf", importance_type = "permutation",
#'   pred_wrapper = pred_wrapper_classif, metric = "RMSE", nsim = 5
#' )
#' print(res_perm$comparison)
#'
#' # SHAP importance (using fastshap)
#' # For classification tasks, this default wrapper returns the probability of the first class.
#' res_shap <- compare_feature_importance(
#'   X, Y, default ~ ., method = "rf", importance_type = "shap",
#'   nsim = 5
#' )
#' print(res_shap$comparison)
#'
#' # SHAP importance (using fastshap)
#' # For classification tasks, this default wrapper returns the probability of the first class.
#' res_shap <- compare_feature_importance(
#'   X, Y, income ~ ., method = "rf", importance_type = "shap",
#'   nsim = 5
#' )
#' print(res_shap$comparison)
#'
#' @importFrom caret train varImp trainControl
#' @importFrom data.table data.table
#' @export
compare_feature_importance <- function(X, Y, formula, method,
                                       importance_type = c("model", "permutation", "shap"),
                                       pred_wrapper = NULL, metric = "Accuracy", nsim = 5) {
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
    if ("Overall" %in% colnames(vi_X)) {
      importance_X <- vi_X[features, "Overall"]
    } else {
      importance_X <- vi_X[features, 1]
    }
    if ("Overall" %in% colnames(vi_Y)) {
      importance_Y <- vi_Y[features, "Overall"]
    } else {
      importance_Y <- vi_Y[features, 1]
    }
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

  comparison <- data.table::data.table(
    feature = features,
    importance_X = importance_X[features],
    dataset_X = "X",
    importance_Y = importance_Y[features],
    dataset_Y = "Y",
    difference = abs(importance_X[features] - importance_Y[features])
  )

  return(list(comparison = comparison, model_X = model_X, model_Y = model_Y))
}

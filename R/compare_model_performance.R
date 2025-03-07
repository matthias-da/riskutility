#' Compare Predictive Model Performance between Two Datasets
#'
#' This function evaluates and compares predictive model performance between an original dataset (X) and a synthetic/anonymized dataset (Y). The model is trained separately within each dataset using cross-validation, and various performance metrics are computed.
#'
#' @param X Original data.frame or data.table.
#' @param Y Synthetic or anonymized data.frame or data.table.
#' @param formula Formula specifying the model structure.
#' @param method A string specifying the model type (e.g., 'rf' for random forest, 'rpart' for CART). See `caret` package documentation.
#' @param metric Metric for model tuning (default 'Accuracy' for classification, 'RMSE' for regression).
#' @param trControl trainControl object from `caret` specifying cross-validation strategy (default 10-fold CV).
#'
#' @return List containing trained models, performance metrics for X and Y, and comparison results.
#'
#' @import caret
#' @importFrom data.table as.data.table
#' @export
#'
#' @examples
#' library(caret)
#' set.seed(123)
#' X <- data.frame(income = rnorm(500, 50000, 10000),
#'                 age = rnorm(500, 40, 10),
#'                 gender = factor(sample(c("M", "F"), 500, replace = TRUE)),
#'                 target = factor(sample(c("Yes", "No"), 500, replace = TRUE)))
#'
#' Y <- data.frame(income = rnorm(500, 48000, 12000),
#'                 age = rnorm(500, 42, 11),
#'                 gender = factor(sample(c("M", "F"), 500, replace = TRUE)),
#'                 target = factor(sample(c("Yes", "No"), 500, replace = TRUE)))
#'
#' result <- compare_model_performance(X, Y,
#'            formula = target ~ income + age + gender,
#'            method = 'rpart',
#'            metric = 'Accuracy')
#'
#' print(result$comparison)

compare_model_performance <- function(X, Y, formula, method = 'rpart',
                                      metric = ifelse(is.factor(X[[as.character(formula[[2]])]]), "Accuracy", "RMSE"),
                                      trControl = trainControl(method = "cv", number = 10)) {

  library(caret)

  # Train and evaluate on X using cross-validation
  model_X <- train(formula, data = X, method = method,
                   trControl = trControl, metric = metric)
  perf_X <- max(model_X$results[[metric]])

  # Train and evaluate on Y using cross-validation
  model_Y <- train(formula, data = Y, method = method,
                   trControl = trControl, metric = metric)
  perf_Y <- max(model_Y$results[[metric]])

  # Comparison
  comparison <- data.frame(
    Dataset = c("Original (X)", "Synthetic (Y)"),
    Metric = metric,
    Performance = c(perf_X, perf_Y)
  )

  # Return results
  return(list(
    model_X = model_X,
    performance_X = perf_X,
    model_Y = model_Y,
    performance_Y = perf_Y,
    comparison = comparison
  ))
}

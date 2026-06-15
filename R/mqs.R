#' Model Quality Score
#'
#' Computes the prediction quality or accuracy of ML models
#'
#' @details
#' It computes the prediction quality (RMSE) when the response is numeric,
#' and the accuracy when the response is a factor variable.
#' The computation is done for both, X (typically a non-anonymized data set)
#' and Y (typically the anonymized or synthetisized version of X) and compares
#' the estimates.
#'
#' The ratio compares these estimates: a value close to 1 indicates comparable
#' prediction quality, a ratio above 1 indicates that the synthetic data have
#' lower prediction quality (they preserve less predictive structure) than the
#' original data, and a ratio below 1 indicates that the synthetic data are at
#' least as predictable as the original.
#'
#' @param ... additional arguments passed to methods
#' @return An object of class "mqs" containing:
#' \itemize{
#'   \item mqs_ratio: the model quality statistics ratio
#'   \item mqs_table: data frame with model performance comparison
#' }
#' @seealso \code{\link{propscore}}, \code{\link{compare_model_performance}}
#'
#' @author Matthias Templ
#' @param X data frame
#' @param Y data frame
#' @param form model formula
#' @param methods classification or regression methods. In principle, all
#' methods supported from R package caret can be used.
#' @param na missing value treatment (remove, impute or stop)
#' @param ntop number of top models considered for the mqs statistics
#' @importFrom utils capture.output
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme
#' @importFrom graphics abline
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom stats glm
#' @importFrom stats binomial
#' @importFrom stats predict
#' @family utility
#' @export
#' @examples
#' \donttest{
#' # Requires the caret, caretEnsemble and rpart packages
#' if (requireNamespace("caret", quietly = TRUE) &&
#'     requireNamespace("caretEnsemble", quietly = TRUE) &&
#'     requireNamespace("rpart", quietly = TRUE)) {
#'   set.seed(123)
#'   X <- data.frame(
#'     y = factor(sample(c("A", "B"), 100, replace = TRUE)),
#'     x1 = rnorm(100),
#'     x2 = rnorm(100)
#'   )
#'   Y <- data.frame(
#'     y = factor(sample(c("A", "B"), 100, replace = TRUE)),
#'     x1 = rnorm(100, 0.1, 1),
#'     x2 = rnorm(100, 0.1, 1)
#'   )
#'   m <- mqs(X, Y, form = y ~ x1 + x2, methods = c("glm", "rpart"))
#'   m
#' }
#' }
mqs <- function(X, ...) {
  UseMethod("mqs")
}

#' @rdname mqs
#' @export
mqs.synth_pair <- function(X, ...) {
  mqs.default(X = X$original, Y = X$synthetic, ...)
}

#' @rdname mqs
#' @export
mqs.default <- function(X, Y, form,
                methods = c("glm", "knn", "simpls", "rpart", "ranger"),
                na = "remove", ntop = length(methods), ...){
  # suggestion of models:
  # categorical: glm, xgbTree, rpart, ranger, knn, naive_bayes, simpls, Linda
  # continuous response: glm, rpart, xgbTree, ranger, lars, knn, simpls

  # Check required packages
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required for mqs(). Please install it.")
  }
  if (!requireNamespace("caretEnsemble", quietly = TRUE)) {
    stop("Package 'caretEnsemble' is required for mqs(). Please install it.")
  }

  # Check if X and Y are data frames
  if (!is.data.frame(X)) {
    stop("X must be a data frame.")
  }
  if (!is.data.frame(Y)) {
    stop("Y must be a data frame.")
  }

  if (is.null(form)) {
    # If no formula is provided, check if X and Y have the same structure
    if (!all(names(X) == names(Y))) {
      stop("X and Y must have the same structure.")
    }
    if (!all(sapply(X, class) == sapply(Y, class))) {
      stop("X and Y must have the same column classes.")
    }
    form <- paste(colnames(X)[1], paste(colnames(X)[2:ncol(X)],
                                        collapse = " + "), collapse = " ~ ")
  } else {
    # Check if form is of class 'formula'
    if (!inherits(form, "formula")) {
      stop("form must be a valid formula.")
    }

    # Extract variable names from the formula
    formula_vars <- all.vars(form)[-1]

    # Check if variables in the formula appear in X (and Y)
    missing_vars <- setdiff(formula_vars, names(X))
    if (length(missing_vars) > 0) {
      stop(paste("The following variables from the formula are missing
                 in X (and Y):", paste(missing_vars, collapse = ", ")))
    }

    # Check if the columns in X and Y for the formula variables
    # have the same structure
    for (var in formula_vars) {
      if (!identical(class(X[[var]]), class(Y[[var]]))) {
        stop(paste("The variable", var, "must have the same class in
                   both X and Y."))
      }
    }
  }

  ## numeric response
  response <- X[, all.vars(form)[1]]
  if(is.numeric(response) & is.null(methods)){
    methods <- c("glm", "rpart", "ranger", "lasso")
  } else if((is.factor(response) | is.character((response))) &
            is.null(methods)){
    methods <- c("glm", "rpart", "ranger", "lasso")
  }

  # # Check the cluster argument
  # if (!is.null(cluster)) {
  #   if (!is.character(cluster) || length(cluster) != 1) {
  #     stop("cluster must be a single string representing a
  #           variable name in X.")
  #   }
  #   if (!(cluster %in% names(X))) {
  #     stop(paste("The cluster variable", cluster, "is not found in X."))
  #   }
  # }

  # Check the na argument
  valid_na_values <- c("impute", "remove", "stop")
  if (!(na %in% valid_na_values)) {
    stop(paste("na must be one of", paste(valid_na_values,
                                          collapse = ", "), "."))
  }

  vars <- attributes(terms.formula(form))$term.labels
  if(any(is.na(X[, vars]))){
    if(na == "stop") stop("data contains missing values")
    if(na == "impute"){
      X[, vars] <- VIM::kNN(X[, vars], imp_var = FALSE)
      Y[, vars] <- VIM::kNN(Y[, vars], imp_var = FALSE)
    }
    if(na == "remove"){
      X <- X[complete.cases(X[, vars]), ]
      Y <- Y[complete.cases(Y[, vars]), ]
    }
  }

  # control for caret
  my_control <- caret::trainControl(
    method="cv",
    number=10,
    savePredictions="final",
    classProbs=TRUE
  )

    # Fit the base learners on the original (X) and synthetic (Y) data.
    model_list_X <- caretEnsemble::caretList(
      form,
      data = X,
      trControl = my_control,
      methodList = methods
    )
    model_list_Y <- caretEnsemble::caretList(
      form,
      data = Y,
      trControl = my_control,
      methodList = methods
    )

    # Per-model cross-validated performance. We read the performance directly
    # from each fitted caret model via getTrainPerf() rather than from
    # caretEnsemble's summary(), whose structure is version-dependent.
    metric <- if (is.numeric(response)) "RMSE" else "Accuracy"
    get_perf <- function(model_list) {
      perf <- vapply(model_list, function(m) {
        tp <- caret::getTrainPerf(m)
        col <- paste0("Train", metric)
        if (col %in% names(tp)) tp[[col]] else NA_real_
      }, numeric(1))
      names(perf) <- names(model_list)
      perf
    }
    perf_X <- get_perf(model_list_X)
    perf_Y <- get_perf(model_list_Y)

    if (is.numeric(response)) {
      # Lower RMSE is better; ratio is computed as Y / X (as in the original).
      mqs_stat <- mean(perf_Y / perf_X)
    } else {
      # Higher accuracy is better; ratio is computed as X / Y (as in the original).
      mqs_stat <- mean(perf_X / perf_Y)
    }
    mqs_table <- rbind(perf_X, perf_Y)
    colnames(mqs_table) <- names(perf_X)
    mqs_table <- data.frame(mqs_table, check.names = FALSE)
    mqs_table <- cbind("data" = c("X", "Y"), mqs_table)
    mqs_table <- cbind(mqs_table, "measure" = c(metric, metric))
    rownames(mqs_table) <- NULL

  result <- list(
    mqs_ratio = mqs_stat,
    mqs_table = mqs_table
  )
  class(result) <- "mqs"
  return(result)

}

#' Print method for mqs objects
#'
#' @param x an object of class "mqs"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.mqs <- function(x, ...) {
  cat("Model Quality Score (MQS)\n")
  cat("=========================\n")
  cat("MQS ratio:", round(x$mqs_ratio, 4), "\n\n")
  cat("Model performance table:\n")
  print(x$mqs_table, row.names = FALSE)
  invisible(x)
}

#' Summary method for mqs objects
#'
#' @param object an object of class "mqs"
#' @param ... additional arguments (ignored)
#' @return An object of class "summary.mqs"
#' @export
summary.mqs <- function(object, ...) {
  summ <- list(
    mqs_ratio = object$mqs_ratio,
    mqs_table = object$mqs_table,
    # mqs_ratio > 1 means the original is more predictable than the synthetic
    # data (Acc_X/Acc_Y for accuracy, RMSE_Y/RMSE_X for error), i.e. the
    # synthetic data have *lower* quality; a ratio < 1 means the reverse.
    interpretation = if (object$mqs_ratio > 1.05) {
      "Synthetic data has lower prediction quality than original"
    } else if (object$mqs_ratio < 0.95) {
      "Synthetic data has better prediction quality than original"
    } else {
      "Synthetic data has comparable prediction quality to original"
    }
  )
  class(summ) <- "summary.mqs"
  summ
}

#' Print method for summary.mqs objects
#'
#' @param x an object of class "summary.mqs"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.mqs <- function(x, ...) {
  cat("Summary: Model Quality Score (MQS)\n")
  cat("===================================\n")
  cat("MQS ratio:", round(x$mqs_ratio, 4), "\n")
  cat("Interpretation:", x$interpretation, "\n\n")
  cat("Model performance table:\n")
  print(x$mqs_table, row.names = FALSE)
  invisible(x)
}


#' Plot method for mqs objects
#'
#' Visualizes the model performance comparison between original (X) and
#' synthetic (Y) data.
#'
#' @param x an object of class \code{"mqs"}
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = grouped bar chart of per-model
#'   performance, 2 = MQS ratio dot plot
#' @importFrom graphics barplot abline legend axis par mtext text
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.mqs <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  if (show[1]) {
    # Grouped bar chart: per-model performance for X and Y
    tbl <- x$mqs_table
    measure <- as.character(tbl$measure[1])
    method_cols <- setdiff(names(tbl), c("data", "measure"))

    mat <- as.matrix(tbl[, method_cols, drop = FALSE])
    rownames(mat) <- c("Original (X)", "Synthetic (Y)")

    bp <- barplot(mat,
                  beside = TRUE,
                  col = c("steelblue", "coral"),
                  ylab = measure,
                  main = paste("Model Performance Comparison:", measure),
                  las = 1,
                  ...)

    legend("topright",
           legend = c("Original (X)", "Synthetic (Y)"),
           fill = c("steelblue", "coral"),
           bty = "n",
           cex = 0.9)

    # Annotate bars
    for (i in seq_along(mat)) {
      text(bp[i], mat[i] + max(mat) * 0.02,
           labels = sprintf("%.3f", mat[i]), cex = 0.7)
    }
  }

  if (show[2]) {
    # MQS ratio dot plot
    ratio <- x$mqs_ratio
    xlim <- c(min(0.5, ratio - 0.1), max(1.5, ratio + 0.1))

    plot(ratio, 1, pch = 19, cex = 2, col = "steelblue",
         xlim = xlim, ylim = c(0.5, 1.5),
         xlab = "MQS Ratio", ylab = "",
         main = "Model Quality Score Ratio",
         yaxt = "n", ...)

    abline(v = 1, lty = 2, col = "grey50", lwd = 1.5)
    abline(v = 0.95, lty = 3, col = "coral", lwd = 1)
    abline(v = 1.05, lty = 3, col = "coral", lwd = 1)

    text(ratio, 1.15, labels = sprintf("%.4f", ratio), cex = 1.1)

    legend("topright",
           legend = c("MQS ratio", "Reference (1.0)", "Tolerance band"),
           col = c("steelblue", "grey50", "coral"),
           lty = c(NA, 2, 3),
           pch = c(19, NA, NA),
           bty = "n",
           cex = 0.8)
  }

  invisible(x)
}


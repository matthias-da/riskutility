#' Compare Missing Value Patterns Between Original and Synthetic Data
#'
#' This function compares missing value patterns in an original dataset (\code{X})
#' and a synthetic/anonymized dataset (\code{Y}). It assesses whether the synthetic data
#' mimics the same percentage and distribution of missingness across variables.
#'
#' Supported methods include:
#' \itemize{
#'   \item \code{"percentage"}: Computes and compares missing value percentages for each variable.
#'   \item \code{"pattern"}: Summarizes the frequency of missingness patterns across rows.
#' }
#'
#' Both \code{X} and \code{Y} should be numeric data frames or matrices with the same structure.
#'
#' @param X A data frame or matrix of original data.
#' @param Y A data frame or matrix of synthetic/anonymized data with the same structure as \code{X}.
#' @param method Character. The method to use for comparison: \code{"percentage"} (default) or \code{"pattern"}.
#'
#' @param ... additional arguments passed to methods
#'
#' @return An object of class \code{"missingCompare"}:
#' \itemize{
#'   \item For \code{"percentage"}: a list with a data frame (\code{summary}) reporting per-variable missing
#'   counts and percentages.
#'   \item For \code{"pattern"}: a list with two tables (\code{patterns_X} and \code{patterns_Y}) reporting the
#'   frequency of missingness patterns.
#' }
#'
#' @examples
#' \dontrun{
#'   set.seed(123)
#'   X <- data.frame(a = c(rnorm(95), rep(NA, 5)), b = c(rnorm(90), rep(NA, 10)))
#'   Y <- data.frame(a = c(rnorm(95), rep(NA, 5)), b = c(rnorm(90), rep(NA, 10)))
#'
#'   # Compare missing value percentages
#'   res_pct <- compare_missing_values(X, Y, method = "percentage")
#'   print(res_pct)
#'   summary(res_pct)
#'   plot(res_pct)
#'
#'   # Compare missing value patterns
#'   res_pat <- compare_missing_values(X, Y, method = "pattern")
#'   print(res_pat)
#'   summary(res_pat)
#'   plot(res_pat)
#' }
#'
#' @importFrom stats na.omit quantile
#' @family comparison
#' @author Matthias Templ
#' @export
compare_missing_values <- function(X, ...) {
  UseMethod("compare_missing_values")
}

#' @rdname compare_missing_values
#' @export
compare_missing_values.synth_pair <- function(X, ...) {
  compare_missing_values.default(X = X$original, Y = X$synthetic, ...)
}

#' @rdname compare_missing_values
#' @export
compare_missing_values.default <- function(X, Y, method = "percentage", ...) {
  # Check inputs
  if (!is.data.frame(X) && !is.matrix(X))
    stop("X must be a data frame or matrix.")
  if (!is.data.frame(Y) && !is.matrix(Y))
    stop("Y must be a data frame or matrix.")
  if (ncol(X) != ncol(Y))
    stop("X and Y must have the same number of columns.")

  method <- match.arg(tolower(method), choices = c("percentage", "pattern"))
  res <- list(method = method)

  if (method == "percentage") {
    summary_df <- data.frame(
      Variable = colnames(X),
      Missing_X = colSums(is.na(X)),
      Missing_Y = colSums(is.na(Y)),
      Total_X = nrow(X),
      Total_Y = nrow(Y),
      Percent_X = colSums(is.na(X)) / nrow(X),
      Percent_Y = colSums(is.na(Y)) / nrow(Y)
    )
    res$summary <- summary_df
  } else if (method == "pattern") {
    pattern_func <- function(df) {
      apply(is.na(df), 1, function(x) paste(ifelse(x, "M", "O"), collapse = ""))
    }
    res$patterns_X <- table(pattern_func(X))
    res$patterns_Y <- table(pattern_func(Y))
  }

  class(res) <- "missingCompare"
  return(res)
}

# S3 print method for missingCompare objects
#' @export
print.missingCompare <- function(x, ...) {
  cat("Missing Value Comparison Result\n")
  cat("Method: ", x$method, "\n")
  if (x$method == "percentage") {
    cat("Per-variable missing value summary:\n")
    print(x$summary)
  } else if (x$method == "pattern") {
    cat("Missing value patterns in Original Data (X):\n")
    print(x$patterns_X)
    cat("\nMissing value patterns in Synthetic Data (Y):\n")
    print(x$patterns_Y)
  }
  invisible(x)
}

# S3 summary method for missingCompare objects
#' @export
summary.missingCompare <- function(object, ...) {
  if (object$method == "percentage") {
    diff <- object$summary$Percent_Y - object$summary$Percent_X
    summ <- cbind(object$summary, Difference = diff)
    object$summary <- summ
  }
  class(object) <- c("summary.missingCompare", class(object))
  return(object)
}

# S3 print method for summary.missingCompare objects
#' @export
print.summary.missingCompare <- function(x, ...) {
  cat("Summary of Missing Value Comparison\n")
  if (x$method == "percentage") {
    print(x$summary)
  } else if (x$method == "pattern") {
    cat("Missing value patterns in Original Data (X):\n")
    print(x$patterns_X)
    cat("\nMissing value patterns in Synthetic Data (Y):\n")
    print(x$patterns_Y)
  }
  invisible(x)
}

# S3 plot method for missingCompare objects
#' @export
plot.missingCompare <- function(x, ...) {
  if (x$method == "percentage") {
    m <- rbind(x$summary$Percent_X, x$summary$Percent_Y)
    colnames(m) <- x$summary$Variable
    barplot(m, beside = TRUE, col = c("skyblue", "orange"),
            legend = c("Original (X)", "Synthetic (Y)"),
            main = "Missing Value Percentage Comparison",
            ylab = "Proportion Missing", ...)
  } else if (x$method == "pattern") {
    op <- par(mfrow = c(1,2))
    barplot(x$patterns_X, main = "Missing Patterns in X", col = "skyblue", ...)
    barplot(x$patterns_Y, main = "Missing Patterns in Y", col = "orange", ...)
    par(op)
  }
}

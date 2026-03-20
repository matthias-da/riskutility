#' Compare Outlier Detection Between Original and Synthetic Data
#'
#' This function compares the prevalence of outliers in an original dataset (\code{X})
#' and an anonymized/synthetic dataset (\code{Y}) using various outlier detection methods.
#' Supported methods include:
#' \itemize{
#'   \item \code{"zscore"}: Uses z-scores (default threshold = 3) to flag outliers.
#'   \item \code{"iqr"}: Flags values outside 1.5 times the interquartile range.
#'   \item \code{"dbscan"}: Uses DBSCAN clustering to identify noise points as outliers.
#'   \item \code{"robust"}: Uses robust Mahalanobis distances based on the Minimum Covariance Determinant.
#' }
#'
#' Input datasets must be numeric and of the same structure. For multivariate methods, the
#' comparison is performed on the entire dataset.
#'
#' @param X A data frame or matrix of original data with numeric columns.
#' @param Y A data frame or matrix of anonymized/synthetic data with the same structure as \code{X}.
#' @param method Character. Outlier detection method: \code{"zscore"}, \code{"iqr"}, \code{"dbscan"}, or \code{"robust"}.
#' @param threshold Numeric. For \code{"zscore"}, the absolute z-score cutoff (default = 3). For \code{"robust"},
#' the significance level for the chi-square cutoff (default = 0.975). Ignored for \code{"iqr"} and \code{"dbscan"}.
#' @param ... additional arguments passed to the underlying method (e.g., \code{eps} and \code{minPts} for DBSCAN).
#' @return A list with:
#' \itemize{
#'   \item \code{summary}: A data frame summarizing outlier counts and proportions in \code{X} and \code{Y}.
#'   \item \code{details}: Method-specific details (e.g., per-variable counts or clustering output).
#' }
#'
#' @examples
#' \dontrun{
#'   set.seed(123)
#'   X <- data.frame(a = rnorm(100), b = rnorm(100))
#'   Y <- data.frame(a = rnorm(100, mean = 0.1), b = rnorm(100, mean = -0.1))
#'
#'   # Using z-score method
#'   res_z <- compare_outliers(X, Y, method = "zscore", threshold = 3)
#'
#'   # Using IQR method
#'   res_iqr <- compare_outliers(X, Y, method = "iqr")
#'
#'   # Using DBSCAN (requires 'dbscan' package)
#'   res_db <- compare_outliers(X, Y, method = "dbscan", eps = 0.5, minPts = 5)
#'
#'   # Using robust Mahalanobis (requires 'robustbase' package)
#'   res_robust <- compare_outliers(X, Y, method = "robust", threshold = 0.975)
#' }
#'
#' @family comparison
#' @author Matthias Templ
#' @export
compare_outliers <- function(X, ...) {
  UseMethod("compare_outliers")
}

#' @rdname compare_outliers
#' @export
compare_outliers.synth_pair <- function(X, ...) {
  compare_outliers.default(X = X$original, Y = X$synthetic, ...)
}

#' @rdname compare_outliers
#' @export
compare_outliers.default <- function(X, Y, method = "zscore", threshold = NULL, ...) {
  # Check inputs
  if (!is.data.frame(X) && !is.matrix(X))
    stop("X must be a data frame or matrix.")
  if (!is.data.frame(Y) && !is.matrix(Y))
    stop("Y must be a data frame or matrix.")
  if (ncol(X) != ncol(Y))
    stop("X and Y must have the same number of columns.")
  if (!all(sapply(X, is.numeric)) || !all(sapply(Y, is.numeric)))
    stop("All columns in X and Y must be numeric.")

  method <- match.arg(tolower(method), choices = c("zscore", "iqr", "dbscan", "robust"))

  summary_df <- data.frame(Variable = character(0),
                           Outliers_X = integer(0),
                           Outliers_Y = integer(0),
                           Total_X = integer(0),
                           Total_Y = integer(0),
                           Prop_X = numeric(0),
                           Prop_Y = numeric(0))
  details <- list()

  if (method %in% c("zscore", "iqr")) {
    for (col in colnames(X)) {
      vecX <- X[[col]]
      vecY <- Y[[col]]

      if (method == "zscore") {
        thres <- ifelse(is.null(threshold), 3, threshold)
        zX <- abs((vecX - mean(vecX, na.rm = TRUE)) / sd(vecX, na.rm = TRUE))
        zY <- abs((vecY - mean(vecY, na.rm = TRUE)) / sd(vecY, na.rm = TRUE))
        outX <- sum(zX > thres, na.rm = TRUE)
        outY <- sum(zY > thres, na.rm = TRUE)
      } else if (method == "iqr") {
        Q1_X <- quantile(vecX, 0.25, na.rm = TRUE)
        Q3_X <- quantile(vecX, 0.75, na.rm = TRUE)
        IQR_X <- Q3_X - Q1_X
        lower_X <- Q1_X - 1.5 * IQR_X
        upper_X <- Q3_X + 1.5 * IQR_X
        outX <- sum(vecX < lower_X | vecX > upper_X, na.rm = TRUE)

        Q1_Y <- quantile(vecY, 0.25, na.rm = TRUE)
        Q3_Y <- quantile(vecY, 0.75, na.rm = TRUE)
        IQR_Y <- Q3_Y - Q1_Y
        lower_Y <- Q1_Y - 1.5 * IQR_Y
        upper_Y <- Q3_Y + 1.5 * IQR_Y
        outY <- sum(vecY < lower_Y | vecY > upper_Y, na.rm = TRUE)
      }

      totalX <- length(vecX)
      totalY <- length(vecY)
      summary_df <- rbind(summary_df, data.frame(Variable = col,
                                                 Outliers_X = outX,
                                                 Outliers_Y = outY,
                                                 Total_X = totalX,
                                                 Total_Y = totalY,
                                                 Prop_X = outX / totalX,
                                                 Prop_Y = outY / totalY))
    }
    details <- summary_df
  } else if (method == "dbscan") {
    if (!requireNamespace("dbscan", quietly = TRUE))
      stop("Package 'dbscan' is required for the 'dbscan' method. Please install it.")
    db_X <- dbscan::dbscan(X, ...)
    db_Y <- dbscan::dbscan(Y, ...)
    outX <- sum(db_X$cluster == 0)
    outY <- sum(db_Y$cluster == 0)
    summary_df <- data.frame(Variable = "Multivariate",
                             Outliers_X = outX,
                             Outliers_Y = outY,
                             Total_X = nrow(X),
                             Total_Y = nrow(Y),
                             Prop_X = outX / nrow(X),
                             Prop_Y = outY / nrow(Y))
    details <- list(dbscan_X = db_X, dbscan_Y = db_Y)
  } else if (method == "robust") {
    if (!requireNamespace("robustbase", quietly = TRUE))
      stop("Package 'robustbase' is required for the 'robust' method. Please install it.")
    p <- ncol(X)
    chi_thresh <- ifelse(is.null(threshold), qchisq(0.975, df = p), qchisq(threshold, df = p))
    covX <- robustbase::covMcd(X)
    covY <- robustbase::covMcd(Y)
    mdX <- mahalanobis(X, center = covX$center, cov = covX$cov)
    mdY <- mahalanobis(Y, center = covY$center, cov = covY$cov)
    outX <- sum(mdX > chi_thresh)
    outY <- sum(mdY > chi_thresh)
    summary_df <- data.frame(variable = "Multivariate",
                             outliers_X = outX,
                             outliers_Y = outY,
                             total_X = nrow(X),
                             total_Y = nrow(Y),
                             prop_X = outX / nrow(X),
                             prop_Y = outY / nrow(Y))
    details <- list(mahalanobis_X = mdX, mahalanobis_Y = mdY, chi_threshold = chi_thresh, mult_outliers_X = outX, mult_outliers_Y = outY)
  }

  return(list(summary = summary_df, details = details))
}

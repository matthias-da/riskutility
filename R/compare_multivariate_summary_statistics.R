#' Compare Multivariate Summary Statistics between Two Datasets
#'
#' This function computes higher-order summary statistics for multiple variables in two datasets (e.g., an original dataset and an anonymized/synthetic dataset) to assess whether the synthetic data replicate the joint properties of the original data. For continuous variables, the function calculates the multivariate mean vector, variance–covariance matrix, and correlation matrix. When sampling weights are provided, weighted estimates are computed. For categorical variables, joint frequency tables (and relative frequencies) are constructed. The function returns the statistics for each dataset along with the differences between them.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the anonymized/synthetic dataset.
#' @param cont_vars A character vector specifying the continuous variables to compare.
#' @param cat_vars A character vector specifying the categorical variables for joint frequency comparison.
#' @param weight_X Optional. A character string specifying the sampling weight variable in X.
#' @param weight_Y Optional. A character string specifying the sampling weight variable in Y.
#'
#' @param ... additional arguments passed to methods
#'
#' @return A list with two elements:
#' \describe{
#'   \item{continuous}{A list containing:
#'      \describe{
#'         \item{mean_X}{The multivariate mean vector for X.}
#'         \item{mean_Y}{The multivariate mean vector for Y.}
#'         \item{mean_diff}{The difference between the mean vectors (X minus Y).}
#'         \item{cov_X}{The variance–covariance matrix for X.}
#'         \item{cov_Y}{The variance–covariance matrix for Y.}
#'         \item{cov_diff}{The difference between the covariance matrices (X minus Y).}
#'         \item{cor_X}{The correlation matrix for X.}
#'         \item{cor_Y}{The correlation matrix for Y.}
#'         \item{cor_diff}{The absolute difference between the correlation matrices.}
#'      }
#'   }
#'   \item{categorical}{A list containing:
#'      \describe{
#'         \item{freq_X}{The joint frequency table for the categorical variables in X.}
#'         \item{freq_Y}{The joint frequency table for the categorical variables in Y.}
#'         \item{rel_freq_X}{The relative frequency table for X.}
#'         \item{rel_freq_Y}{The relative frequency table for Y.}
#'      }
#'   }
#' }
#'
#' @details For continuous variables, if sampling weights are provided the function computes weighted means using \code{weighted.mean} and weighted covariance matrices using a custom function. The correlation matrices are derived from the covariance matrices. For categorical variables, if weights are provided the frequencies are computed as the sum of weights; otherwise, raw counts are used. Joint frequency tables are normalized to obtain relative frequencies.
#'
#' @importFrom data.table as.data.table rbindlist
#' @family comparison
#' @export
#'
#' @examples
#' set.seed(123)
#' X <- data.frame(
#'   income = rnorm(500, mean = 50000, sd = 10000),
#'   age = rnorm(500, mean = 40, sd = 10),
#'   gender = factor(sample(c("Male", "Female"), 500, replace = TRUE)),
#'   region = factor(sample(c("North", "South", "East", "West"), 500, replace = TRUE)),
#'   weight = runif(500, 0.5, 1.5)
#' )
#'
#' Y <- data.frame(
#'   income = rnorm(1000, mean = 48000, sd = 12000),
#'   age = rnorm(1000, mean = 42, sd = 11),
#'   gender = factor(sample(c("Male", "Female"), 1000, replace = TRUE)),
#'   region = factor(sample(c("North", "South", "East", "West"), 1000, replace = TRUE)),
#'   weight = runif(1000, 0.5, 1.5)
#' )
#'
#' result <- compare_multivariate_summary_statistics(X, Y,
#'               cont_vars = c("income", "age"),
#'               cat_vars = c("gender", "region"),
#'               weight_X = "weight", weight_Y = "weight")
#'
#' print(result$continuous)
#' print(result$categorical)
compare_multivariate_summary_statistics <- function(X, ...) {
  UseMethod("compare_multivariate_summary_statistics")
}

#' @rdname compare_multivariate_summary_statistics
#' @export
compare_multivariate_summary_statistics.synth_pair <- function(X, ...) {
  compare_multivariate_summary_statistics.default(X = X$original, Y = X$synthetic, ...)
}

#' @rdname compare_multivariate_summary_statistics
#' @export
compare_multivariate_summary_statistics.default <- function(X, Y, cont_vars, cat_vars,
                                                    weight_X = NULL, weight_Y = NULL, ...) {
  # Convert datasets to data.table if needed
  X <- as.data.table(X)
  Y <- as.data.table(Y)

  ## Continuous variables: Compute multivariate summary statistics
  # Helper function for weighted covariance matrix
  weighted_cov <- function(mat, w) {
    # Compute weighted mean for centering
    wm <- sapply(as.data.frame(mat), weighted.mean, w = w, na.rm = TRUE)
    centered <- sweep(mat, 2, wm)
    cov_mat <- crossprod(centered, centered * w) / (sum(w) - 1)
    return(cov_mat)
  }

  # For dataset X
  if (!is.null(weight_X) && weight_X %in% names(X)) {
    wX <- X[[weight_X]]
  } else {
    wX <- rep(1, nrow(X))
  }
  X_cont <- as.matrix(X[, ..cont_vars])
  mean_X <- sapply(as.data.frame(X_cont), function(x) weighted.mean(x, wX, na.rm = TRUE))
  cov_X <- weighted_cov(X_cont, wX)
  cor_X <- cov2cor(cov_X)

  # For dataset Y
  if (!is.null(weight_Y) && weight_Y %in% names(Y)) {
    wY <- Y[[weight_Y]]
  } else {
    wY <- rep(1, nrow(Y))
  }
  Y_cont <- as.matrix(Y[, ..cont_vars])
  mean_Y <- sapply(as.data.frame(Y_cont), function(x) weighted.mean(x, wY, na.rm = TRUE))
  cov_Y <- weighted_cov(Y_cont, wY)
  cor_Y <- cov2cor(cov_Y)

  continuous <- list(
    mean_X = mean_X,
    mean_Y = mean_Y,
    mean_diff = mean_X - mean_Y,
    cov_X = cov_X,
    cov_Y = cov_Y,
    cov_diff = cov_X - cov_Y,
    cor_X = cor_X,
    cor_Y = cor_Y,
    cor_diff = abs(cor_X - cor_Y)
  )

  ## Categorical variables: Compute joint frequency tables
  # For dataset X
  if (!is.null(weight_X) && weight_X %in% names(X)) {
    freq_X <- X[, .(freq = sum(get(weight_X), na.rm = TRUE)), by = cat_vars]
  } else {
    freq_X <- X[, .N, by = cat_vars]
    setnames(freq_X, "N", "freq")
  }
  freq_X[, rel_freq := freq / sum(freq)]

  # For dataset Y
  if (!is.null(weight_Y) && weight_Y %in% names(Y)) {
    freq_Y <- Y[, .(freq = sum(get(weight_Y), na.rm = TRUE)), by = cat_vars]
  } else {
    freq_Y <- Y[, .N, by = cat_vars]
    setnames(freq_Y, "N", "freq")
  }
  freq_Y[, rel_freq := freq / sum(freq)]

  categorical <- list(
    freq_X = freq_X,
    freq_Y = freq_Y
  )

  return(list(continuous = continuous,
              categorical = categorical))
}

#' Compare Multivariate Distributions using Mahalanobis Distance or Jensen-Shannon Divergence
#'
#' This function compares the joint distribution of a set of variables between two datasets (X and Y)
#' using one of two methods. For continuous variables, it computes the weighted (if weights are provided)
#' Mahalanobis distance between the mean vectors, using a pooled covariance matrix (or its generalized inverse).
#' For nominal variables, it discretizes the data (if necessary) and computes a divergence based on the Jensen-Shannon
#' divergence of the joint frequency distributions, which serves as a symmetric measure of distributional difference.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the anonymized/synthetic dataset.
#' @param vars A character vector specifying the variables to compare.
#' @param method A character string indicating which method to use. Either "mahalanobis" (default) or "mutual_information".
#' @param weight_X Optional. A character string specifying the sampling weight variable in X (used only for method "mahalanobis").
#' @param weight_Y Optional. A character string specifying the sampling weight variable in Y (used only for method "mahalanobis").
#' @param n_bins Numeric. The number of bins to use when discretizing variables for the mutual_information method. Default is 10.
#'
#' @param ... additional arguments passed to methods
#'
#' @return A list containing:
#' \describe{
#'   \item{distance}{The computed distance measure (Mahalanobis distance or Jensen-Shannon divergence).}
#'   \item{method}{The method used ("mahalanobis" or "mutual_information").}
#'   \item{additional}{A list of additional computed quantities (e.g., weighted means and pooled covariance for Mahalanobis, or probability vectors for mutual_information).}
#' }
#'
#' @details For the "mahalanobis" method, if sampling weights are provided, the function computes the weighted mean vector for each dataset
#' and a weighted pooled covariance matrix. The Mahalanobis distance is then computed as
#' \deqn{D_M = \sqrt{(\mu_X - \mu_Y)^T S^{-1} (\mu_X - \mu_Y)}}
#' where \eqn{S} is the pooled covariance matrix (or its generalized inverse if singular). For the "mutual_information" method,
#' numeric variables are discretized into \code{n_bins} bins using \code{cut()}, and joint frequency tables are constructed.
#' The tables are normalized into probability vectors, and the Jensen-Shannon divergence is computed as
#' \deqn{JS(p,q) = \frac{1}{2} KL(p||m) + \frac{1}{2} KL(q||m)}
#' where \eqn{m = \frac{1}{2}(p+q)} and \eqn{KL} is the Kullback-Leibler divergence.
#'
#' @importFrom data.table as.data.table rbindlist
#' @importFrom MASS ginv
#' @family comparison
#' @export
#'
#' @examples
#' # Continuous example using Mahalanobis distance (weighted)
#' set.seed(123)
#' X <- data.frame(
#'   income = rnorm(500, mean = 50000, sd = 10000),
#'   age = rnorm(500, mean = 40, sd = 10),
#'   weight = runif(500, 0.5, 1.5)
#' )
#' Y <- data.frame(
#'   income = rnorm(1000, mean = 48000, sd = 12000),
#'   age = rnorm(1000, mean = 42, sd = 11),
#'   weight = runif(1000, 0.5, 1.5)
#' )
#' result_mahal <- compare_multivariate_distribution(X, Y, vars = c("income", "age"),
#'                                                   method = "mahalanobis",
#'                                                   weight_X = "weight", weight_Y = "weight")
#' print(result_mahal)
#'
#' # Nominal example using mutual information (Jensen-Shannon divergence)
#' set.seed(456)
#' X_nom <- data.frame(
#'   gender = factor(sample(c("Male", "Female"), 500, replace = TRUE)),
#'   region = factor(sample(c("North", "South", "East", "West"), 500, replace = TRUE))
#' )
#' Y_nom <- data.frame(
#'   gender = factor(sample(c("Male", "Female"), 1000, replace = TRUE)),
#'   region = factor(sample(c("North", "South", "East", "West"), 1000, replace = TRUE))
#' )
#' result_js <- compare_multivariate_distribution(X_nom, Y_nom, vars = c("gender", "region"),
#'                                                method = "mutual_information", n_bins = 5)
#' print(result_js)
compare_multivariate_distribution <- function(X, ...) {
  UseMethod("compare_multivariate_distribution")
}

#' @rdname compare_multivariate_distribution
#' @export
compare_multivariate_distribution.synth_pair <- function(X, ...) {
  compare_multivariate_distribution.default(X = X$original, Y = X$synthetic, ...)
}

#' @rdname compare_multivariate_distribution
#' @export
compare_multivariate_distribution.default <- function(X, Y, vars, method = "mahalanobis",
                                              weight_X = NULL, weight_Y = NULL,
                                              n_bins = 10, ...) {
  # Convert to data.table if necessary
  X <- as.data.table(X)
  Y <- as.data.table(Y)

  # Subset only the variables of interest and weights if provided
  if (!is.null(weight_X)) {
    if (!(weight_X %in% names(X))) stop("weight_X not found in X")
    X_sub <- X[, c(vars, weight_X), with = FALSE]
  } else {
    X_sub <- X[, ..vars]
  }
  if (!is.null(weight_Y)) {
    if (!(weight_Y %in% names(Y))) stop("weight_Y not found in Y")
    Y_sub <- Y[, c(vars, weight_Y), with = FALSE]
  } else {
    Y_sub <- Y[, ..vars]
  }

  if (method == "mahalanobis") {
    # Continuous method: Compute weighted means and pooled covariance matrix
    if (!is.null(weight_X)) {
      wX <- X_sub[[weight_X]]
    } else {
      wX <- rep(1, nrow(X_sub))
    }
    if (!is.null(weight_Y)) {
      wY <- Y_sub[[weight_Y]]
    } else {
      wY <- rep(1, nrow(Y_sub))
    }

    # Convert X_sub and Y_sub to matrices (only the variables in vars)
    X_mat <- as.matrix(X_sub[, ..vars])
    Y_mat <- as.matrix(Y_sub[, ..vars])

    # Function to compute weighted mean for a vector
    weighted_mean <- function(x, w) {
      weighted.mean(x, w, na.rm = TRUE)
    }
    mu_X <- sapply(as.data.frame(X_mat), weighted_mean, w = wX)
    mu_Y <- sapply(as.data.frame(Y_mat), weighted_mean, w = wY)

    # Function to compute weighted covariance matrix for a matrix
    weighted_cov <- function(mat, w) {
      # Compute weighted mean for centering
      wm <- sapply(as.data.frame(mat), weighted.mean, w = w, na.rm = TRUE)
      centered <- sweep(mat, 2, wm)
      cov_mat <- crossprod(centered, centered * w) / (sum(w) - 1)
      return(cov_mat)
    }
    cov_X <- if (!is.null(weight_X)) weighted_cov(X_mat, wX) else cov(X_mat, use = "complete.obs")
    cov_Y <- if (!is.null(weight_Y)) weighted_cov(Y_mat, wY) else cov(Y_mat, use = "complete.obs")

    n_X <- nrow(X_mat)
    n_Y <- nrow(Y_mat)

    # Compute pooled covariance matrix
    pooled_cov <- ((n_X - 1) * cov_X + (n_Y - 1) * cov_Y) / (n_X + n_Y - 2)
    inv_cov <- tryCatch(solve(pooled_cov), error = function(e) MASS::ginv(pooled_cov))

    mahal_sq <- t(mu_X - mu_Y) %*% inv_cov %*% (mu_X - mu_Y)
    mahal_dist <- sqrt(mahal_sq)

    result <- list(distance = as.numeric(mahal_dist),
                   method = "mahalanobis",
                   additional = list(mu_X = mu_X, mu_Y = mu_Y, pooled_cov = pooled_cov))
    return(result)

  } else if (method == "mutual_information") {
    # Nominal method: Discretize numeric variables if needed.
    discretize_vars <- function(data, n_bins) {
      data_copy <- copy(data)
      for (v in names(data_copy)) {
        if (is.numeric(data_copy[[v]])) {
          data_copy[[v]] <- cut(data_copy[[v]], breaks = n_bins, include.lowest = TRUE)
        }
      }
      return(data_copy)
    }

    X_disc <- discretize_vars(X_sub[, ..vars], n_bins)
    Y_disc <- discretize_vars(Y_sub[, ..vars], n_bins)

    # Compute joint frequency tables and convert to data.table
    tab_X <- as.data.table(table(X_disc))
    tab_Y <- as.data.table(table(Y_disc))

    # Determine frequency column name (either "Freq" or "N")
    freq_col <- if ("Freq" %in% names(tab_X)) "Freq" else if ("N" %in% names(tab_X)) "N" else stop("Frequency column not found in tabulated data.")

    # Create complete grid of all combinations, force it to be a data.table
    complete_grid <- as.data.table(do.call(expand.grid, lapply(X_disc, function(x) levels(as.factor(x)))))

    # Merge complete grid with frequency tables to fill missing combinations with zeros
    complete_X <- merge(complete_grid, tab_X, by = names(complete_grid), all.x = TRUE)
    complete_Y <- merge(complete_grid, tab_Y, by = names(complete_grid), all.x = TRUE)
    complete_X[is.na(get(freq_col)), (freq_col) := 0]
    complete_Y[is.na(get(freq_col)), (freq_col) := 0]

    p_X <- complete_X[[freq_col]] / sum(complete_X[[freq_col]])
    p_Y <- complete_Y[[freq_col]] / sum(complete_Y[[freq_col]])

    m <- 0.5 * (p_X + p_Y)
    KL <- function(p, q) {
      valid <- p > 0
      sum(p[valid] * log(p[valid] / q[valid]))
    }
    js_div <- 0.5 * KL(p_X, m) + 0.5 * KL(p_Y, m)

    result <- list(distance = js_div,
                   method = "mutual_information (Jensen-Shannon divergence)",
                   additional = list(p_X = p_X, p_Y = p_Y))
    return(result)

  } else {
    stop("Invalid method. Choose 'mahalanobis' or 'mutual_information'.")
  }
}

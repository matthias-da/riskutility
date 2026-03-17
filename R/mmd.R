#' Maximum Mean Discrepancy (MMD) for Multivariate Data
#'
#' Computes the Maximum Mean Discrepancy between two datasets using kernel-based
#' two-sample testing. MMD measures whether two samples come from the same
#' distribution by comparing their embeddings in a reproducing kernel Hilbert
#' space (RKHS). It equals zero if and only if the distributions are identical.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the synthetic/anonymized dataset.
#' @param vars Character vector of numeric variable names to compare.
#'   If NULL (default), all common numeric variables are used.
#' @param kernel Character, kernel function to use. One of \code{"gaussian"}
#'   (default) or \code{"rational_quadratic"}.
#' @param method Character, computation method. One of \code{"exact"} (default,
#'   O(n^2) complexity) or \code{"rff"} (Random Fourier Features, O(nD)
#'   approximation for large datasets).
#' @param n_features Integer, number of random Fourier features when
#'   \code{method = "rff"}. Default 500.
#' @param n_perm Integer or NULL, number of permutations for a permutation test.
#'   If NULL (default), no permutation test is performed.
#' @param standardize Logical, whether to standardize variables before computing.
#'   Default TRUE (recommended for variables on different scales).
#' @param na.rm Logical, whether to remove rows with NA values. Default TRUE.
#' @param seed Integer, random seed for reproducibility. Default NULL.
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class \code{"mmd"} containing:
#' \itemize{
#'   \item \code{mmd2}: the unbiased MMD^2 statistic
#'   \item \code{sigma}: the kernel bandwidth (median heuristic)
#'   \item \code{kernel}: the kernel used
#'   \item \code{method}: the computation method used
#'   \item \code{n_X}, \code{n_Y}: sample sizes
#'   \item \code{n_vars}: number of variables
#'   \item \code{vars}: variable names used
#'   \item \code{standardized}: whether standardization was applied
#'   \item \code{utility_score}: transformed score in \eqn{[0,1]} (higher = better)
#'   \item \code{perm_pvalue}: permutation test p-value (if \code{n_perm} was set)
#'   \item \code{perm_null}: vector of permutation MMD^2 values (if \code{n_perm} was set)
#' }
#'
#' @details
#' The unbiased MMD^2 estimator with kernel \eqn{k} is:
#' \deqn{MMD^2 = \frac{1}{n(n-1)} \sum_{i \neq j} k(x_i, x_j)
#'   + \frac{1}{m(m-1)} \sum_{i \neq j} k(y_i, y_j)
#'   - \frac{2}{nm} \sum_i \sum_j k(x_i, y_j)}
#'
#' The Gaussian kernel is \eqn{k(x, y) = \exp(-\|x - y\|^2 / (2\sigma^2))}.
#' The rational quadratic kernel is \eqn{k(x, y) = (1 + \|x - y\|^2 / (2\sigma^2))^{-1}}.
#'
#' The bandwidth \eqn{\sigma} is chosen via the median heuristic: the median of
#' all pairwise Euclidean distances (subsampled to 2000 points for efficiency).
#'
#' For large datasets (\eqn{n > 5000}), consider using \code{method = "rff"}
#' which approximates the kernel via Random Fourier Features at O(nD) cost
#' instead of O(n^2).
#'
#' The utility score is computed as \eqn{\exp(-MMD^2 / \sigma^2)}, mapping to
#' \eqn{[0,1]} where 1 indicates identical distributions.
#'
#' @seealso \code{\link{energy_distance}} for energy distance,
#'   \code{\link{compare_wasserstein}} for univariate Wasserstein distance,
#'   \code{\link{gower}} for mixed-type data
#'
#' @references
#' Gretton, A., Borgwardt, K. M., Rasch, M. J., Schoelkopf, B., Smola, A. (2012).
#' A Kernel Two-Sample Test. Journal of Machine Learning Research, 13, 723-773.
#'
#' Rahimi, A. and Recht, B. (2007). Random Features for Large-Scale Kernel
#' Machines. Advances in Neural Information Processing Systems, 20.
#'
#' @author Matthias Templ
#' @family utility
#' @export
#' @importFrom stats complete.cases sd dist median rnorm runif
#'
#' @examples
#' set.seed(123)
#' # Original data
#' X <- data.frame(
#'   income = rnorm(200, mean = 50000, sd = 10000),
#'   age = rnorm(200, mean = 40, sd = 10),
#'   score = rnorm(200, mean = 100, sd = 15)
#' )
#'
#' # Good synthetic data (similar distribution)
#' Y_good <- data.frame(
#'   income = rnorm(200, mean = 50000, sd = 10000),
#'   age = rnorm(200, mean = 40, sd = 10),
#'   score = rnorm(200, mean = 100, sd = 15)
#' )
#'
#' # Poor synthetic data (shifted distribution)
#' Y_poor <- data.frame(
#'   income = rnorm(200, mean = 60000, sd = 15000),
#'   age = rnorm(200, mean = 50, sd = 15),
#'   score = rnorm(200, mean = 80, sd = 25)
#' )
#'
#' result_good <- mmd(X, Y_good, seed = 42)
#' print(result_good)
#'
#' result_poor <- mmd(X, Y_poor, seed = 42)
#' print(result_poor)
#'
#' \donttest{
#' # With permutation test
#' result_perm <- mmd(X, Y_poor, n_perm = 200, seed = 42)
#' summary(result_perm)
#'
#' # Using Random Fourier Features
#' result_rff <- mmd(X, Y_good, method = "rff", seed = 42)
#' print(result_rff)
#' }
mmd <- function(X, ...) {
  UseMethod("mmd")
}

#' @rdname mmd
#' @export
mmd.synth_pair <- function(X, ...) {
  mmd.default(
    X = X$original,
    Y = X$synthetic,
    vars = X$num_vars,
    ...
  )
}

#' @rdname mmd
#' @export
mmd.default <- function(X, Y,
                        vars = NULL,
                        kernel = c("gaussian", "rational_quadratic"),
                        method = c("exact", "rff"),
                        n_features = 500L,
                        n_perm = NULL,
                        standardize = TRUE,
                        na.rm = TRUE,
                        seed = NULL,
                        ...) {

  kernel <- match.arg(kernel)
  method <- match.arg(method)
  n_features <- as.integer(n_features)

  # Convert to data.frame if needed
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)

  # Auto-detect numeric variables if not specified
  if (is.null(vars)) {
    num_vars_X <- names(X)[sapply(X, is.numeric)]
    num_vars_Y <- names(Y)[sapply(Y, is.numeric)]
    vars <- intersect(num_vars_X, num_vars_Y)
  }

  if (length(vars) == 0) {
    stop("No numeric variables found or specified for comparison.")
  }

  # Check variables exist
  missing_X <- setdiff(vars, names(X))
  missing_Y <- setdiff(vars, names(Y))
  if (length(missing_X) > 0) {
    stop(paste("Variables missing in X:", paste(missing_X, collapse = ", ")))
  }
  if (length(missing_Y) > 0) {
    stop(paste("Variables missing in Y:", paste(missing_Y, collapse = ", ")))
  }

  # Check all specified variables are numeric
  non_num_X <- vars[!sapply(X[, vars, drop = FALSE], is.numeric)]
  non_num_Y <- vars[!sapply(Y[, vars, drop = FALSE], is.numeric)]
  if (length(non_num_X) > 0) {
    stop(paste("Non-numeric variables in X:", paste(non_num_X, collapse = ", ")))
  }
  if (length(non_num_Y) > 0) {
    stop(paste("Non-numeric variables in Y:", paste(non_num_Y, collapse = ", ")))
  }

  # Subset to selected variables
  X <- X[, vars, drop = FALSE]
  Y <- Y[, vars, drop = FALSE]

  # Handle missing values
  if (na.rm) {
    X <- X[complete.cases(X), , drop = FALSE]
    Y <- Y[complete.cases(Y), , drop = FALSE]
  }

  if (nrow(X) == 0) stop("No complete cases in X after removing NAs.")
  if (nrow(Y) == 0) stop("No complete cases in Y after removing NAs.")

  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  n_X <- nrow(X)
  n_Y <- nrow(Y)

  # Suggest RFF for large datasets
  if (method == "exact" && (n_X > 5000 || n_Y > 5000)) {
    message("Dataset has >5000 rows. Consider method='rff' for faster computation.")
  }

  # Standardize if requested (using pooled mean/sd)
  if (standardize && length(vars) > 0) {
    combined <- rbind(X, Y)
    for (v in vars) {
      m <- mean(combined[[v]], na.rm = TRUE)
      s <- sd(combined[[v]], na.rm = TRUE)
      if (s > 0) {
        X[[v]] <- (X[[v]] - m) / s
        Y[[v]] <- (Y[[v]] - m) / s
      }
    }
  }

  # Convert to matrices
  X_mat <- as.matrix(X)
  Y_mat <- as.matrix(Y)

  # Compute bandwidth via median heuristic
  sigma <- .mmd_median_heuristic(X_mat, Y_mat)

  # Compute MMD^2
  if (method == "exact") {
    mmd2 <- .mmd_exact(X_mat, Y_mat, sigma = sigma, kernel = kernel)
  } else {
    mmd2 <- .mmd_rff(X_mat, Y_mat, sigma = sigma, n_features = n_features)
  }

  # Permutation test
  perm_pvalue <- NULL
  perm_null <- NULL
  if (!is.null(n_perm) && n_perm > 0) {
    n_perm <- as.integer(n_perm)
    combined <- rbind(X_mat, Y_mat)
    n_total <- nrow(combined)
    perm_null <- numeric(n_perm)
    for (p in seq_len(n_perm)) {
      idx <- sample(n_total)
      X_perm <- combined[idx[seq_len(n_X)], , drop = FALSE]
      Y_perm <- combined[idx[(n_X + 1):n_total], , drop = FALSE]
      if (method == "exact") {
        perm_null[p] <- .mmd_exact(X_perm, Y_perm, sigma = sigma, kernel = kernel)
      } else {
        perm_null[p] <- .mmd_rff(X_perm, Y_perm, sigma = sigma, n_features = n_features)
      }
    }
    perm_pvalue <- (sum(perm_null >= mmd2) + 1) / (n_perm + 1)
  }

  # Utility score: exp(-mmd2 / sigma^2), maps to [0,1], 1 = identical
  if (sigma > 0) {
    utility_score <- exp(-mmd2 / sigma^2)
  } else {
    utility_score <- if (mmd2 == 0) 1 else 0
  }
  # Clamp to [0, 1]
  utility_score <- max(0, min(1, utility_score))

  result <- list(
    mmd2 = mmd2,
    sigma = sigma,
    kernel = kernel,
    method = method,
    n_features = if (method == "rff") n_features else NULL,
    n_X = n_X,
    n_Y = n_Y,
    n_vars = length(vars),
    vars = vars,
    standardized = standardize,
    utility_score = utility_score,
    perm_pvalue = perm_pvalue,
    perm_null = perm_null,
    n_perm = n_perm
  )

  class(result) <- "mmd"
  return(result)
}


# ---- Internal helpers --------------------------------------------------------

# Median heuristic for kernel bandwidth
# Computes sigma = median of pairwise Euclidean distances.
# Subsamples to 2000 points for efficiency.
# @param X_mat numeric matrix
# @param Y_mat numeric matrix
# @return scalar bandwidth sigma
# @keywords internal
.mmd_median_heuristic <- function(X_mat, Y_mat) {
  combined <- rbind(X_mat, Y_mat)
  n <- nrow(combined)
  if (n > 2000) {
    idx <- sample(n, 2000)
    combined <- combined[idx, , drop = FALSE]
  }
  d <- as.vector(dist(combined))
  sigma <- median(d)
  # Avoid zero bandwidth
  if (sigma == 0) sigma <- 1
  return(sigma)
}


# Exact MMD^2 computation
# Computes the unbiased MMD^2 estimator using full kernel matrix.
# @param X_mat numeric matrix (n x d)
# @param Y_mat numeric matrix (m x d)
# @param sigma bandwidth
# @param kernel character, "gaussian" or "rational_quadratic"
# @return scalar MMD^2
# @keywords internal
.mmd_exact <- function(X_mat, Y_mat, sigma, kernel) {
  n <- nrow(X_mat)
  m <- nrow(Y_mat)

  # Compute squared distance matrices
  # XX
  X_sq <- rowSums(X_mat^2)
  D_XX_sq <- outer(X_sq, X_sq, "+") - 2 * X_mat %*% t(X_mat)
  D_XX_sq[D_XX_sq < 0] <- 0

  # YY
  Y_sq <- rowSums(Y_mat^2)
  D_YY_sq <- outer(Y_sq, Y_sq, "+") - 2 * Y_mat %*% t(Y_mat)
  D_YY_sq[D_YY_sq < 0] <- 0

  # XY
  D_XY_sq <- outer(X_sq, Y_sq, "+") - 2 * X_mat %*% t(Y_mat)
  D_XY_sq[D_XY_sq < 0] <- 0

  gamma <- 1 / (2 * sigma^2)

  if (kernel == "gaussian") {
    K_XX <- exp(-gamma * D_XX_sq)
    K_YY <- exp(-gamma * D_YY_sq)
    K_XY <- exp(-gamma * D_XY_sq)
  } else {
    # rational quadratic: k(x,y) = (1 + ||x-y||^2 / (2*sigma^2))^{-1}
    K_XX <- 1 / (1 + gamma * D_XX_sq)
    K_YY <- 1 / (1 + gamma * D_YY_sq)
    K_XY <- 1 / (1 + gamma * D_XY_sq)
  }

  # Unbiased estimator: exclude diagonal for within-sample terms
  diag(K_XX) <- 0
  diag(K_YY) <- 0

  term_xx <- sum(K_XX) / (n * (n - 1))
  term_yy <- sum(K_YY) / (m * (m - 1))
  term_xy <- sum(K_XY) / (n * m)

  mmd2 <- term_xx + term_yy - 2 * term_xy
  return(mmd2)
}


# RFF-based MMD^2 approximation
# Approximates MMD^2 using Random Fourier Features (Rahimi & Recht, 2007).
# @param X_mat numeric matrix (n x d)
# @param Y_mat numeric matrix (m x d)
# @param sigma bandwidth
# @param n_features integer, number of random features
# @return scalar MMD^2 (approximate)
# @keywords internal
.mmd_rff <- function(X_mat, Y_mat, sigma, n_features) {
  d <- ncol(X_mat)
  # Random frequencies: omega ~ N(0, 1/sigma^2 * I_d)
  Omega <- matrix(rnorm(d * n_features, mean = 0, sd = 1 / sigma),
                  nrow = d, ncol = n_features)
  # Random offsets: b ~ Uniform(0, 2*pi)
  b <- runif(n_features, min = 0, max = 2 * pi)

  # Feature maps: Z = sqrt(2/D) * cos(X %*% Omega + b)
  Z_X <- sqrt(2 / n_features) * cos(X_mat %*% Omega + matrix(b, nrow = nrow(X_mat), ncol = n_features, byrow = TRUE))
  Z_Y <- sqrt(2 / n_features) * cos(Y_mat %*% Omega + matrix(b, nrow = nrow(Y_mat), ncol = n_features, byrow = TRUE))

  # MMD^2 ~ ||mean(Z_X) - mean(Z_Y)||^2
  diff_means <- colMeans(Z_X) - colMeans(Z_Y)
  mmd2 <- sum(diff_means^2)
  return(mmd2)
}


# ---- S3 methods --------------------------------------------------------------

#' Print method for mmd objects
#'
#' @param x an object of class \code{"mmd"}
#' @param ... additional arguments (ignored)
#' @export
print.mmd <- function(x, ...) {
  cat("Maximum Mean Discrepancy (MMD)\n")
  cat("==============================\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original (X):", x$n_X, "\n")
  cat("  Synthetic (Y):", x$n_Y, "\n")
  cat("  Variables:", x$n_vars, "\n")
  cat("  Standardized:", x$standardized, "\n\n")

  cat("Settings:\n")
  cat("  Kernel:", x$kernel, "\n")
  cat("  Method:", x$method, "\n")
  if (x$method == "rff") {
    cat("  Features:", x$n_features, "\n")
  }
  cat("  Bandwidth (sigma):", sprintf("%.4f", x$sigma), "\n\n")

  cat("Results:\n")
  cat("  MMD^2:    ", sprintf("%.6f", x$mmd2), "\n")
  cat("  Utility:  ", sprintf("%.4f", x$utility_score), "(exp(-MMD^2/sigma^2), higher=better)\n")

  if (!is.null(x$perm_pvalue)) {
    cat("  Perm. p:  ", sprintf("%.4f", x$perm_pvalue),
        sprintf("(%d permutations)", x$n_perm), "\n")
  }
  cat("\n")

  cat("Interpretation:\n")
  if (x$utility_score > 0.95) {
    cat("  EXCELLENT: Distributions are very similar.\n")
  } else if (x$utility_score > 0.80) {
    cat("  GOOD: Distributions are reasonably similar.\n")
  } else if (x$utility_score > 0.50) {
    cat("  MODERATE: Some differences in distribution.\n")
  } else {
    cat("  POOR: Significant distributional differences.\n")
  }

  invisible(x)
}


#' Summary method for mmd objects
#'
#' @param object an object of class \code{"mmd"}
#' @param ... additional arguments (ignored)
#' @export
summary.mmd <- function(object, ...) {
  summ <- list(
    mmd2 = object$mmd2,
    sigma = object$sigma,
    kernel = object$kernel,
    method = object$method,
    n_features = object$n_features,
    n_X = object$n_X,
    n_Y = object$n_Y,
    n_vars = object$n_vars,
    vars = object$vars,
    standardized = object$standardized,
    utility_score = object$utility_score,
    perm_pvalue = object$perm_pvalue,
    perm_null = object$perm_null,
    n_perm = object$n_perm
  )

  class(summ) <- "summary.mmd"
  return(summ)
}


#' Print method for summary.mmd objects
#'
#' @param x an object of class \code{"summary.mmd"}
#' @param ... additional arguments (ignored)
#' @export
print.summary.mmd <- function(x, ...) {
  cat("Summary: Maximum Mean Discrepancy (MMD)\n")
  cat("========================================\n\n")

  cat("Variables (", x$n_vars, "):", paste(x$vars, collapse = ", "), "\n\n")

  cat("Kernel:", x$kernel, " | Method:", x$method)
  if (!is.null(x$n_features)) cat(" | Features:", x$n_features)
  cat("\n")
  cat("Bandwidth (sigma):", sprintf("%.4f", x$sigma), "\n\n")

  cat("MMD^2:         ", sprintf("%.6f", x$mmd2), "\n")
  cat("Utility score: ", sprintf("%.4f", x$utility_score), "\n")

  if (!is.null(x$perm_pvalue)) {
    cat("\nPermutation Test:\n")
    cat("  Permutations:", x$n_perm, "\n")
    cat("  p-value:     ", sprintf("%.4f", x$perm_pvalue), "\n")
    cat("  Null mean:   ", sprintf("%.6f", mean(x$perm_null)), "\n")
    cat("  Null sd:     ", sprintf("%.6f", sd(x$perm_null)), "\n")
  }

  cat("\nSample Sizes: X =", x$n_X, ", Y =", x$n_Y, "\n")
  cat("Standardized:", x$standardized, "\n")

  invisible(x)
}


#' Plot method for mmd objects
#'
#' @param x an object of class \code{"mmd"}
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = permutation null distribution histogram
#'   (requires \code{n_perm} to have been set)
#' @importFrom graphics hist abline legend
#' @export
plot.mmd <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 1)
  show[which] <- TRUE

  if (show[1]) {
    if (is.null(x$perm_null)) {
      message("No permutation results. Run mmd() with n_perm > 0 to enable this plot.")
      return(invisible(x))
    }

    hist(x$perm_null,
         main = "MMD Permutation Test",
         xlab = expression(MMD^2),
         col = "grey80",
         border = "white",
         freq = FALSE,
         ...)
    abline(v = x$mmd2, col = "red", lty = 2, lwd = 2)
    legend("topright",
           legend = c(
             sprintf("Observed MMD^2 = %.4f", x$mmd2),
             sprintf("p-value = %.4f", x$perm_pvalue)
           ),
           col = c("red", NA),
           lty = c(2, NA),
           lwd = c(2, NA),
           bty = "n")
  }

  invisible(x)
}

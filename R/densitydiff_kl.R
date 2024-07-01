#' Kullback-Leibler divergence
#'
#' Kullback-Leibler divergence between two variables or two matrices/data.frames.
#'
#' @details
#' First the probability distributions are estimated before the Kullback-Leiber
#' divergence is calculated.
#'
#' The Kullback-Leibler divergence is defined as:
#' \deqn{KLD(X,Y) = \sum{X \times \log\left(\frac{X}{Y}\right)}}
#'
#' If the data are multivariate, an Kullback-Leibler inspired divergence
#' is calculated through
#'
#' \deqn{KLD(X,Y) = \sum{X \times \log\left(\frac{X}{Y}\right)}}
#' \deqn{KLD(X,Y) = \sum_{i,j}{X_{i,j} ( \log\left(\frac{X_{i,j}{Y_{i,j}} +
#' log(Y) - log(X) \right))}}
#' with \deqn{X = \sum{X}_{i,j}} and \deqn{Y = \sum{Y}_{i,j}}
#'
#' @param Either a numeric vector or a matrix or data frame with numeric entries.
#' @param Either a numeric vector or a matrix or data frame with numeric entries.
#' @author Matthias Templ
#' @return The Kullback-Leibler divergence of \code{X} and \code{Y}.
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#'
#' densitydiff_kl_num(x, y)
#'
#' X <- MASS::mvrnorm(100, mu = c(0,0), Sigma = diag(2))
#' Y <- MASS::mvrnorm(100, mu = c(0,0), Sigma = diag(2))
#' densitydiff_kl_num(X, Y)
#'
#' X <- MASS::mvrnorm(100, mu = c(0,0,0), Sigma = diag(3))
#' Y <- MASS::mvrnorm(100, mu = c(0,0,0), Sigma = diag(3))
#' densitydiff_kl_num(X, Y)
#'
densitydiff_kl_num <- function(X, Y, stepsize = 1000) {
  # Check if X and Y are either numeric vectors, matrices or data frames with numeric values
  is_valid_input <- function(data) {
    if (is.numeric(data) && is.vector(data)) {
      return(TRUE)
    } else if (is.matrix(data) && is.numeric(data)) {
      return(TRUE)
    } else if (is.data.frame(data) && all(sapply(data, is.numeric))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }

  # Check if X and Y have the same structure
  has_same_structure <- function(X, Y) {
    if (is.vector(X) && is.vector(Y)) {
      return(length(X) == length(Y))
    } else if (is.matrix(X) && is.matrix(Y)) {
      return(all(dim(X) == dim(Y)))
    } else if (is.data.frame(X) && is.data.frame(Y)) {
      return(all(dim(X) == dim(Y)) && all(names(X) == names(Y)))
    } else {
      return(FALSE)
    }
  }

  # Apply the checks

  if (!is_valid_input(X)) {
    stop("X should be either a numeric vector, matrix or a data frame with only numeric values.")
  }

  if (!is_valid_input(Y)) {
    stop("Y should be either a numeric vector, matrix or a data frame with only numeric values.")
  }

  if (!has_same_structure(X, Y)) {
    stop("X and Y should have the same structure
         (dimensions and column names).")
  }

  # univariate
  if(is.numeric(X) & is.vector((X))){
    # Kernel Density Estimation
    kde_x <- density(x)
    kde_y <- density(y)

    # Define a sequence of points
    points <- seq(min(kde_x$x, kde_y$x),
                  max(kde_x$x, kde_y$x),
                  length.out = stepsize)

    # Interpolate densities at the sequence of points
    density_X <- approx(kde_x$x, kde_x$y, xout = points)$y
    density_Y <- approx(kde_y$x, kde_y$y, xout = points)$y

    return(sum(density_X * log(density_X / density_Y), na.rm = TRUE))
  }

  # bivariate
  if(ncol(X) == 0){
    dX <- kde2d(x = X[, 1], y = X[, 2])
    dY <- kde2d(x = Y[, 1], y = Y[, 2])
    return(dX$z * log(dX$z / dY$z))
  }

  # tri-variate
  if(ncol(X) == 3){
    dX <- misc3d::kde3d(x = X[, 1], y = X[, 2], z = X[, 3])
    dY <- misc3d::kde3d(x = Y[, 1], y = Y[, 2], z = Y[, 3])
    return(sum(dX$d * log(dX$d / dY$d)))
  }
  # >3 variate
  if(ncol(X) > 3){
   stop("not implemented for more than three-dimensional data")
  }
}

#' Computes the Jensen-Shannon divergence between two probabiliy distributions. TBD
#'
#' @param P A probility distribution (vector summing to one).
#' @param Q A probility distribution (vector summing to one).
#' @return The JSD of \code{P} and \code{Q}.
#'
#' @examples
#' P = prop.table(sample(1:10, 20, replace = TRUE))
#' Q = prop.table(sample(5:15, 20, replace = TRUE))
#'
#' JSD(P,Q)
#' JSD(Q,P)
JSD = function(P, Q) {
  M = (P + Q)/2
  jsd = 0.5 * KLD(P, M) + 0.5 * KLD(Q, M)
  return(jsd)
}

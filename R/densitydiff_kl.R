#' Kullback-Leibler divergence
#'
#' Kullback-Leibler divergence between two variables or two matrices/data.frames.
#'
#' @details First the probability distributions are estimated before the
#' Kullback-Leiber divergence is calculated.
#'
#' The Kullback-Leibler divergence is defined as:
#' \deqn{KLD(X,Y) = \sum{X \times \log\left(\frac{X}{Y}\right)}}
#'
#' If the data are multivariate, a Kullback-Leibler inspired divergence
#' is calculated through
#'
#' \deqn{KLD(X,Y) = \sum{X \times \log\left(\frac{X}{Y}\right)}}
#' \deqn{KLD(X,Y) = \sum_{i,j}{X_{i,j} \left( \log\left(\frac{X_{i,j}}{Y_{i,j}}\right) +
#' \log(Y) - \log(X) \right)}}
#' with \deqn{X = \sum{X_{i,j}}} and \deqn{Y = \sum{Y_{i,j}}}
#'
#' @param X a numeric vector or a matrix or data frame with numeric entries.
#' @param Y a numeric vector or a matrix or data frame with numeric entries.
#' @param stepsize number of interval points where the density is evaluated.
#' @export
#' @importFrom MASS kde2d
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
  # TBD
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
    kde_x <- density(X)
    kde_y <- density(Y)

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
  if(ncol(X) == 2){
    dX <- MASS::kde2d(x = X[, 1], y = X[, 2])
    dY <- MASS::kde2d(x = Y[, 1], y = Y[, 2])
    return(dX$z * log(dX$z / dY$z))
  }

  # tri-variate
  if(ncol(X) == 3){
    if (!requireNamespace("misc3d", quietly = TRUE)) {
      stop("Package 'misc3d' is required for 3D kernel density estimation. Please install it.")
    }
    dX <- misc3d::kde3d(x = X[, 1], y = X[, 2], z = X[, 3])
    dY <- misc3d::kde3d(x = Y[, 1], y = Y[, 2], z = Y[, 3])
    return(sum(dX$d * log(dX$d / dY$d)))
  }
  # >3 variate
  if(ncol(X) > 3){
   stop("not implemented for more than three-dimensional data")
  }
}

#' evaluation statistics
#'
#' Several kinds of evaluation statistics (MAPE, MAE, MSE, RMSE, AIT)
#'
#' @param x numeric vector
#' @param y numeric vector
#' @author Matthias Templ
#' @rdname evaluation_stats
#' @name evaluation_stats
#' @export
#' @aliases mape mae mse rmse ait
#' @return The MAPE, MAE, MSE, RMSE, and (normalized) AIT (Aitchison distance)
#' between two vectors.
#' @examples
#' x <- rnorm(10)
#' y <- rnorm(10)
#' mape(x, y)
#' mae(x, y)
#' mse(x, y)
#' rmse(x, y)
#' ait(x, y)
NULL

#' @rdname evaluation_stats
#' @export
mape <- function(x, y) {
  return(mean(abs((x - y) / x), na.rm = TRUE) * 100)
}

#' @rdname evaluation_stats
#' @export
mae <- function(x, y) {
  return(mean(abs(x - y), na.rm = TRUE))
}

#' @rdname evaluation_stats
#' @export
mse <- function(x, y) {
  return(mean((x - y)^2, na.rm = TRUE))
}

#' @rdname evaluation_stats
#' @export
rmse <- function(x, y) {
  return(sqrt(mse(x, y)))
}

#' @rdname evaluation_stats
#' @export
ait <- function(x, y) {
  n <- length(c(x))
  x <- c(x)
  y <- c(y)
  index <- x != 0 & y != 0
  x <- x[index]
  y <- y[index]
  return(robCompositions::aDist(c(x), c(y)) * (1 / (2 * n)) )
}

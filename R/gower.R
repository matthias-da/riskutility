#' Gower distance between two data frames
#'
#' Average (per observation) Gower distance
#'
#' @param X data frame
#' @param Y data frame
#' @param ... additional arguments passed to methods
#' @return An object of class "gower" containing:
#' \itemize{
#'   \item gower_distance: the average Gower distance
#'   \item n: number of observations
#'   \item n_vars: number of variables
#' }
#'
#' @references
#' Gower, J.C. (1971). A General Coefficient of Similarity and Some of Its
#' Properties. \emph{Biometrics}, 27(4), 857--871.
#'
#' @seealso \code{\link{dcr}}, \code{\link{ims}}
#'
#' @author Matthias Templ. Based on the gowerD function from Alexander Kowarik
#' in the VIM package.
#' @family utility
#' @export
#' @importFrom VIM gowerD
#' @examples
#' # Simple example with mixed data types
#' X <- data.frame(
#'   age = c(25, 30, 35, 40),
#'   income = c(30000, 45000, 50000, 60000),
#'   gender = factor(c("M", "F", "M", "F"))
#' )
#' Y <- data.frame(
#'   age = c(26, 31, 34, 42),
#'   income = c(32000, 44000, 52000, 58000),
#'   gender = factor(c("M", "F", "M", "F"))
#' )
#' gower(X, Y)
#'
gower <- function(X, ...) {
  UseMethod("gower")
}

#' @rdname gower
#' @export
gower.synth_pair <- function(X, ...) {
  gower.default(X = X$original, Y = X$synthetic, ...)
}

#' @rdname gower
#' @export
gower.default <- function(X, Y, ...){
  gd <- VIM::gowerD(X, Y)
  avg_distance <- sum(abs(gd)) / nrow(X)

  result <- list(
    gower_distance = avg_distance,
    n = nrow(X),
    n_vars = ncol(X)
  )
  class(result) <- "gower"
  return(result)
}

#' Print method for gower objects
#'
#' @param x an object of class "gower"
#' @param ... additional arguments (ignored)
#' @export
print.gower <- function(x, ...) {
  cat("Average Gower Distance:", round(x$gower_distance, 6), "\n")
  cat("  Observations:", x$n, "| Variables:", x$n_vars, "\n")
  invisible(x)
}

#' Summary method for gower objects
#'
#' @param object an object of class "gower"
#' @param ... additional arguments (ignored)
#' @return An object of class "summary.gower"
#' @export
summary.gower <- function(object, ...) {
  summ <- list(
    gower_distance = object$gower_distance,
    n = object$n,
    n_vars = object$n_vars
  )
  class(summ) <- "summary.gower"
  summ
}

#' Print method for summary.gower objects
#'
#' @param x an object of class "summary.gower"
#' @param ... additional arguments (ignored)
#' @export
print.summary.gower <- function(x, ...) {
  cat("Summary: Average Gower Distance\n")
  cat("================================\n")
  cat("Distance:", round(x$gower_distance, 6), "\n")
  cat("Observations:", x$n, "\n")
  cat("Variables:", x$n_vars, "\n")
  invisible(x)
}

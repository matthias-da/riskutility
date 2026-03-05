#' Gower distance between two data frames
#'
#' Average (per observation) Gower distance
#'
#' @param X data frame
#' @param Y data frame
#' @param ... additional arguments passed to methods
#' @return A single numeric value: the average Gower distance between the
#'   rows of \code{X} and the corresponding rows of \code{Y}, computed as
#'   the sum of absolute pairwise Gower distances divided by the number of
#'   observations in \code{X}.
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
  return(sum(abs(gd)) / nrow(X))
}

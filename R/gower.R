#' Gower distance between two data frames
#'
#' Average (per observation) Gower distance
#'
#' @param X data frame
#' @param Y data frame
#' @author Matthias Templ. Based on the gowerD function from Alexander Kowarik
#' in the VIM package.
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
gower <- function(X, Y){
  gd <- VIM::gowerD(X, Y)
  return(sum(abs(gd)) / nrow(X))
}

#' Difference of densities with categorical information
#'
#' Computes the compositional Aitchison distance between densities of two categorical variables
#' @param x a factor variable
#' @param y a factor variable
#' @importFrom robCompositions aDist
#' @author Matthias Templ
#' @return The Aitchison distance between two densities
#' @examples
#' # example code
#'
#' x <- sample(1:5, 100, replace = TRUE)
#' y <- sample(1:5, 100, replace = TRUE)

densitydiff_1d_cat <- function(x, y){
  # estimate densities
  tabx <- prop.table(table(x))
  taby <- prop.table(table(y))
  tabx <- table(x)
  taby <- table(y)
  aDist(as.numeric(tabx), as.numeric(taby))
}

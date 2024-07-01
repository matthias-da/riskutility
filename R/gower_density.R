#' Gower density difference
#'
#' Computes the Gower distance based on densities
#' @param X data frame
#' @param Y data frame
#' @author Matthias Templ
gower_density <- function(X, Y){
  ## wie bei kategorischen?
  # ...
  gd <- VIM::gowerD(X, Y)
  return(sum(abs(gd)) / nrow(X))
}

#' Gower distance between two data frames
#'
#' @param X data frame
#' @param Y data frame
#' @author Matthias Templ
#' @examples
#' # example code
#' data(eusilc13puf, package = "simPop")
#' X <- Y <- eusilc
#' sdc <- createSdcObj(X,
#'                     keyVars = c("db040", "hsize", "pb220a",
#'                                 "rb090", "pl030", "age"),
#'                     numVars = c("eqIncome"),
#'                     pramVars = "db040",
#'                     weightVar = "rb050", hhId = "db030")
#' sdc <- pram(sdc)
#' sdc <- microaggregation(sdc, strata="db040")
#' Y <- extractManipData(sdc)
#' gower(X[, c("eqIncome","db040")], Y[, c("eqIncome","db040")])
#'
gower <- function(X, Y){
  gd <- VIM::gowerD(X, Y)
  return(sum(abs(gd)) / nrow(X))
}

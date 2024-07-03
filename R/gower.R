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
#' @importFrom sdcMicro createSdcObj pram microaggregation extractManipData
#' @examples
#' # example code
#' data(eusilc13puf, package = "simPop")
#' X <- Y <- eusilc13puf
#' sdc <- sdcMicro::createSdcObj(X,
#'                     keyVars = c("db040", "hsize", "pb220a",
#'                                 "rb090", "pl031", "age"),
#'                     numVars = c("pgrossIncome"),
#'                     pramVars = "db040",
#'                     weightVar = "rb050", hhId = "db030")
#' sdc <- sdcMicro::pram(sdc)
#' sdc <- sdcMicro::microaggregation(sdc, strata="db040")
#' Y <- sdcMicro::extractManipData(sdc)
#' gower(X[, c("pgrossIncome","db040")], Y[, c("pgrossIncome","db040")])
#'
gower <- function(X, Y){
  gd <- VIM::gowerD(X, Y)
  return(sum(abs(gd)) / nrow(X))
}

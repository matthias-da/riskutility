#' Gower distance between two data frames
#'
#' Average (per observation) Gower distance between an original and a
#' synthetic/anonymized data frame.  Returns an S3 object of class
#' \code{"gower"} with print, summary, and plot methods.
#'
#' @param X data frame (or a \code{synth_pair} object)
#' @param Y data frame
#' @param ... additional arguments passed to methods
#' @return An object of class \code{"gower"} containing:
#' \itemize{
#'   \item \code{gower_distance}: the average Gower distance (lower = more similar)
#'   \item \code{utility_score}: \code{1 - gower_distance}, in \eqn{[0, 1]}
#'     (higher = better utility)
#'   \item \code{n_records}: number of observations
#'   \item \code{n_variables}: number of variables
#'   \item \code{n}: alias for \code{n_records} (for backward compatibility)
#'   \item \code{n_vars}: alias for \code{n_variables} (for backward compatibility)
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
#' result <- gower(X, Y)
#' result
#' summary(result)
#' \donttest{
#' plot(result)
#' }
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

  nr <- nrow(X)
  nv <- ncol(X)

  result <- list(
    gower_distance = avg_distance,
    utility_score  = 1 - avg_distance,
    n_records      = nr,
    n_variables    = nv,
    n              = nr,
    n_vars         = nv
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
  cat("  Utility score:", round(x$utility_score, 6),
      "(1 - distance, higher = better)\n")
  cat("  Observations:", x$n_records, "| Variables:", x$n_variables, "\n")
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
    utility_score  = object$utility_score,
    n_records      = object$n_records,
    n_variables    = object$n_variables,
    n              = object$n_records,
    n_vars         = object$n_variables
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
  cat("Distance:      ", round(x$gower_distance, 6), "\n")
  cat("Utility score: ", round(x$utility_score, 6),
      "(1 - distance, higher = better)\n")
  cat("Observations:  ", x$n_records, "\n")
  cat("Variables:     ", x$n_variables, "\n")
  invisible(x)
}

#' Plot method for gower objects
#'
#' Visualizes the average Gower distance as a simple bar plot with a reference
#' scale from 0 (identical) to 1 (maximally different).
#'
#' @param x an object of class \code{"gower"}
#' @param y not used
#' @param ... additional arguments passed to \code{\link[graphics]{barplot}}
#' @param which integer, which plot: 1 = bar chart of distance and utility score
#' @importFrom graphics barplot abline axis mtext
#' @export
plot.gower <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 1)
  show[which] <- TRUE

  if (show[1]) {
    vals <- c(Distance = x$gower_distance, Utility = x$utility_score)
    cols <- c("coral", "steelblue")

    bp <- barplot(vals,
                  col = cols,
                  ylim = c(0, 1),
                  ylab = "Score",
                  main = "Average Gower Distance",
                  las = 1,
                  ...)

    # reference lines
    abline(h = 0.5, lty = 2, col = "grey50")

    # annotate bar values
    text(bp, vals + 0.04, labels = sprintf("%.4f", vals), cex = 0.9)

    mtext(paste0("n = ", x$n_records, " obs, ", x$n_variables, " vars"),
          side = 1, line = 2.5, cex = 0.8)
  }

  invisible(x)
}

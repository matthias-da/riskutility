#' Gower distance between two data frames
#'
#' Mean per-observation Gower distance between an original data frame and a
#' \emph{row-corresponding} anonymized version (e.g. a perturbative method such
#' as noise addition, microaggregation, or rank swapping, where record \eqn{i}
#' of the anonymized data corresponds to record \eqn{i} of the original).
#' Returns an S3 object of class \code{"gower"} with print, summary, and plot
#' methods.
#'
#' @details The two data frames must have the same number of rows; rows are
#' compared pairwise, so the reported \code{gower_distance} is
#' \eqn{\frac{1}{n}\sum_{i=1}^{n} d_G(X_i, Y_i)}, the mean Gower distance between
#' corresponding records, in \eqn{[0, 1]}. Identical data therefore yields
#' distance \code{0} and \code{utility_score} \code{1}. This measure assumes
#' record correspondence (perturbative anonymization); for fully synthetic data
#' without such correspondence use a distributional measure such as
#' \code{\link{energy_distance}}, \code{\link{mmd}}, or \code{compare_wasserstein}.
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
  if (nrow(X) != nrow(Y)) {
    stop("gower() compares row-corresponding data frames; X and Y must have ",
         "the same number of rows. For synthetic data without record ",
         "correspondence, use energy_distance(), mmd(), or compare_wasserstein().")
  }
  # Pairwise (diagonal) Gower distance: record i of X vs record i of Y.
  # gowerD() returns the full n x n cross-distance matrix; the diagonal holds
  # the corresponding-record distances, whose mean is the per-observation loss.
  gd <- VIM::gowerD(X, Y)
  avg_distance <- mean(diag(as.matrix(gd)))

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
#' @return The input object, invisibly.
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
#' @return The input object, invisibly.
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
#' @return No return value; called for the side effect of producing a plot.
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

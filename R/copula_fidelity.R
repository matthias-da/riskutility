#' Copula Fidelity for Dependence Structure Comparison
#'
#' Compares the dependence structure between original and synthetic data using
#' empirical copulas. For each pair of numeric variables, the bivariate empirical
#' copula is estimated via rank-transformation and compared using the
#' Cramer-von Mises (CvM) statistic on a grid. The mean pairwise CvM distance
#' provides an overall measure of how well the synthetic data preserves
#' the joint dependence structure.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the synthetic/anonymized dataset.
#' @param vars Character vector of numeric variable names to compare.
#'   If NULL (default), all common numeric variables are used. Categorical
#'   variables are skipped with a message.
#' @param n_grid Integer, grid resolution for CvM evaluation. Default 50.
#'   Higher values give more accurate CvM estimates but increase computation.
#' @param na.rm Logical, whether to remove rows with NA values. Default TRUE.
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class \code{"copula_fidelity"} containing:
#' \itemize{
#'   \item \code{mean_cvm}: mean of all pairwise CvM distances
#'   \item \code{utility_score}: transformed score in \eqn{[0,1]} (higher = better),
#'     computed as \code{1 / (1 + mean_cvm * 100)}
#'   \item \code{pairwise}: data.frame with columns \code{var1}, \code{var2},
#'     \code{cvm_distance} for each pair of variables
#'   \item \code{n_vars}: number of numeric variables used
#'   \item \code{vars}: names of numeric variables used
#'   \item \code{n_X}: number of rows in X (after NA removal)
#'   \item \code{n_Y}: number of rows in Y (after NA removal)
#'   \item \code{n_grid}: grid resolution used
#' }
#'
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item Each numeric variable is rank-transformed to \eqn{[0, 1]}:
#'     \code{rank(x) / (length(x) + 1)}.
#'   \item For each pair of variables \eqn{(i, j)}, the bivariate empirical CDF
#'     is computed for both original and synthetic data:
#'     \deqn{F(u_1, u_2) = \frac{1}{n} \sum_{k=1}^{n} 1(U_{k,i} \le u_1, U_{k,j} \le u_2)}
#'   \item The Cramer-von Mises statistic measures the integrated squared
#'     difference between the two empirical copulas on a regular grid:
#'     \deqn{CvM_{ij} = \text{mean}\left( (F_{\text{orig}}(u_1, u_2) - F_{\text{syn}}(u_1, u_2))^2 \right)}
#'   \item The overall fidelity is summarized as the mean of all pairwise CvM
#'     distances. A utility score is computed as \code{1 / (1 + mean_cvm * 100)}.
#' }
#'
#' This measure is invariant to monotone transformations of individual variables
#' (since only ranks matter) and focuses purely on the dependence structure.
#' It complements marginal distribution comparisons (e.g., Wasserstein, Hellinger)
#' by assessing whether the joint structure is preserved.
#'
#' @seealso \code{\link{compare_correlation_matrices}} for linear dependence,
#'   \code{\link{energy_distance}} for multivariate distribution comparison,
#'   \code{\link{mmd}} for kernel-based comparison
#'
#' @references
#' Genest, C. and Remillard, B. (2008). Validity of the parametric bootstrap
#' for goodness-of-fit testing in semiparametric models. Annales de l'Institut
#' Henri Poincare, Probabilites et Statistiques, 44(6), 1096-1127.
#'
#' Deheuvels, P. (1979). La fonction de dependance empirique et ses proprietes.
#' Un test non parametrique d'independance. Academie Royale de Belgique.
#' Bulletin de la Classe des Sciences, 5e Serie, 65, 274-292.
#'
#' @author Matthias Templ
#' @family utility
#' @export
#' @importFrom stats complete.cases
#'
#' @examples
#' set.seed(123)
#' # Original data with dependence
#' n <- 300
#' x1 <- rnorm(n)
#' x2 <- 0.8 * x1 + rnorm(n, sd = 0.6)
#' x3 <- rnorm(n)
#' X <- data.frame(x1 = x1, x2 = x2, x3 = x3)
#'
#' # Good synthetic data (preserves dependence)
#' s1 <- rnorm(n)
#' s2 <- 0.8 * s1 + rnorm(n, sd = 0.6)
#' s3 <- rnorm(n)
#' Y_good <- data.frame(x1 = s1, x2 = s2, x3 = s3)
#'
#' # Poor synthetic data (independent variables)
#' Y_poor <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
#'
#' result_good <- copula_fidelity(X, Y_good)
#' print(result_good)
#'
#' result_poor <- copula_fidelity(X, Y_poor)
#' print(result_poor)
#'
#' \donttest{
#' # Heatmap of pairwise copula distances
#' plot(result_poor)
#' }
copula_fidelity <- function(X, ...) {
  UseMethod("copula_fidelity")
}

#' @rdname copula_fidelity
#' @export
copula_fidelity.synth_pair <- function(X, ...) {
  copula_fidelity.default(
    X = X$original,
    Y = X$synthetic,
    vars = X$num_vars,  # Use auto-detected numeric variables
    ...
  )
}

#' @rdname copula_fidelity
#' @export
copula_fidelity.default <- function(X, Y,
                                    vars = NULL,
                                    n_grid = 50L,
                                    na.rm = TRUE,
                                    ...) {

  # Convert to data.frame if needed
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)

  n_grid <- as.integer(n_grid)

  # Auto-detect numeric variables if not specified
  if (is.null(vars)) {
    num_vars_X <- names(X)[sapply(X, is.numeric)]
    num_vars_Y <- names(Y)[sapply(Y, is.numeric)]
    vars <- intersect(num_vars_X, num_vars_Y)

    # Check for skipped categorical variables
    all_common <- intersect(names(X), names(Y))
    cat_vars <- setdiff(all_common, vars)
    if (length(cat_vars) > 0) {
      message("Skipping non-numeric variables: ", paste(cat_vars, collapse = ", "))
    }
  }

  if (length(vars) < 2) {
    stop("At least 2 numeric variables are required for copula fidelity (need pairwise comparisons).")
  }

  # Check variables exist
  missing_X <- setdiff(vars, names(X))
  missing_Y <- setdiff(vars, names(Y))
  if (length(missing_X) > 0) {
    stop(paste("Variables missing in X:", paste(missing_X, collapse = ", ")))
  }
  if (length(missing_Y) > 0) {
    stop(paste("Variables missing in Y:", paste(missing_Y, collapse = ", ")))
  }

  # Subset to selected variables
  X <- X[, vars, drop = FALSE]
  Y <- Y[, vars, drop = FALSE]

  # Handle missing values
  if (na.rm) {
    X <- X[complete.cases(X), , drop = FALSE]
    Y <- Y[complete.cases(Y), , drop = FALSE]
  }

  if (nrow(X) == 0) stop("No complete cases in X after removing NAs.")
  if (nrow(Y) == 0) stop("No complete cases in Y after removing NAs.")

  n_X <- nrow(X)
  n_Y <- nrow(Y)

  # Rank-transform each variable to [0, 1]
  U_X <- as.data.frame(lapply(X, function(col) rank(col) / (length(col) + 1)))
  U_Y <- as.data.frame(lapply(Y, function(col) rank(col) / (length(col) + 1)))

  # Build the evaluation grid
  grid_pts <- seq(0, 1, length.out = n_grid)

  # Compute pairwise CvM distances
  n_vars <- length(vars)
  pairs <- combn(n_vars, 2)
  n_pairs <- ncol(pairs)

  pairwise_results <- data.frame(
    var1 = character(n_pairs),
    var2 = character(n_pairs),
    cvm_distance = numeric(n_pairs)
  )

  for (p in seq_len(n_pairs)) {
    i <- pairs[1, p]
    j <- pairs[2, p]

    cvm_val <- .bivariate_cvm(
      u1_X = U_X[[i]], u2_X = U_X[[j]],
      u1_Y = U_Y[[i]], u2_Y = U_Y[[j]],
      grid = grid_pts
    )

    pairwise_results$var1[p] <- vars[i]
    pairwise_results$var2[p] <- vars[j]
    pairwise_results$cvm_distance[p] <- cvm_val
  }

  # Summary statistics
  mean_cvm <- mean(pairwise_results$cvm_distance)
  utility_score <- 1 / (1 + mean_cvm * 100)

  result <- list(
    mean_cvm = mean_cvm,
    utility_score = utility_score,
    pairwise = pairwise_results,
    n_vars = n_vars,
    vars = vars,
    n_X = n_X,
    n_Y = n_Y,
    n_grid = n_grid
  )

  class(result) <- "copula_fidelity"
  return(result)
}


# ---- Internal helpers --------------------------------------------------------

# Bivariate Cramer-von Mises statistic
#
# Computes the CvM distance between two bivariate empirical copulas
# evaluated on a grid.
#
# @param u1_X rank-transformed variable 1 from X (in [0,1])
# @param u2_X rank-transformed variable 2 from X (in [0,1])
# @param u1_Y rank-transformed variable 1 from Y (in [0,1])
# @param u2_Y rank-transformed variable 2 from Y (in [0,1])
# @param grid numeric vector of grid points in [0,1]
# @return scalar CvM statistic (mean of squared differences)
# @keywords internal
.bivariate_cvm <- function(u1_X, u2_X, u1_Y, u2_Y, grid) {
  n_X <- length(u1_X)
  n_Y <- length(u1_Y)
  n_g <- length(grid)

  # Vectorized computation: for each grid point (g1, g2),

  # compute the empirical CDF as proportion of points <= (g1, g2)
  # We use outer products for efficiency

  # Pre-compute indicator matrices: rows = observations, cols = grid points
  # ind1_X[k, g] = 1 if u1_X[k] <= grid[g]
  ind1_X <- outer(u1_X, grid, "<=")  # n_X x n_g

  ind2_X <- outer(u2_X, grid, "<=")  # n_X x n_g
  ind1_Y <- outer(u1_Y, grid, "<=")  # n_Y x n_g
  ind2_Y <- outer(u2_Y, grid, "<=")  # n_Y x n_g

  # For each grid combination (g1, g2), count joint indicator
  # F_X(g1, g2) = mean(u1_X <= g1 & u2_X <= g2)
  # = (1/n_X) * sum_k ind1_X[k, g1] * ind2_X[k, g2]
  # = (1/n_X) * t(ind1_X) %*% ind2_X  evaluated at (g1, g2)

  F_X <- crossprod(ind1_X, ind2_X) / n_X  # n_g x n_g matrix
  F_Y <- crossprod(ind1_Y, ind2_Y) / n_Y  # n_g x n_g matrix

  # CvM = mean of squared differences
  cvm <- mean((F_X - F_Y)^2)
  return(cvm)
}


# ---- S3 methods --------------------------------------------------------------

#' Print method for copula_fidelity objects
#'
#' @param x an object of class \code{"copula_fidelity"}
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.copula_fidelity <- function(x, ...) {
  cat("Copula Fidelity - Empirical Copula Dependence Comparison\n")
  cat("========================================================\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original (X):", x$n_X, "\n")
  cat("  Synthetic (Y):", x$n_Y, "\n")
  cat("  Variables:", x$n_vars, "\n")
  cat("  Grid resolution:", x$n_grid, "\n\n")

  cat("Results:\n")
  cat("  Mean CvM distance:", sprintf("%.6f", x$mean_cvm), "\n")
  cat("  Utility score:    ", sprintf("%.4f", x$utility_score),
      "(1/(1+CvM*100), higher=better)\n\n")

  cat("Pairwise CvM Distances:\n")
  for (i in seq_len(nrow(x$pairwise))) {
    cat(sprintf("  %-12s vs %-12s: %.6f\n",
                x$pairwise$var1[i], x$pairwise$var2[i],
                x$pairwise$cvm_distance[i]))
  }
  cat("\n")

  cat("Interpretation:\n")
  if (x$utility_score > 0.95) {
    cat("  EXCELLENT: Dependence structure is very well preserved.\n")
  } else if (x$utility_score > 0.80) {
    cat("  GOOD: Dependence structure is reasonably preserved.\n")
  } else if (x$utility_score > 0.50) {
    cat("  MODERATE: Some differences in dependence structure.\n")
  } else {
    cat("  POOR: Significant differences in dependence structure.\n")
  }

  invisible(x)
}


#' Summary method for copula_fidelity objects
#'
#' @param object an object of class \code{"copula_fidelity"}
#' @param ... additional arguments (ignored)
#' @return A list of summary statistics for the corresponding object.
#' @export
summary.copula_fidelity <- function(object, ...) {
  summ <- list(
    mean_cvm = object$mean_cvm,
    utility_score = object$utility_score,
    pairwise = object$pairwise,
    max_cvm = max(object$pairwise$cvm_distance),
    min_cvm = min(object$pairwise$cvm_distance),
    sd_cvm = if (nrow(object$pairwise) > 1) sd(object$pairwise$cvm_distance) else 0,
    n_pairs = nrow(object$pairwise),
    n_vars = object$n_vars,
    vars = object$vars,
    n_X = object$n_X,
    n_Y = object$n_Y,
    n_grid = object$n_grid
  )

  class(summ) <- "summary.copula_fidelity"
  return(summ)
}


#' Print method for summary.copula_fidelity objects
#'
#' @param x an object of class \code{"summary.copula_fidelity"}
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.copula_fidelity <- function(x, ...) {
  cat("Summary: Copula Fidelity\n")
  cat("========================\n\n")

  cat("Variables (", x$n_vars, "):", paste(x$vars, collapse = ", "), "\n")
  cat("Pairs compared:", x$n_pairs, "\n")
  cat("Grid resolution:", x$n_grid, "\n\n")

  cat("CvM Distance Summary:\n")
  cat("  Mean:", sprintf("%.6f", x$mean_cvm), "\n")
  cat("  Min: ", sprintf("%.6f", x$min_cvm), "\n")
  cat("  Max: ", sprintf("%.6f", x$max_cvm), "\n")
  cat("  SD:  ", sprintf("%.6f", x$sd_cvm), "\n\n")

  cat("Utility score:", sprintf("%.4f", x$utility_score), "\n\n")

  cat("Pairwise Details:\n")
  print(x$pairwise, row.names = FALSE)
  cat("\n")

  cat("Sample Sizes: X =", x$n_X, ", Y =", x$n_Y, "\n")

  invisible(x)
}


#' Plot method for copula_fidelity objects
#'
#' @param x an object of class \code{"copula_fidelity"}
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = heatmap of pairwise copula distances
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2
#'   labs theme_minimal theme element_text
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.copula_fidelity <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 1)
  show[which] <- TRUE

  if (show[1]) {
    # Build a symmetric matrix of CvM distances
    vars <- x$vars
    n_v <- length(vars)
    mat <- matrix(0, nrow = n_v, ncol = n_v,
                  dimnames = list(vars, vars))

    for (i in seq_len(nrow(x$pairwise))) {
      v1 <- x$pairwise$var1[i]
      v2 <- x$pairwise$var2[i]
      mat[v1, v2] <- x$pairwise$cvm_distance[i]
      mat[v2, v1] <- x$pairwise$cvm_distance[i]
    }

    # Use reshape2::melt() explicitly to avoid data.table conflict
    df <- reshape2::melt(mat, varnames = c("Var1", "Var2"), value.name = "CvM")

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Var1, y = .data$Var2,
                                           fill = .data$CvM)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.4f", .data$CvM)),
                         size = 3.5) +
      ggplot2::scale_fill_gradient2(
        low = "steelblue", mid = "white", high = "firebrick",
        midpoint = median(x$pairwise$cvm_distance),
        name = "CvM Distance"
      ) +
      ggplot2::labs(
        title = "Pairwise Copula CvM Distances",
        subtitle = sprintf("Mean CvM = %.4f | Utility = %.4f",
                           x$mean_cvm, x$utility_score),
        x = NULL, y = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(face = "bold")
      )

    print(p)
  }

  invisible(x)
}

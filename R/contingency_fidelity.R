#' Contingency Table Fidelity for Categorical Dependence Comparison
#'
#' Compares bivariate contingency tables (joint distributions) for all pairs of
#' categorical variables between original and synthetic data. For each pair, the
#' total variation (TV) distance between the proportion tables is computed. The
#' mean TV distance across all pairs summarizes how well the synthetic data
#' preserves the bivariate categorical dependence structure.
#'
#' Use this when the data contains multiple categorical variables and you want
#' to check whether their bivariate dependence structure is preserved. This
#' complements \code{\link{copula_fidelity}}, which handles numeric pairs.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the synthetic/anonymized dataset.
#' @param vars Character vector of categorical (factor/character) variable names
#'   to compare. If NULL (default), all common categorical variables are used.
#'   Numeric variables are skipped with a message.
#' @param na.rm Logical, whether to remove rows with NA values. Default TRUE.
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class \code{"contingency_fidelity"} containing:
#' \itemize{
#'   \item \code{mean_tv}: mean of all pairwise total variation distances
#'   \item \code{utility_score}: \code{1 - mean_tv}, in \eqn{[0, 1]} (higher = better)
#'   \item \code{pairwise}: data.frame with columns \code{var1}, \code{var2},
#'     \code{tv_distance} for each pair of variables
#'   \item \code{n_vars}: number of categorical variables used
#'   \item \code{vars}: names of categorical variables used
#'   \item \code{n_X}: number of rows in X (after NA removal)
#'   \item \code{n_Y}: number of rows in Y (after NA removal)
#' }
#'
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item Auto-detect categorical variables (factor or character). Numeric
#'     variables are skipped with an informational message.
#'   \item For each pair of categorical variables \eqn{(i, j)}:
#'     \itemize{
#'       \item Factor levels are aligned: the union of levels from both datasets
#'         is used so that cells present in one dataset but not the other receive
#'         a count of zero.
#'       \item The bivariate contingency table (proportions) is computed for both
#'         original and synthetic data.
#'       \item The total variation distance is computed as
#'         \deqn{TV = \frac{1}{2} \sum_{c} |P_{\mathrm{orig}}(c) - P_{\mathrm{syn}}(c)|}
#'         where the sum runs over all cells \eqn{c} in the cross-tabulation.
#'     }
#'   \item The overall fidelity is summarized as \code{mean_tv}, the mean of all
#'     pairwise TV distances. The utility score is \code{1 - mean_tv}.
#' }
#'
#' Total variation distance is bounded in \eqn{[0, 1]}: 0 means identical
#' joint distributions and 1 means completely non-overlapping distributions.
#' Pairs where a variable has only a single level in both datasets are skipped
#' and receive \code{NA}.
#'
#' \strong{Interpretation (heuristic thresholds):}
#' \itemize{
#'   \item utility_score > 0.95: EXCELLENT -- dependence structure very well preserved
#'   \item utility_score > 0.80: GOOD -- reasonably preserved
#'   \item utility_score > 0.50: MODERATE -- some differences
#'   \item utility_score <= 0.50: POOR -- significant differences
#' }
#'
#' Note that this measure only captures bivariate (pairwise) categorical
#' associations. Higher-order interactions (3-way or more) are not assessed.
#'
#' @seealso \code{\link{copula_fidelity}} for numeric dependence comparison,
#'   \code{\link{compare_chisq_gof}} for univariate categorical comparison,
#'   \code{\link{hellinger}} for univariate Hellinger distance,
#'   \code{\link{regression_fidelity}} for analysis-specific fidelity,
#'   \code{\link{subgroup_utility}} for stratified utility assessment
#'
#' @references
#' Snoke, J., Raab, G. M., Nowok, B., Dibben, C., and Slavkovic, A. (2018).
#' General and Specific Utility Measures for Synthetic Data. Journal of the
#' Royal Statistical Society: Series A, 181(3), 663-688.
#'
#' @author Matthias Templ
#' @family utility
#' @export
#' @importFrom stats complete.cases
#'
#' @examples
#' set.seed(42)
#' n <- 500
#' # Original data with dependence between cat variables
#' gender <- sample(c("M", "F"), n, replace = TRUE)
#' region <- sample(c("N", "S", "E", "W"), n, replace = TRUE)
#' edu <- ifelse(gender == "M",
#'               sample(c("low", "mid", "high"), n, replace = TRUE,
#'                      prob = c(0.2, 0.5, 0.3)),
#'               sample(c("low", "mid", "high"), n, replace = TRUE,
#'                      prob = c(0.4, 0.4, 0.2)))
#' X <- data.frame(gender = gender, region = region, education = edu,
#'                  stringsAsFactors = TRUE)
#'
#' # Good synthetic data (preserves dependence)
#' gender2 <- sample(c("M", "F"), n, replace = TRUE)
#' edu2 <- ifelse(gender2 == "M",
#'                sample(c("low", "mid", "high"), n, replace = TRUE,
#'                       prob = c(0.2, 0.5, 0.3)),
#'                sample(c("low", "mid", "high"), n, replace = TRUE,
#'                       prob = c(0.4, 0.4, 0.2)))
#' Y_good <- data.frame(gender = gender2, region = sample(region),
#'                       education = edu2, stringsAsFactors = TRUE)
#'
#' result <- contingency_fidelity(X, Y_good)
#' print(result)
#' summary(result)
#'
#' # Using synth_pair
#' pair <- synth_pair(X, Y_good)
#' result2 <- contingency_fidelity(pair)
#'
#' # Selecting specific variables
#' result3 <- contingency_fidelity(X, Y_good, vars = c("gender", "education"))
#'
#' \donttest{
#' plot(result)
#' }
contingency_fidelity <- function(X, ...) {
  UseMethod("contingency_fidelity")
}

#' @rdname contingency_fidelity
#' @export
contingency_fidelity.synth_pair <- function(X, ...) {
  contingency_fidelity.default(
    X = X$original,
    Y = X$synthetic,
    vars = X$cat_vars,  # Use auto-detected categorical variables
    ...
  )
}

#' @rdname contingency_fidelity
#' @export
contingency_fidelity.default <- function(X, Y,
                                         vars = NULL,
                                         na.rm = TRUE,
                                         ...) {

  # Convert to data.frame if needed
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)

  # Auto-detect categorical variables if not specified
  if (is.null(vars)) {
    is_cat_X <- sapply(X, function(col) is.factor(col) || is.character(col))
    is_cat_Y <- sapply(Y, function(col) is.factor(col) || is.character(col))
    cat_vars_X <- names(X)[is_cat_X]
    cat_vars_Y <- names(Y)[is_cat_Y]
    vars <- intersect(cat_vars_X, cat_vars_Y)

    # Check for skipped numeric variables
    all_common <- intersect(names(X), names(Y))
    num_vars <- setdiff(all_common, vars)
    if (length(num_vars) > 0) {
      message("Skipping non-categorical variables: ", paste(num_vars, collapse = ", "))
    }
  }

  if (length(vars) < 2) {
    stop("At least 2 categorical variables are required for contingency fidelity (need pairwise comparisons).")
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

  # Ensure all columns are factors with aligned levels
  for (v in vars) {
    lvls <- sort(union(
      if (is.factor(X[[v]])) levels(X[[v]]) else unique(as.character(X[[v]])),
      if (is.factor(Y[[v]])) levels(Y[[v]]) else unique(as.character(Y[[v]]))
    ))
    X[[v]] <- factor(X[[v]], levels = lvls)
    Y[[v]] <- factor(Y[[v]], levels = lvls)
  }

  # Compute pairwise TV distances
  n_vars <- length(vars)
  pairs <- combn(n_vars, 2)
  n_pairs <- ncol(pairs)

  pairwise_results <- data.frame(
    var1 = character(n_pairs),
    var2 = character(n_pairs),
    tv_distance = numeric(n_pairs)
  )

  for (p in seq_len(n_pairs)) {
    i <- pairs[1, p]
    j <- pairs[2, p]

    tv_val <- .contingency_tv(X[[i]], X[[j]], Y[[i]], Y[[j]])

    pairwise_results$var1[p] <- vars[i]
    pairwise_results$var2[p] <- vars[j]
    pairwise_results$tv_distance[p] <- tv_val
  }

  # Summary statistics (ignore NA pairs from single-level skips)
  valid_tv <- pairwise_results$tv_distance[!is.na(pairwise_results$tv_distance)]
  mean_tv <- if (length(valid_tv) > 0) mean(valid_tv) else NA_real_
  utility_score <- if (!is.na(mean_tv)) 1 - mean_tv else NA_real_

  result <- list(
    mean_tv = mean_tv,
    utility_score = utility_score,
    pairwise = pairwise_results,
    n_vars = n_vars,
    vars = vars,
    n_X = n_X,
    n_Y = n_Y
  )

  class(result) <- "contingency_fidelity"
  return(result)
}


# ---- Internal helpers --------------------------------------------------------

# Total variation distance between two bivariate contingency tables.
# Computes TV = 0.5 * sum(|P_orig - P_synth|) where P is the proportion
# table (cross-tabulation normalized to sum to 1).
# Returns NA if a variable has < 2 observed levels in both datasets.
.contingency_tv <- function(x1_orig, x2_orig, x1_synth, x2_synth) {
  # Skip pairs where a variable has < 2 levels in both datasets
  if (nlevels(droplevels(x1_orig)) < 2 && nlevels(droplevels(x1_synth)) < 2) return(NA_real_)
  if (nlevels(droplevels(x2_orig)) < 2 && nlevels(droplevels(x2_synth)) < 2) return(NA_real_)

  # Build contingency tables (counts)
  tab_orig <- table(x1_orig, x2_orig)
  tab_synth <- table(x1_synth, x2_synth)

  # Normalize to proportions
  p_orig <- tab_orig / sum(tab_orig)
  p_synth <- tab_synth / sum(tab_synth)

  # Total variation distance
  tv <- 0.5 * sum(abs(p_orig - p_synth))
  return(tv)
}


# ---- S3 methods --------------------------------------------------------------

#' Print method for contingency_fidelity objects
#'
#' @param x an object of class \code{"contingency_fidelity"}
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.contingency_fidelity <- function(x, ...) {
  cat("Contingency Fidelity - Categorical Dependence Comparison\n")
  cat("========================================================\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original (X):", x$n_X, "\n")
  cat("  Synthetic (Y):", x$n_Y, "\n")
  cat("  Categorical variables:", x$n_vars, "\n\n")

  cat("Results:\n")
  cat("  Mean TV distance:", sprintf("%.6f", x$mean_tv), "\n")
  cat("  Utility score:   ", sprintf("%.4f", x$utility_score),
      "(1 - mean_tv, higher=better)\n\n")

  n_valid <- sum(!is.na(x$pairwise$tv_distance))
  n_total <- nrow(x$pairwise)
  cat("Pairwise TV Distances (", n_valid, "/", n_total, " pairs):\n", sep = "")
  for (i in seq_len(nrow(x$pairwise))) {
    tv_str <- if (is.na(x$pairwise$tv_distance[i])) "NA (skipped)" else sprintf("%.6f", x$pairwise$tv_distance[i])
    cat(sprintf("  %-12s vs %-12s: %s\n",
                x$pairwise$var1[i], x$pairwise$var2[i], tv_str))
  }
  cat("\n")

  cat("Interpretation:\n")
  if (is.na(x$utility_score)) {
    cat("  No valid pairs to evaluate.\n")
  } else if (x$utility_score > 0.95) {
    cat("  EXCELLENT: Categorical dependence structure is very well preserved.\n")
  } else if (x$utility_score > 0.80) {
    cat("  GOOD: Categorical dependence structure is reasonably preserved.\n")
  } else if (x$utility_score > 0.50) {
    cat("  MODERATE: Some differences in categorical dependence structure.\n")
  } else {
    cat("  POOR: Significant differences in categorical dependence structure.\n")
  }

  invisible(x)
}


#' Summary method for contingency_fidelity objects
#'
#' @param object an object of class \code{"contingency_fidelity"}
#' @param ... additional arguments (ignored)
#' @return A list of summary statistics for the corresponding object.
#' @export
summary.contingency_fidelity <- function(object, ...) {
  valid_tv <- object$pairwise$tv_distance[!is.na(object$pairwise$tv_distance)]
  summ <- list(
    mean_tv = object$mean_tv,
    utility_score = object$utility_score,
    pairwise = object$pairwise,
    max_tv = if (length(valid_tv) > 0) max(valid_tv) else NA_real_,
    min_tv = if (length(valid_tv) > 0) min(valid_tv) else NA_real_,
    sd_tv = if (length(valid_tv) > 1) sd(valid_tv) else NA_real_,
    n_pairs = nrow(object$pairwise),
    n_valid_pairs = length(valid_tv),
    n_vars = object$n_vars,
    vars = object$vars,
    n_X = object$n_X,
    n_Y = object$n_Y
  )

  class(summ) <- "summary.contingency_fidelity"
  return(summ)
}


#' Print method for summary.contingency_fidelity objects
#'
#' @param x an object of class \code{"summary.contingency_fidelity"}
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.contingency_fidelity <- function(x, ...) {
  cat("Summary: Contingency Fidelity\n")
  cat("==============================\n\n")

  cat("Variables (", x$n_vars, "):", paste(x$vars, collapse = ", "), "\n")
  cat("Pairs compared:", x$n_pairs, "(", x$n_valid_pairs, "valid )\n")
  cat("\n")

  cat("TV Distance Summary:\n")
  cat("  Mean:", sprintf("%.6f", x$mean_tv), "\n")
  cat("  Min: ", sprintf("%.6f", x$min_tv), "\n")
  cat("  Max: ", sprintf("%.6f", x$max_tv), "\n")
  cat("  SD:  ", sprintf("%.6f", x$sd_tv), "\n\n")

  cat("Utility score:", sprintf("%.4f", x$utility_score), "\n\n")

  cat("Pairwise Details:\n")
  print(x$pairwise, row.names = FALSE)
  cat("\n")

  cat("Sample Sizes: X =", x$n_X, ", Y =", x$n_Y, "\n")

  invisible(x)
}


#' Plot method for contingency_fidelity objects
#'
#' @param x an object of class \code{"contingency_fidelity"}
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = heatmap of pairwise TV distances
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2
#'   labs theme_minimal theme element_text
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.contingency_fidelity <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 1)
  show[which] <- TRUE

  if (show[1]) {
    # Build a symmetric matrix of TV distances
    vars <- x$vars
    n_v <- length(vars)
    mat <- matrix(NA_real_, nrow = n_v, ncol = n_v,
                  dimnames = list(vars, vars))

    for (i in seq_len(nrow(x$pairwise))) {
      v1 <- x$pairwise$var1[i]
      v2 <- x$pairwise$var2[i]
      val <- x$pairwise$tv_distance[i]
      if (!is.na(val)) {
        mat[v1, v2] <- val
        mat[v2, v1] <- val
      }
    }

    # Use reshape2::melt() explicitly to avoid data.table conflict
    df <- reshape2::melt(mat, varnames = c("Var1", "Var2"), value.name = "TV")

    valid_tv <- x$pairwise$tv_distance[!is.na(x$pairwise$tv_distance)]
    mid <- if (length(valid_tv) > 0) median(valid_tv) else 0

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Var1, y = Var2,
                                           fill = TV)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(TV), "",
                                                      sprintf("%.4f", TV))),
                         size = 3.5) +
      ggplot2::scale_fill_gradient2(
        low = "steelblue", mid = "white", high = "firebrick",
        midpoint = mid,
        name = "TV Distance"
      ) +
      ggplot2::labs(
        title = "Pairwise Contingency Table TV Distances",
        subtitle = sprintf("Mean TV = %.4f | Utility = %.4f",
                           x$mean_tv, x$utility_score),
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

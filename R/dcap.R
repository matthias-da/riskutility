#' Correct Attribution Probability (CAP/DCAP)
#'
#' Computes the Correct Attribution Probability for synthetic data disclosure risk.
#' Measures the probability that an adversary can correctly infer a sensitive
#' target variable using known quasi-identifier (key) attributes.
#'
#' @param X data frame of original data, or a \code{\link{synth_pair}} object
#' @param Y data frame of synthetic/anonymized data (not needed if X is a synth_pair)
#' @param key_vars character vector of quasi-identifier variable names
#' @param target_var character, name of the sensitive target variable
#' @param method character, matching method: "exact" or "gower" (default: "exact")
#' @param gower_threshold numeric, maximum Gower distance for a match (default: 0.1)
#' @param cont_bins integer, number of bins for continuous target variables (default: 10)
#' @param na.rm logical, remove records with NA in key or target (default: TRUE)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "dcap" containing:
#' \itemize{
#'   \item cap_scores: per-record CAP values
#'   \item dcap: mean CAP (overall disclosure risk)
#'   \item dcap_median: median CAP
#'   \item n_matched: records with at least one synthetic match
#'   \item n_unmatched: records with no matches
#'   \item baseline: expected CAP under random guessing
#'   \item key_vars, target_var, method: input parameters
#' }
#'
#' @details
#' The Correct Attribution Probability (CAP) measures attribute disclosure risk
#' in synthetic data. For each original record, it finds matching synthetic
#' records based on key attributes and computes the probability that the
#' correct target value appears among the matches.
#'
#' \deqn{CAP_i = \frac{|\{j \in Y : keys_j = keys_i \land target_j = target_i\}|}{|\{j \in Y : keys_j = keys_i\}|}}
#'
#' Higher CAP values indicate higher disclosure risk. A CAP of 1 means the
#' intruder can perfectly infer the target; a CAP equal to baseline means
#' no information leakage beyond random guessing.
#'
#' @section Baseline computation:
#' The baseline is computed as the maximum target category frequency in the
#' synthetic data (Y), representing the best guess an intruder could make
#' without any key information.
#'
#' \strong{Note:} This differs from synthpop's baseline computation, which uses
#' the original data (X). The riskutility approach is more conservative: it
#' measures what an intruder learns from synthetic data specifically, rather
#' than from the original distribution. Both approaches are valid but answer
#' slightly different questions.
#'
#' @section Gower threshold selection:
#' When \code{method = "gower"}, the \code{gower_threshold} parameter determines
#' how close records must be to be considered "matching." Guidelines:
#' \itemize{
#'   \item \strong{0.05}: Very strict matching, few matches, may miss disclosure
#'   \item \strong{0.1} (default): Balanced threshold, suitable for most cases
#'   \item \strong{0.2}: More permissive, catches near-matches
#'   \item \strong{0.3+}: Very permissive, may inflate risk estimates
#' }
#' Consider the scale and nature of your quasi-identifiers when choosing.
#' For datasets with many categorical variables, lower thresholds work well.
#' For continuous-heavy data, slightly higher thresholds may be appropriate.
#'
#' @section synth_pair method:
#' If \code{X} is a \code{\link{synth_pair}} object, the function extracts
#' original, synthetic, key_vars, and target_var from the object. The synth_pair
#' must have \code{key_vars} and \code{target_var} set.
#'
#' @references
#' Taub, J., Elliot, M., Pampaka, M., & Smith, D. (2018). Differential Correct
#' Attribution Probability for Synthetic Data: An Exploration.
#' \emph{Privacy in Statistical Databases}, 122-137.
#'
#' @seealso \code{\link{synth_pair}} for creating a comparison pair,
#'   \code{\link{tcap}} for per-record CAP, \code{\link{disco}} for disclosive records
#'
#' @author Matthias Templ
#' @export
#' @importFrom VIM gowerD
#' @importFrom data.table data.table merge.data.table setorder .N
#' @importFrom stats complete.cases quantile median sd
#' @importFrom graphics hist abline par legend
#' @examples
#' # Create example data
#' set.seed(123)
#' X <- data.frame(
#'   age = sample(20:60, 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 100, replace = TRUE),
#'   income = sample(c("low", "medium", "high"), 100, replace = TRUE)
#' )
#' # Synthetic version with some noise
#' Y <- X
#' Y$income <- sample(Y$income)  # Shuffle target
#'
#' # Compute DCAP (traditional API)
#' result <- dcap(X, Y,
#'                key_vars = c("age", "gender", "region"),
#'                target_var = "income")
#' print(result)
#' summary(result)
#'
#' # Using synth_pair (container API)
#' pair <- synth_pair(X, Y,
#'                    key_vars = c("age", "gender", "region"),
#'                    target_var = "income")
#' result2 <- dcap(pair)
dcap <- function(X, ...) {
  UseMethod("dcap")
}

#' @rdname dcap
#' @export
dcap.synth_pair <- function(X, ...) {
  if (is.null(X$key_vars)) {
    stop("synth_pair must have 'key_vars' set for dcap()")
  }
  if (is.null(X$target_var)) {
    stop("synth_pair must have 'target_var' set for dcap()")
  }

  dcap.default(
    X = X$original,
    Y = X$synthetic,
    key_vars = X$key_vars,
    target_var = X$target_var,
    ...
  )
}

#' @rdname dcap
#' @export
dcap.default <- function(X, Y,
                         key_vars,
                         target_var,
                         method = c("exact", "gower"),
                         gower_threshold = 0.1,
                         cont_bins = 10,
                         na.rm = TRUE,
                         ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

  method <- match.arg(method)

  # Check variables exist
  all_vars <- c(key_vars, target_var)
  missing_X <- setdiff(all_vars, names(X))
  missing_Y <- setdiff(all_vars, names(Y))
  if (length(missing_X) > 0) {
    stop(paste("Variables missing in X:", paste(missing_X, collapse = ", ")))
  }
  if (length(missing_Y) > 0) {
    stop(paste("Variables missing in Y:", paste(missing_Y, collapse = ", ")))
  }

  # Check variable types match
  for (var in all_vars) {
    if (!identical(class(X[[var]]), class(Y[[var]]))) {
      stop(paste("Variable", var, "has different class in X and Y."))
    }
  }

  # Handle missing values
  if (na.rm) {
    complete_X <- complete.cases(X[, all_vars])
    complete_Y <- complete.cases(Y[, all_vars])
    X <- X[complete_X, , drop = FALSE]
    Y <- Y[complete_Y, , drop = FALSE]
  }

  if (nrow(X) == 0) stop("No complete cases in X after removing NAs.")
  if (nrow(Y) == 0) stop("No complete cases in Y after removing NAs.")

  # Prepare target variable (bin continuous)
  target_X <- X[[target_var]]
  target_Y <- Y[[target_var]]

  if (is.numeric(target_X) && !is.factor(target_X)) {
    # Bin continuous target
    breaks <- quantile(c(target_X, target_Y),
                       probs = seq(0, 1, length.out = cont_bins + 1),
                       na.rm = TRUE)
    breaks <- unique(breaks)
    target_X <- cut(target_X, breaks = breaks, include.lowest = TRUE)
    target_Y <- cut(target_Y, breaks = breaks, include.lowest = TRUE)
  }

  # Compute baseline (random guessing probability)
  target_freq <- table(target_Y) / length(target_Y)
  baseline <- max(target_freq)  # Best guess without keys

  # Initialize CAP scores
  n <- nrow(X)
  cap_scores <- numeric(n)
  n_matches <- integer(n)

  if (method == "exact") {
    # Vectorized exact matching using data.table for O(n + m) performance
    # Create key signature using fast vectorized paste
    make_key <- function(df, vars) {
      do.call(paste, c(df[, vars, drop = FALSE], sep = "|"))
    }

    keys_X <- make_key(X, key_vars)
    keys_Y <- make_key(Y, key_vars)

    # Build lookup tables using data.table for fast grouping
    dt_Y <- data.table::data.table(
      key_sig = keys_Y,
      target = target_Y,
      idx = seq_along(keys_Y)
    )

    # Pre-compute: for each key in Y, count total and count per target value
    key_counts <- dt_Y[, .(n_total = .N), by = key_sig]
    key_target_counts <- dt_Y[, .(n_correct = .N), by = .(key_sig, target)]

    # Create lookup for original data
    dt_X <- data.table::data.table(
      orig_idx = seq_len(n),
      key_sig = keys_X,
      target = target_X
    )

    # Join to get match counts
    dt_X <- data.table::merge.data.table(dt_X, key_counts, by = "key_sig", all.x = TRUE)
    dt_X <- data.table::merge.data.table(dt_X, key_target_counts,
                                          by = c("key_sig", "target"), all.x = TRUE)

    # Fill NAs with 0
    dt_X[is.na(n_total), n_total := 0L]
    dt_X[is.na(n_correct), n_correct := 0L]

    # Restore original order
    data.table::setorder(dt_X, orig_idx)

    # Extract results
    n_matches <- dt_X$n_total
    cap_scores <- ifelse(n_matches == 0, NA_real_, dt_X$n_correct / n_matches)

  } else if (method == "gower") {
    # Gower distance-based matching (vectorized)
    # Compute pairwise Gower distances
    gd <- VIM::gowerD(X[, key_vars, drop = FALSE],
                      Y[, key_vars, drop = FALSE])

    # Vectorized: create match matrix (TRUE where distance <= threshold)
    match_matrix <- gd <= gower_threshold

    # Count matches per original record (row sums)
    n_matches <- rowSums(match_matrix)

    # Vectorized: create target match matrix
    # target_match[i,j] = TRUE if target_X[i] == target_Y[j]
    target_match <- outer(target_X, target_Y, "==")

    # Count correct matches: sum where both match_matrix AND target_match are TRUE
    n_correct <- rowSums(match_matrix & target_match)

    # Compute CAP scores
    cap_scores <- ifelse(n_matches == 0, NA_real_, n_correct / n_matches)
  }

  # Aggregate results
  n_matched <- sum(!is.na(cap_scores))
  n_unmatched <- sum(is.na(cap_scores))

  results <- list(
    cap_scores = cap_scores,
    n_matches = n_matches,
    dcap = mean(cap_scores, na.rm = TRUE),
    dcap_median = median(cap_scores, na.rm = TRUE),
    n_matched = n_matched,
    n_unmatched = n_unmatched,
    n_total = n,
    baseline = baseline,
    key_vars = key_vars,
    target_var = target_var,
    method = method,
    gower_threshold = if (method == "gower") gower_threshold else NULL
  )

  class(results) <- "dcap"
  return(results)
}

#' Print method for dcap objects
#'
#' @param x an object of class "dcap"
#' @param ... additional arguments passed to the print method
#' @export
print.dcap <- function(x, ...) {
  cat("Correct Attribution Probability (DCAP) Analysis\n")
  cat("================================================\n")
  cat("Method:", x$method, "\n")
  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Target variable:", x$target_var, "\n\n")
  cat("DCAP (mean CAP):", round(x$dcap, 4), "\n")
  cat("Median CAP:", round(x$dcap_median, 4), "\n")
  cat("Baseline (random):", round(x$baseline, 4), "\n")
  cat("Records matched:", x$n_matched, "/", x$n_total,
      sprintf("(%.1f%%)", 100 * x$n_matched / x$n_total), "\n")
  invisible(x)
}

#' Summary method for dcap objects
#'
#' @param object an object of class "dcap"
#' @param ... additional arguments passed to the summary method
#' @export
summary.dcap <- function(object, ...) {
  cap <- object$cap_scores[!is.na(object$cap_scores)]
  summ <- list(
    dcap = object$dcap,
    dcap_median = object$dcap_median,
    baseline = object$baseline,
    risk_ratio = object$dcap / object$baseline,
    cap_quantiles = quantile(cap, probs = c(0, 0.25, 0.5, 0.75, 1)),
    cap_sd = sd(cap),
    n_matched = object$n_matched,
    n_unmatched = object$n_unmatched,
    mean_matches = mean(object$n_matches[object$n_matches > 0]),
    method = object$method,
    key_vars = object$key_vars,
    target_var = object$target_var
  )
  class(summ) <- "summary.dcap"
  return(summ)
}

#' Print method for summary.dcap objects
#'
#' @param x an object of class "summary.dcap"
#' @param ... additional arguments passed to the print method
#' @export
print.summary.dcap <- function(x, ...) {
  cat("Summary: Correct Attribution Probability (DCAP)\n")
  cat("================================================\n")
  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Target variable:", x$target_var, "\n")
  cat("Method:", x$method, "\n\n")

  cat("Risk Assessment:\n")
  cat("  DCAP (mean):", round(x$dcap, 4), "\n")
  cat("  Baseline:", round(x$baseline, 4), "\n")
  cat("  Risk ratio:", round(x$risk_ratio, 2),
      ifelse(x$risk_ratio > 1.5, " (elevated risk)",
             ifelse(x$risk_ratio > 1, " (slightly elevated)", " (low risk)")), "\n\n")

  cat("CAP Distribution:\n")
  print(round(x$cap_quantiles, 4))
  cat("  SD:", round(x$cap_sd, 4), "\n\n")

  cat("Matching Statistics:\n")
  cat("  Records matched:", x$n_matched, "\n")
  cat("  Records unmatched:", x$n_unmatched, "\n")
  cat("  Avg matches per record:", round(x$mean_matches, 1), "\n")
  invisible(x)
}

#' Plot method for dcap objects
#'
#' @param x an object of class "dcap"
#' @param y not used
#' @param ... additional arguments passed to the plot method
#' @param which which plot(s) to show: 1 for CAP histogram, 2 for matches histogram
#' @export
plot.dcap <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  cap <- x$cap_scores[!is.na(x$cap_scores)]

  if (sum(show) > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Histogram of CAP scores
    hist(cap, breaks = 20, main = "Distribution of CAP Scores",
         xlab = "CAP Score", col = "steelblue", border = "white", ...)
    abline(v = x$dcap, col = "red", lwd = 2, lty = 2)
    abline(v = x$baseline, col = "darkgreen", lwd = 2, lty = 3)
    legend("topright", legend = c("Mean DCAP", "Baseline"),
           col = c("red", "darkgreen"), lty = c(2, 3), lwd = 2)
  }

  if (show[2]) {
    # Matches histogram
    matches <- x$n_matches[x$n_matches > 0]
    hist(matches, breaks = 30, main = "Number of Matches per Record",
         xlab = "Number of Synthetic Matches", col = "coral", border = "white", ...)
    abline(v = mean(matches), col = "red", lwd = 2, lty = 2)
  }
}

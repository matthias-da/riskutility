#' Disclosive in Synthetic Correct Original (DiSCO)
#'
#' Identifies synthetic records that could potentially disclose information about
#' original records. A DiSCO record is one that matches an original record on
#' key variables AND has the same target value, making it potentially disclosive.
#'
#' @param X data frame of original data
#' @param Y data frame of synthetic data
#' @param key_vars character vector of quasi-identifier variable names
#' @param target_var character, name of the sensitive target variable
#' @param disclosure_type character, type of disclosure to measure:
#'   \itemize{
#'     \item \code{"potential"} (default): Counts any synthetic record that matches
#'       an original record on key+target. This is the broader measure.
#'     \item \code{"certain"}: Only counts disclosures where the key combination
#'       uniquely determines the target in the original data (i.e., all original
#'       records with that key have the same target value). This is compatible
#'       with synthpop's DiSCO measure.
#'   }
#' @param method character, matching method: "exact" or "gower" (default: "exact")
#' @param gower_threshold numeric, maximum Gower distance for a match when
#'   method="gower" (default: 0.1)
#' @param na.rm logical, remove records with NA in key or target (default: TRUE)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "disco" containing:
#' \itemize{
#'   \item disco_idx: indices of DiSCO records in Y (synthetic data)
#'   \item disco_records: the actual DiSCO records from Y
#'   \item n_disco: number of DiSCO records (based on disclosure_type)
#'   \item pct_disco: percentage of synthetic records that are DiSCO
#'   \item matched_original_idx: for each DiSCO record, which original record(s) it matches
#'   \item n_disco_potential: count using "potential" method (always computed)
#'   \item n_disco_certain: count using "certain" method (always computed)
#'   \item pct_original_disclosed: percentage of original records disclosed (for "certain" method)
#'   \item baseline_disco: expected DiSCO count under random target assignment
#'   \item key_vars, target_var, method, disclosure_type: input parameters
#' }
#'
#' @details
#' DiSCO (Disclosive in Synthetic Correct Original) is a measure developed for
#' the synthpop package to identify attribute disclosure risk. It identifies
#' synthetic records that:
#' \enumerate{
#'   \item Match an original record on all key variables
#'   \item Have the same target value as that original record
#' }
#'
#' Two disclosure types are available:
#'
#' \strong{Potential disclosure} (\code{disclosure_type = "potential"}):
#' Counts any synthetic record that matches an original on key+target.
#' This is a broader measure that identifies all records that \emph{could}
#' potentially leak information. Use this when you want to identify all
#' possible disclosure risks.
#'
#' \strong{Certain disclosure} (\code{disclosure_type = "certain"}):
#' Only counts disclosures where the key combination uniquely determines
#' the target in the original data. This means an intruder who matches
#' a synthetic record to the original would learn the target with certainty.
#' This is compatible with synthpop's DiSCO measure and is more conservative.
#'
#' The function always computes both metrics for comparison, but the primary
#' output (n_disco, pct_disco, disco_idx) reflects the chosen disclosure_type.
#'
#' Interpretation:
#' \itemize{
#'   \item High DiSCO count: Many synthetic records could leak original information
#'   \item DiSCO/Baseline ratio > 1: Synthetic data leaks more than expected by chance
#'   \item DiSCO = 0: No exact key+target matches (lowest risk scenario)
#'   \item pct_original_disclosed: For "certain" method, shows what fraction of
#'     original records have their target certainly exposed
#' }
#'
#' @seealso \code{\link{tcap}} for targeted correct attribution probability,
#'   \code{\link{weap}} for within equivalence class attribution probability,
#'   \code{\link{dcap}} for differential correct attribution probability
#'
#' @references
#' Raab, G.M., Nowok, B., & Dibben, C. (2021). Assessing, visualizing and
#' improving the utility of synthetic data. \emph{arXiv preprint arXiv:2109.12717}.
#'
#' @author Matthias Templ
#' @family attribution-risk
#' @export
#' @importFrom VIM gowerD
#' @importFrom stats complete.cases
#' @importFrom graphics barplot pie
#' @examples
#' # Create example data
#' set.seed(42)
#' X <- data.frame(
#'   age = sample(20:40, 50, replace = TRUE),
#'   gender = sample(c("M", "F"), 50, replace = TRUE),
#'   disease = sample(c("none", "A", "B"), 50, replace = TRUE,
#'                    prob = c(0.7, 0.2, 0.1))
#' )
#' # Synthetic version - some records will match by chance
#' Y <- data.frame(
#'   age = sample(20:40, 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   disease = sample(c("none", "A", "B"), 100, replace = TRUE,
#'                    prob = c(0.7, 0.2, 0.1))
#' )
#'
#' # Find DiSCO records using "potential" disclosure (default)
#' result_potential <- disco(X, Y,
#'                           key_vars = c("age", "gender"),
#'                           target_var = "disease",
#'                           disclosure_type = "potential")
#' print(result_potential)
#'
#' # Find DiSCO records using "certain" disclosure (synthpop-compatible)
#' result_certain <- disco(X, Y,
#'                         key_vars = c("age", "gender"),
#'                         target_var = "disease",
#'                         disclosure_type = "certain")
#' print(result_certain)
#'
#' # Compare both methods
#' cat("Potential DiSCO:", result_potential$n_disco_potential, "\n")
#' cat("Certain DiSCO:", result_potential$n_disco_certain, "\n")
#'
#' # View the disclosive records
#' if (result_potential$n_disco > 0) {
#'   head(result_potential$disco_records)
#' }
disco <- function(X, ...) {
  UseMethod("disco")
}

#' @rdname disco
#' @export
disco.synth_pair <- function(X, ...) {
  if (is.null(X$key_vars)) {
    stop("synth_pair must have 'key_vars' set for disco()")
  }
  if (is.null(X$target_var)) {
    stop("synth_pair must have 'target_var' set for disco()")
  }

  disco.default(
    X = X$original,
    Y = X$synthetic,
    key_vars = X$key_vars,
    target_var = X$target_var,
    ...
  )
}

#' @rdname disco
#' @export
disco.default <- function(X, Y,
                          key_vars,
                          target_var,
                          disclosure_type = c("potential", "certain"),
                          method = c("exact", "gower"),
                          gower_threshold = 0.1,
                          na.rm = TRUE,
                          ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

  disclosure_type <- match.arg(disclosure_type)
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
    complete_X <- complete.cases(X[, all_vars, drop = FALSE])
    complete_Y <- complete.cases(Y[, all_vars, drop = FALSE])
    X <- X[complete_X, , drop = FALSE]
    Y <- Y[complete_Y, , drop = FALSE]
  }

  if (nrow(X) == 0) stop("No complete cases in X after removing NAs.")
  if (nrow(Y) == 0) stop("No complete cases in Y after removing NAs.")

  n_X <- nrow(X)
  n_Y <- nrow(Y)
  target_X <- X[[target_var]]
  target_Y <- Y[[target_var]]

  # Create key signatures
  make_key <- function(df, vars) {
    apply(df[, vars, drop = FALSE], 1, paste, collapse = "|")
  }
  keys_X <- make_key(X, key_vars)
  keys_Y <- make_key(Y, key_vars)

  # For "certain" disclosure: identify keys that uniquely determine target
  # A key uniquely determines target if all original records with that key

  # have the same target value
  key_target_unique <- tapply(target_X, keys_X, function(t) length(unique(t)) == 1)
  unique_keys <- names(key_target_unique)[key_target_unique]

  # Track DiSCO records for both methods
  is_disco_potential <- logical(n_Y)
  is_disco_certain <- logical(n_Y)
  matched_original <- vector("list", n_Y)

  if (method == "exact") {
    # Create combined key+target signatures
    key_target_X <- paste(keys_X, target_X, sep = "||")
    key_target_Y <- paste(keys_Y, target_Y, sep = "||")

    # Find matches
    for (j in seq_len(n_Y)) {
      # Which original records have same key AND target?
      matches <- which(key_target_X == key_target_Y[j])
      if (length(matches) > 0) {
        is_disco_potential[j] <- TRUE
        matched_original[[j]] <- matches

        # For "certain" disclosure: check if this key uniquely determines target
        if (keys_Y[j] %in% unique_keys) {
          is_disco_certain[j] <- TRUE
        }
      }
    }

  } else if (method == "gower") {
    # Gower distance-based matching
    gd <- VIM::gowerD(Y[, key_vars, drop = FALSE],
                      X[, key_vars, drop = FALSE])

    for (j in seq_len(n_Y)) {
      # Find original records within threshold
      distances <- gd[j, ]
      key_matches <- which(distances <= gower_threshold)

      if (length(key_matches) > 0) {
        # Check if target also matches
        target_matches <- key_matches[target_X[key_matches] == target_Y[j]]
        if (length(target_matches) > 0) {
          is_disco_potential[j] <- TRUE
          matched_original[[j]] <- target_matches

          # For "certain": check if all matched records have same target
          all_targets_same <- length(unique(target_X[key_matches])) == 1
          if (all_targets_same) {
            is_disco_certain[j] <- TRUE
          }
        }
      }
    }
  }

  # Extract DiSCO counts for both methods
  disco_idx_potential <- which(is_disco_potential)
  disco_idx_certain <- which(is_disco_certain)

  n_disco_potential <- length(disco_idx_potential)
  n_disco_certain <- length(disco_idx_certain)

  # Select primary results based on disclosure_type
  if (disclosure_type == "potential") {
    disco_idx <- disco_idx_potential
    n_disco <- n_disco_potential
    is_disco <- is_disco_potential
  } else {
    disco_idx <- disco_idx_certain
    n_disco <- n_disco_certain
    is_disco <- is_disco_certain
  }

  pct_disco <- 100 * n_disco / n_Y

  # Calculate percentage of original records disclosed (for "certain" method)
  # These are original records in key groups where key uniquely determines target
  # AND there's a matching synthetic record
  matched_originals_certain <- unique(unlist(matched_original[disco_idx_certain]))
  n_originals_disclosed_certain <- length(matched_originals_certain)
  pct_original_disclosed <- 100 * n_originals_disclosed_certain / n_X

  # Count key matches
  n_key_matches <- sum(keys_Y %in% keys_X)

  # Expected proportion of correct targets given a key match
  target_freq <- table(target_X) / n_X
  expected_target_match <- sum(target_freq^2)  # Prob of random match

  baseline_disco <- n_key_matches * expected_target_match

  results <- list(
    disco_idx = disco_idx,
    disco_records = if (n_disco > 0) Y[disco_idx, , drop = FALSE] else Y[0, , drop = FALSE],
    n_disco = n_disco,
    pct_disco = pct_disco,
    n_synthetic = n_Y,
    n_original = n_X,
    matched_original_idx = matched_original[disco_idx],
    # Both methods always computed
    n_disco_potential = n_disco_potential,
    pct_disco_potential = 100 * n_disco_potential / n_Y,
    n_disco_certain = n_disco_certain,
    pct_disco_certain = 100 * n_disco_certain / n_Y,
    # Original records disclosed (certain method)
    n_originals_disclosed = n_originals_disclosed_certain,
    pct_original_disclosed = pct_original_disclosed,
    # Other metrics
    n_key_matches = n_key_matches,
    baseline_disco = baseline_disco,
    disco_ratio = if (baseline_disco > 0) n_disco / baseline_disco else NA,
    key_vars = key_vars,
    target_var = target_var,
    method = method,
    disclosure_type = disclosure_type,
    gower_threshold = if (method == "gower") gower_threshold else NULL
  )

  class(results) <- "disco"
  return(results)
}

#' Print method for disco objects
#'
#' @param x an object of class "disco"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.disco <- function(x, ...) {
  cat("Disclosive in Synthetic Correct Original (DiSCO)\n")
  cat("================================================\n")
  cat("Matching method:", x$method, "\n")
  cat("Disclosure type:", x$disclosure_type, "\n")
  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Target variable:", x$target_var, "\n\n")

  cat("Data Summary:\n")
  cat("  Original records:", x$n_original, "\n")
  cat("  Synthetic records:", x$n_synthetic, "\n")
  cat("  Synthetic with key match:", x$n_key_matches,
      sprintf("(%.1f%%)", 100 * x$n_key_matches / x$n_synthetic), "\n\n")

  cat("DiSCO Results (", x$disclosure_type, " disclosure):\n", sep = "")
  cat("  DiSCO records:", x$n_disco,
      sprintf("(%.1f%% of synthetic)", x$pct_disco), "\n")
  cat("  Expected by chance:", round(x$baseline_disco, 1), "\n")

  if (!is.na(x$disco_ratio)) {
    cat("  DiSCO/Baseline ratio:", round(x$disco_ratio, 2))
    if (x$disco_ratio > 1.5) {
      cat(" (elevated risk)\n")
    } else if (x$disco_ratio > 1.0) {
      cat(" (slightly elevated)\n")
    } else {
      cat(" (low risk)\n")
    }
  }

  cat("\nComparison of Both Methods:\n")
  cat("  Potential DiSCO:", x$n_disco_potential,
      sprintf("(%.1f%% of synthetic)", x$pct_disco_potential), "\n")
  cat("  Certain DiSCO:", x$n_disco_certain,
      sprintf("(%.1f%% of synthetic)", x$pct_disco_certain), "\n")
  cat("  Original records disclosed:", x$n_originals_disclosed,
      sprintf("(%.1f%% of original)", x$pct_original_disclosed), "\n")

  invisible(x)
}

#' Summary method for disco objects
#'
#' @param object an object of class "disco"
#' @param ... additional arguments (ignored)
#' @return A list of summary statistics for the corresponding object.
#' @export
summary.disco <- function(object, ...) {
  # Count how many original records are matched by DiSCO records
  all_matched <- unique(unlist(object$matched_original_idx))
  n_originals_exposed <- length(all_matched)

  # Distribution of matches per DiSCO record
  matches_per_disco <- sapply(object$matched_original_idx, length)

  summ <- list(
    n_disco = object$n_disco,
    pct_disco = object$pct_disco,
    n_synthetic = object$n_synthetic,
    n_original = object$n_original,
    n_key_matches = object$n_key_matches,
    pct_key_matches = 100 * object$n_key_matches / object$n_synthetic,
    baseline_disco = object$baseline_disco,
    disco_ratio = object$disco_ratio,
    # Both methods
    n_disco_potential = object$n_disco_potential,
    pct_disco_potential = object$pct_disco_potential,
    n_disco_certain = object$n_disco_certain,
    pct_disco_certain = object$pct_disco_certain,
    n_originals_disclosed = object$n_originals_disclosed,
    pct_original_disclosed = object$pct_original_disclosed,
    # Impact analysis
    n_originals_exposed = n_originals_exposed,
    pct_originals_exposed = 100 * n_originals_exposed / object$n_original,
    matches_per_disco_mean = if (length(matches_per_disco) > 0) mean(matches_per_disco) else NA,
    matches_per_disco_max = if (length(matches_per_disco) > 0) max(matches_per_disco) else NA,
    key_vars = object$key_vars,
    target_var = object$target_var,
    method = object$method,
    disclosure_type = object$disclosure_type
  )

  class(summ) <- "summary.disco"
  return(summ)
}

#' Print method for summary.disco objects
#'
#' @param x an object of class "summary.disco"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.disco <- function(x, ...) {
  cat("Summary: Disclosive in Synthetic Correct Original (DiSCO)\n")
  cat("=========================================================\n")
  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Target variable:", x$target_var, "\n")
  cat("Matching method:", x$method, "\n")
  cat("Disclosure type:", x$disclosure_type, "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original records:", x$n_original, "\n")
  cat("  Synthetic records:", x$n_synthetic, "\n\n")

  cat("Key Matching:\n")
  cat("  Synthetic records with key match:", x$n_key_matches,
      sprintf("(%.1f%%)", x$pct_key_matches), "\n\n")

  cat("DiSCO Analysis (primary: ", x$disclosure_type, "):\n", sep = "")
  cat("  DiSCO records found:", x$n_disco,
      sprintf("(%.1f%% of synthetic)", x$pct_disco), "\n")
  cat("  Expected by chance:", round(x$baseline_disco, 1), "\n")
  if (!is.na(x$disco_ratio)) {
    cat("  DiSCO/Baseline ratio:", round(x$disco_ratio, 2), "\n")
  }
  cat("\n")

  cat("Comparison of Both Methods:\n")
  cat("  Potential DiSCO:", x$n_disco_potential,
      sprintf("(%.1f%% of synthetic)", x$pct_disco_potential), "\n")
  cat("  Certain DiSCO:", x$n_disco_certain,
      sprintf("(%.1f%% of synthetic)", x$pct_disco_certain), "\n")
  cat("  Original records disclosed (certain):", x$n_originals_disclosed,
      sprintf("(%.1f%% of original)", x$pct_original_disclosed), "\n\n")

  cat("Impact on Original Data:\n")
  cat("  Original records exposed:", x$n_originals_exposed,
      sprintf("(%.1f%% of original)", x$pct_originals_exposed), "\n")
  if (!is.na(x$matches_per_disco_mean)) {
    cat("  Avg original matches per DiSCO:", round(x$matches_per_disco_mean, 1), "\n")
    cat("  Max original matches per DiSCO:", x$matches_per_disco_max, "\n")
  }

  invisible(x)
}

#' Plot method for disco objects
#'
#' @param x an object of class "disco"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot(s) to show: 1 = DiSCO vs non-DiSCO bar chart,
#'   2 = comparison with baseline, 3 = comparison of both methods
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.disco <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 3)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    if (n_plots == 2) {
      op <- par(mfrow = c(1, 2))
    } else {
      op <- par(mfrow = c(1, 3))
    }
    on.exit(par(op))
  }

  if (show[1]) {
    # Bar chart of DiSCO vs non-DiSCO
    counts <- c(x$n_synthetic - x$n_disco, x$n_disco)
    names(counts) <- c("Non-DiSCO", "DiSCO")
    barplot(counts,
            main = paste("Synthetic Record Classification\n(",
                         x$disclosure_type, " disclosure)", sep = ""),
            ylab = "Number of Records",
            col = c("forestgreen", "firebrick"), ...)
  }

  if (show[2]) {
    # Observed vs expected DiSCO
    counts <- c(x$baseline_disco, x$n_disco)
    names(counts) <- c("Expected\n(baseline)", "Observed\n(DiSCO)")
    barplot(counts,
            main = "DiSCO: Observed vs Expected",
            ylab = "Count",
            col = c("steelblue", "coral"), ...)
    abline(h = x$baseline_disco, col = "steelblue", lwd = 2, lty = 2)
  }

  if (show[3]) {
    # Comparison of both methods
    counts <- c(x$n_disco_potential, x$n_disco_certain)
    names(counts) <- c("Potential", "Certain")
    barplot(counts,
            main = "DiSCO: Potential vs Certain",
            ylab = "Number of DiSCO Records",
            col = c("darkorange", "darkgreen"), ...)
  }
}

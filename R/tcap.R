#' Targeted Correct Attribution Probability (TCAP)
#'
#' Computes the Targeted Correct Attribution Probability for assessing attribute
#' disclosure risk in synthetic data. TCAP measures, for each original record,
#' the probability that an intruder correctly infers a sensitive target variable
#' by matching on quasi-identifier (key) attributes.
#'
#' @param X data frame of original data
#' @param Y data frame of synthetic data
#' @param key_vars character vector of quasi-identifier variable names
#' @param target_var character, name of the sensitive target variable
#' @param method character, matching method: "exact" or "gower" (default: "exact")
#' @param gower_threshold numeric, maximum Gower distance for a match when
#'   method="gower" (default: 0.1)
#' @param cont_bins integer, number of quantile bins used to discretise a
#'   continuous target variable (default: 10). Matching \code{dcap()}, a
#'   numeric target is binned before attribution is assessed; without
#'   binning an exact match on a continuous value essentially never occurs
#'   and every TCAP score collapses to zero.
#' @param na.rm logical, remove records with NA in key or target (default: TRUE)
#' @param kind character, which TCAP variant \code{print()} and
#'   \code{summary()} highlight as the selected headline: \code{"certain"}
#'   (default), \code{"matched"}, or \code{"conditional"}. All three variants
#'   are always computed and returned regardless of \code{kind}; see Details
#'   for their definitions and the mapping to synthpop versions.
#' @param ... additional arguments passed to methods
#'
#' @return An object of class "tcap" containing:
#' \itemize{
#'   \item tcap_scores: numeric vector of per-record CAP values
#'   \item tcap_mean: mean CAP across matched records. synthpop's DCAP is the
#'     same sum divided by the number of all original records (unmatched
#'     records counting as zero), so the two coincide exactly when every
#'     original record has a key match; see \code{\link{dcap}} for the
#'     synthpop-DCAP estimand itself.
#'   \item tcap_max: maximum CAP (worst-case individual risk)
#'   \item tcap_median: median CAP
#'   \item tcap_certain: percentage of matched records at "certain" disclosure
#'     risk: the key uniquely determines the target in both the original and
#'     the synthetic data, and the attribution is correct
#'   \item tcap_matched: percentage of matched records that are disclosive in
#'     the synthetic data and correctly attributed (CAP = 1). Reproduces
#'     synthpop's TCAP as defined up to synthpop 1.9-2.
#'   \item tcap_conditional: the same numerator divided by the number of
#'     records matched to a disclosive synthetic key class only
#'     (Little et al., 2025). This is the TCAP definition synthpop adopts
#'     from version 1.9-3. Following the reference implementation, 0 when no
#'     record is matched to a disclosive class.
#'   \item n_certain: count of records at certain disclosure risk
#'   \item n_disclosive_correct: count of matched records with CAP = 1 (the
#'     shared numerator of tcap_matched and tcap_conditional)
#'   \item n_disclosive: count of records matched to a disclosive synthetic
#'     key class (the denominator of tcap_conditional)
#'   \item is_certain, is_disclosive: per-record logical flags
#'   \item kind: the selected headline variant
#'   \item n_matches: number of synthetic matches per original record
#'   \item n_matched: number of original records with at least one match
#'   \item n_unmatched: number of original records with no matches
#'   \item baseline: baseline probability (modal target frequency in Y)
#'   \item key_vars, target_var, method: input parameters
#' }
#'
#' @details
#' This function computes Correct Attribution Probability (CAP) metrics developed
#' by Elliot (2014) and Taub et al. (2018) to assess attribute disclosure risk
#' in synthetic data.
#'
#' For each original record i with key values k_i and target value t_i:
#' \deqn{CAP_i = \frac{|\{j \in Y : keys_j = k_i \land target_j = t_i\}|}{|\{j \in Y : keys_j = k_i\}|}}
#'
#' The literature has never used a single denominator for the "targeted" CAP,
#' so \code{tcap()} computes the common variants side by side. Call a
#' synthetic key class \emph{disclosive} when all its synthetic records share
#' one target value; a matched original record has CAP = 1 exactly when its
#' class is disclosive and the attribution is correct.
#' \itemize{
#'   \item \strong{tcap_certain} (default headline): percentage of matched
#'     records whose key uniquely determines the target in BOTH the original
#'     AND the synthetic data, with correct attribution. The strictest
#'     variant: the intruder's inference is unambiguous in both datasets.
#'   \item \strong{tcap_matched}: percentage of matched records with CAP = 1.
#'     Reproduces synthpop's TCAP as defined up to synthpop 1.9-2 (DiSCO-type
#'     numerator over all original records whose key occurs in the synthetic
#'     data).
#'   \item \strong{tcap_conditional}: the same numerator divided by the
#'     number of original records matched to a disclosive synthetic key
#'     class, following Little et al. (2025); synthpop adopts this
#'     definition from version 1.9-3 (Raab et al., 2024). Records without
#'     such a match are excluded, and the value is 0 when no record is
#'     matched to a disclosive class (as in the reference implementation).
#' }
#'
#' On any data the variants order as
#' \code{tcap_certain <= tcap_matched <= tcap_conditional}: from
#' \code{tcap_matched} to \code{tcap_certain} the numerator shrinks (records
#' whose key is ambiguous in the original data drop out), and from
#' \code{tcap_matched} to \code{tcap_conditional} the denominator shrinks.
#'
#' For \code{method = "gower"}, equivalence classes are replaced by fuzzy
#' match sets: \code{is_certain} simplifies to CAP = 1 (so
#' \code{tcap_certain} equals \code{tcap_matched}), and a match set counts
#' as disclosive when all matched synthetic records share one target value.
#'
#' In addition, \strong{tcap_mean} gives the mean CAP across matched records
#' (the standard aggregate disclosure-risk measure) and \strong{tcap_max}
#' the maximum CAP (worst-case individual risk): a value of 1 means at least
#' one record has all its synthetic matches with the correct target value.
#'
#' Interpretation:
#' \itemize{
#'   \item CAP = 0: No synthetic matches have the correct target (low risk)
#'   \item CAP = 1: All synthetic matches have the correct target (high risk)
#'   \item CAP close to baseline: No additional risk beyond random guessing
#' }
#'
#' @section Baseline computation:
#' The baseline is computed as the maximum target category frequency in the
#' synthetic data, representing the best guess without key information.
#' \strong{Note:} This differs from synthpop's approach which uses the original data.
#' See \code{\link{dcap}} for detailed discussion.
#'
#' @section Gower threshold selection:
#' When \code{method = "gower"}, records within \code{gower_threshold} distance
#' are considered matches. The default of 0.1 works well for most datasets.
#' See \code{\link{dcap}} for guidelines on choosing this value.
#'
#' @seealso \code{\link{dcap}} for the differential measure comparing to baseline,
#'   \code{\link{weap}} for within equivalence class attribution probability
#'
#' @references
#' Elliot, M. (2014). Final Report on the Disclosure Risk Associated with the
#' Synthetic Data Produced by the SYLLS Team. Report 2015-2.
#'
#' Taub, J., Elliot, M., Pampaka, M., & Smith, D. (2018). Differential Correct
#' Attribution Probability for Synthetic Data: An Exploration.
#' \emph{Privacy in Statistical Databases}, 122-137.
#'
#' Little, C., Allmendinger, R., & Elliot, M. (2025). Synthetic Census
#' Microdata Generation: A Comparative Study of Synthesis Methods Examining
#' the Trade-Off Between Disclosure Risk and Utility.
#' \emph{Journal of Official Statistics}, 41(1), 255-308.
#' \doi{10.1177/0282423X241266523}
#'
#' Raab, G. M., Nowok, B., & Dibben, C. (2024). Practical Privacy Metrics for
#' Synthetic Data. arXiv:2406.16826.
#'
#' @author Matthias Templ
#' @family attribution-risk
#' @export
#' @importFrom VIM gowerD
#' @importFrom data.table data.table merge.data.table setorder .N
#' @importFrom stats complete.cases quantile median sd
#' @importFrom graphics hist abline legend barplot
#' @examples
#' # Create example data
#' set.seed(42)
#' X <- data.frame(
#'   age = sample(20:60, 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 100, replace = TRUE),
#'   disease = sample(c("none", "A", "B"), 100, replace = TRUE,
#'                    prob = c(0.7, 0.2, 0.1))
#' )
#' # Synthetic version
#' Y <- X
#' Y$disease <- sample(Y$disease)  # Shuffle target
#'
#' # Compute TCAP
#' result <- tcap(X, Y,
#'                key_vars = c("age", "gender", "region"),
#'                target_var = "disease")
#' print(result)
#' summary(result)
#' plot(result)
#'
#' # The three TCAP variants (always returned; `kind` picks the headline
#' # reported by print()/summary()); certain <= matched <= conditional
#' c(certain = result$tcap_certain,
#'   matched = result$tcap_matched,          # synthpop <= 1.9-2 TCAP
#'   conditional = result$tcap_conditional)  # Little et al. 2025; synthpop >= 1.9-3
#'
#' \donttest{
#' # Example 2: Using SD2011 dataset with synthpop-generated synthetic data
#' # This follows the synthpop::disclosure example pattern
#' # Requires: synthpop package
#'
#' if (requireNamespace("synthpop", quietly = TRUE)) {
#'
#' # Load SD2011 data and select variables (as in ?synthpop::disclosure)
#' ods <- synthpop::SD2011[, c("sex", "age", "edu", "marital", "income")]
#'
#' # Convert income to categorical (7 groups), handling special NA code
#' odsF <- synthpop::numtocat.syn(ods, numtocat = "income",
#'                                 catgroups = 7,
#'                                 cont.na = list(income = -8))
#' original_data <- odsF$data
#'
#' # Generate synthetic data using CART method
#' synth_obj <- synthpop::syn(original_data,
#'                            method = "ctree",
#'                            seed = 75,
#'                            m = 1,
#'                            k = 1000,
#'                            print.flag = FALSE)
#' synthetic_data <- synth_obj$syn
#'
#' # Compute TCAP: Can we infer income from demographic keys?
#' tcap_result <- tcap(original_data, synthetic_data,
#'                     key_vars = c("sex", "age", "edu", "marital"),
#'                     target_var = "income")
#'
#' # View results
#' print(tcap_result)
#' summary(tcap_result)
#'
#' # Compare with synthpop's disclosure function
#' disc_sp <- synthpop::disclosure(synth_obj, original_data,
#'                                  keys = c("sex", "age", "edu", "marital"),
#'                                  target = "income",
#'                                  print.flag = FALSE)
#'
#' # Compare with synthpop disclosure metrics
#' # Note: Core CAP algorithm matches exactly (verified with controlled data).
#' # Differences with real data are due to NA handling:
#' # - riskutility: removes NA records, computes mean over matched only
#' # - synthpop: uses all records with package-specific NA handling
#' cat("\nComparison with synthpop:\n")
#' cat("riskutility mean TCAP:", round(tcap_result$tcap_mean * 100, 2), "%\n")
#' cat("  (over", tcap_result$n_matched, "matched of", tcap_result$n_total, "complete cases)\n")
#' cat("synthpop CAPd:", round(disc_sp$allCAPs$CAPd, 2), "%\n")
#'
#' # Plot the TCAP distribution
#' plot(tcap_result, which = 1:2)
#' }
#' }
tcap <- function(X, ...) {
  UseMethod("tcap")
}

#' @rdname tcap
#' @export
tcap.synth_pair <- function(X, ...) {
  if (is.null(X$key_vars)) {
    stop("synth_pair must have 'key_vars' set for tcap()")
  }
  if (is.null(X$target_var)) {
    stop("synth_pair must have 'target_var' set for tcap()")
  }

  tcap.default(
    X = X$original,
    Y = X$synthetic,
    key_vars = X$key_vars,
    target_var = X$target_var,
    ...
  )
}

#' @rdname tcap
#' @export
tcap.default <- function(X, Y,
                         key_vars,
                         target_var,
                         method = c("exact", "gower"),
                         gower_threshold = 0.1,
                         cont_bins = 10,
                         na.rm = TRUE,
                         kind = c("certain", "matched", "conditional"),
                         ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

  method <- match.arg(method)
  kind <- match.arg(kind)

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

  # Get target values
  target_X <- X[[target_var]]
  target_Y <- Y[[target_var]]

  # Bin a continuous target, exactly as dcap() does, so that the members of the
  # CAP family are computed on the same discretisation of the target.
  if (is.numeric(target_X) && !is.factor(target_X)) {
    breaks <- quantile(c(target_X, target_Y),
                       probs = seq(0, 1, length.out = cont_bins + 1),
                       na.rm = TRUE)
    breaks <- unique(breaks)
    target_X <- cut(target_X, breaks = breaks, include.lowest = TRUE)
    target_Y <- cut(target_Y, breaks = breaks, include.lowest = TRUE)
  }

  # Compute baseline (modal target frequency in synthetic data)
  target_freq <- table(target_Y) / length(target_Y)
  baseline <- max(target_freq)

  # Initialize results
  n <- nrow(X)
  tcap_scores <- numeric(n)
  n_matches <- integer(n)
  is_certain <- logical(n)     # Track records at "certain" disclosure risk
  is_disclosive <- logical(n)  # Matched to a disclosive synthetic key class

  if (method == "exact") {
    # Vectorized exact matching using data.table for O(n + m) performance
    # Create key signature using fast vectorized paste
    make_key <- function(df, vars) {
      do.call(paste, c(df[, vars, drop = FALSE], sep = "|"))
    }

    keys_X <- make_key(X, key_vars)
    keys_Y <- make_key(Y, key_vars)

    # Build lookup tables using data.table
    dt_Y <- data.table::data.table(
      key_sig = keys_Y,
      target = target_Y,
      idx = seq_along(keys_Y)
    )

    # Pre-compute: for each key in Y, count total and count per target value
    key_counts <- dt_Y[, .(n_total = .N), by = key_sig]
    key_target_counts <- dt_Y[, .(n_correct = .N), by = .(key_sig, target)]

    # Pre-compute: for each key, is target unique in synthetic?
    key_unique_syn <- dt_Y[, .(unique_syn = length(unique(target)) == 1), by = key_sig]

    # Create lookup for original data
    dt_X <- data.table::data.table(
      orig_idx = seq_len(n),
      key_sig = keys_X,
      target = target_X
    )

    # Pre-compute: for each key, is target unique in original?
    key_unique_orig <- dt_X[, .(unique_orig = length(unique(target)) == 1), by = key_sig]

    # Join to get match counts
    dt_X <- data.table::merge.data.table(dt_X, key_counts, by = "key_sig", all.x = TRUE)
    dt_X <- data.table::merge.data.table(dt_X, key_target_counts,
                                          by = c("key_sig", "target"), all.x = TRUE)
    dt_X <- data.table::merge.data.table(dt_X, key_unique_orig, by = "key_sig", all.x = TRUE)
    dt_X <- data.table::merge.data.table(dt_X, key_unique_syn, by = "key_sig", all.x = TRUE)

    # Fill NAs with appropriate defaults
    dt_X[is.na(n_total), n_total := 0L]
    dt_X[is.na(n_correct), n_correct := 0L]
    dt_X[is.na(unique_orig), unique_orig := FALSE]
    dt_X[is.na(unique_syn), unique_syn := FALSE]

    # Restore original order
    data.table::setorder(dt_X, orig_idx)

    # Extract results
    n_matches <- dt_X$n_total
    tcap_scores <- ifelse(n_matches == 0, NA_real_, dt_X$n_correct / n_matches)

    # Certain disclosure: unique in both orig and syn, and CAP = 1
    is_certain <- dt_X$unique_orig & dt_X$unique_syn &
                  !is.na(tcap_scores) & tcap_scores == 1

    # Disclosive-matched: the record's synthetic key class has a single
    # target value (denominator of the conditional TCAP variant)
    is_disclosive <- dt_X$unique_syn & !is.na(tcap_scores)

  } else if (method == "gower") {
    # Gower distance-based matching (vectorized)
    # Note: "certain" disclosure is less meaningful with fuzzy matching,
    # so we use a simplified criterion (CAP = 1)
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

    # Compute TCAP scores
    tcap_scores <- ifelse(n_matches == 0, NA_real_, n_correct / n_matches)

    # For Gower matching, "certain" means CAP = 1 (all matches have correct target)
    # If CAP = 1, all matches share the same target value as the original
    is_certain <- !is.na(tcap_scores) & tcap_scores == 1

    # Disclosive-matched analogue: all matched synthetic records share one
    # target value, regardless of whether it is the correct one
    max_target_count <- rep(0, n)
    for (t in unique(target_Y)) {
      cnt_t <- rowSums(match_matrix[, target_Y == t, drop = FALSE])
      max_target_count <- pmax(max_target_count, cnt_t)
    }
    is_disclosive <- n_matches > 0 & max_target_count == n_matches
  }

  # Aggregate results
  n_matched <- sum(!is.na(tcap_scores))
  n_unmatched <- sum(is.na(tcap_scores))

  # Count records at "certain" disclosure risk (among matched records)
  n_certain <- sum(is_certain, na.rm = TRUE)
  tcap_certain <- if (n_matched > 0) 100 * n_certain / n_matched else NA

  # Shared numerator of the matched/conditional variants: disclosive in the
  # synthetic data and correctly attributed (CAP = 1)
  n_disclosive_correct <- sum(!is.na(tcap_scores) & tcap_scores == 1)
  n_disclosive <- sum(is_disclosive, na.rm = TRUE)

  # synthpop <= 1.9-2 TCAP: numerator over all matched records
  tcap_matched <- if (n_matched > 0) 100 * n_disclosive_correct / n_matched else NA
  # Little et al. (2025) / synthpop >= 1.9-3 TCAP: numerator over records
  # matched to a disclosive class; 0 when that set is empty, following the
  # reference implementation
  tcap_conditional <- if (n_disclosive > 0) {
    100 * n_disclosive_correct / n_disclosive
  } else {
    0
  }

  # Handle case when no matches exist (all NA scores)
  valid_scores <- tcap_scores[!is.na(tcap_scores)]
  tcap_max_val <- if (length(valid_scores) > 0) max(valid_scores) else NA_real_

  results <- list(
    tcap_scores = tcap_scores,
    tcap_mean = mean(tcap_scores, na.rm = TRUE),
    tcap_max = tcap_max_val,
    tcap_median = median(tcap_scores, na.rm = TRUE),
    tcap_certain = tcap_certain,
    tcap_matched = tcap_matched,
    tcap_conditional = tcap_conditional,
    n_certain = n_certain,
    n_disclosive_correct = n_disclosive_correct,
    n_disclosive = n_disclosive,
    is_certain = is_certain,
    is_disclosive = is_disclosive,
    kind = kind,
    n_matches = n_matches,
    n_matched = n_matched,
    n_unmatched = n_unmatched,
    n_total = n,
    baseline = baseline,
    key_vars = key_vars,
    target_var = target_var,
    method = method,
    gower_threshold = if (method == "gower") gower_threshold else NULL
  )

  class(results) <- "tcap"
  return(results)
}

#' Print method for tcap objects
#'
#' @param x an object of class "tcap"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.tcap <- function(x, ...) {
  cat("Correct Attribution Probability (CAP) Analysis\n")
  cat("===============================================\n")
  cat("Method:", x$method, "\n")
  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Target variable:", x$target_var, "\n\n")

  sel <- function(k) if (identical(x$kind, k)) " (selected)" else ""
  pct <- function(v) sprintf("(%.1f%%)", v)

  cat("CAP Metrics:\n")
  cat("  Mean CAP:", round(x$tcap_mean, 4), "\n")
  cat("  Max CAP:", round(x$tcap_max, 4), "\n")
  cat("  Median CAP:", round(x$tcap_median, 4), "\n")
  cat("  Baseline (random):", round(x$baseline, 4), "\n")
  cat("  Certain disclosure:", x$n_certain, "/", x$n_matched,
      paste0(pct(x$tcap_certain), sel("certain")), "\n")
  cat("  TCAP matched (synthpop <= 1.9-2):", x$n_disclosive_correct, "/",
      x$n_matched, paste0(pct(x$tcap_matched), sel("matched")), "\n")
  cat("  TCAP conditional (Little et al. 2025):", x$n_disclosive_correct, "/",
      x$n_disclosive, paste0(pct(x$tcap_conditional), sel("conditional")), "\n")
  cat("  Records matched:", x$n_matched, "/", x$n_total,
      sprintf("(%.1f%%)", 100 * x$n_matched / x$n_total), "\n")

  # Risk assessment
  risk_ratio <- x$tcap_mean / x$baseline
  cat("\nRisk Assessment:\n")
  cat("  CAP/Baseline ratio:", round(risk_ratio, 2))
  if (!is.finite(risk_ratio)) {
    cat(" (undefined - no matched records)\n")
  } else if (risk_ratio > 1.5) {
    cat(" (elevated risk)\n")
  } else if (risk_ratio > 1.0) {
    cat(" (slightly elevated)\n")
  } else {
    cat(" (low risk)\n")
  }

  if (isTRUE(x$tcap_certain > 10)) {
    cat("  Warning: >10% of records at certain disclosure risk\n")
  }

  invisible(x)
}

#' Summary method for tcap objects
#'
#' @param object an object of class "tcap"
#' @param ... additional arguments (ignored)
#' @return A list of summary statistics for the corresponding object.
#' @export
summary.tcap <- function(object, ...) {
  scores <- object$tcap_scores[!is.na(object$tcap_scores)]

  # Categorize risk levels
  high_risk <- sum(scores >= 0.8)
  medium_risk <- sum(scores >= 0.5 & scores < 0.8)
  low_risk <- sum(scores < 0.5)

  summ <- list(
    tcap_mean = object$tcap_mean,
    tcap_max = object$tcap_max,
    tcap_median = object$tcap_median,
    tcap_sd = sd(scores),
    tcap_quantiles = quantile(scores, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1)),
    tcap_certain = object$tcap_certain,
    n_certain = object$n_certain,
    tcap_matched = object$tcap_matched,
    tcap_conditional = object$tcap_conditional,
    n_disclosive_correct = object$n_disclosive_correct,
    n_disclosive = object$n_disclosive,
    kind = object$kind,
    baseline = object$baseline,
    risk_ratio = object$tcap_mean / object$baseline,
    n_high_risk = high_risk,
    n_medium_risk = medium_risk,
    n_low_risk = low_risk,
    pct_high_risk = 100 * high_risk / length(scores),
    n_matched = object$n_matched,
    n_unmatched = object$n_unmatched,
    n_total = object$n_total,
    mean_matches = mean(object$n_matches[object$n_matches > 0]),
    key_vars = object$key_vars,
    target_var = object$target_var,
    method = object$method
  )

  class(summ) <- "summary.tcap"
  return(summ)
}

#' Print method for summary.tcap objects
#'
#' @param x an object of class "summary.tcap"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.tcap <- function(x, ...) {
  cat("Summary: Correct Attribution Probability (CAP) Analysis\n")
  cat("=======================================================\n")
  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Target variable:", x$target_var, "\n")
  cat("Method:", x$method, "\n\n")

  cat("CAP Statistics:\n")
  cat("  Mean (CAPd):", round(x$tcap_mean, 4), "\n")
  cat("  Max:", round(x$tcap_max, 4), "\n")
  cat("  Median:", round(x$tcap_median, 4), "\n")
  cat("  SD:", round(x$tcap_sd, 4), "\n")
  cat("  Baseline:", round(x$baseline, 4), "\n")
  cat("  Risk ratio:", round(x$risk_ratio, 2), "\n\n")

  mark <- function(k) if (identical(x$kind, k)) "* " else "  "
  pctv <- function(v) sprintf("%.1f%%", v)

  cat("Certain Disclosure and TCAP Variants (* = selected kind):\n")
  cat("  ", mark("certain"), "certain:     ", pctv(x$tcap_certain),
      " (", x$n_certain, " / ", x$n_matched,
      " matched; unique key->target in original & synthetic)\n", sep = "")
  cat("  ", mark("matched"), "matched:     ", pctv(x$tcap_matched),
      " (", x$n_disclosive_correct, " / ", x$n_matched,
      " matched; synthpop <= 1.9-2 TCAP)\n", sep = "")
  cat("  ", mark("conditional"), "conditional: ", pctv(x$tcap_conditional),
      " (", x$n_disclosive_correct, " / ", x$n_disclosive,
      " disclosive-matched; Little et al. 2025; synthpop >= 1.9-3)\n\n", sep = "")

  cat("CAP Distribution (quantiles):\n")
  print(round(x$tcap_quantiles, 4))
  cat("\n")

  cat("Risk Categories:\n")
  cat("  High risk (CAP >= 0.8):", x$n_high_risk,
      sprintf("(%.1f%%)", x$pct_high_risk), "\n")
  cat("  Medium risk (0.5 <= CAP < 0.8):", x$n_medium_risk, "\n")
  cat("  Low risk (CAP < 0.5):", x$n_low_risk, "\n\n")

  cat("Matching Statistics:\n")
  cat("  Records matched:", x$n_matched, "/", x$n_total, "\n")
  cat("  Records unmatched:", x$n_unmatched, "\n")
  cat("  Avg matches per record:", round(x$mean_matches, 1), "\n")

  invisible(x)
}

#' Plot method for tcap objects
#'
#' @param x an object of class "tcap"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot(s) to show: 1 = TCAP histogram, 2 = risk categories,
#'   3 = matches distribution
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.tcap <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 3)
  show[which] <- TRUE

  scores <- x$tcap_scores[!is.na(x$tcap_scores)]

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
    # Histogram of CAP scores
    hist(scores, breaks = 20,
         main = "Distribution of CAP Scores",
         xlab = "CAP Score", col = "steelblue", border = "white", ...)
    abline(v = x$tcap_mean, col = "red", lwd = 2, lty = 2)
    abline(v = x$tcap_max, col = "purple", lwd = 2, lty = 4)
    abline(v = x$baseline, col = "darkgreen", lwd = 2, lty = 3)
    legend("topright", legend = c("Mean CAP", "Max CAP", "Baseline"),
           col = c("red", "purple", "darkgreen"), lty = c(2, 4, 3), lwd = 2, cex = 0.8)
  }

  if (show[2]) {
    # Risk categories bar plot
    high_risk <- sum(scores >= 0.8)
    medium_risk <- sum(scores >= 0.5 & scores < 0.8)
    low_risk <- sum(scores < 0.5)

    counts <- c(low_risk, medium_risk, high_risk)
    names(counts) <- c("Low\n(<0.5)", "Medium\n(0.5-0.8)", "High\n(>=0.8)")
    barplot(counts, main = "Risk Categories",
            ylab = "Number of Records",
            col = c("forestgreen", "orange", "firebrick"), ...)
  }

  if (show[3]) {
    # Matches distribution
    matches <- x$n_matches[x$n_matches > 0]
    hist(matches, breaks = 30,
         main = "Matches per Record",
         xlab = "Number of Synthetic Matches",
         col = "coral", border = "white", ...)
    abline(v = mean(matches), col = "red", lwd = 2, lty = 2)
  }
}

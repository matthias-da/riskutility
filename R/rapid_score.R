# RAPID scoring core: tolerance band, baseline, gain and the rule classes
#
# A *rule* says how one record is scored and when it is flagged. Each rule is a
# subclass of `rapid_rule` with one `.rule_score()` method. Adding a rule means
# adding one branch in `.rapid_rule()` and one method here; nothing outside this
# file may branch on the name of a rule.

# ---- Error metrics and tolerance band ----------------------------------------

.rapid_metrics <- c("symmetric", "stabilised_relative", "absolute")

#' Error of a value against the true value
#' @noRd
.rapid_error <- function(y, v, metric, delta) {
  switch(metric,
    symmetric = 2 * abs(y - v) / (abs(y) + abs(v) + 2 * delta),
    stabilised_relative = abs(y - v) / (abs(y) + delta),
    absolute = abs(y - v),
    .rapid_abort(c("Unknown error metric '", metric, "'. Use one of: ",
                   paste(.rapid_metrics, collapse = ", "), "."), "input")
  )
}

#' Tolerance band around each true value
#'
#' The open interval of values `v` with `e(y, v) < epsilon`, solved in closed
#' form for each error metric.
#'
#' @param y Numeric vector of true values.
#' @param metric One of `"symmetric"`, `"stabilised_relative"`, `"absolute"`.
#' @param epsilon Tolerance on the internal scale (a proportion for the two
#'   relative metrics).
#' @param delta Smoothing constant of the relative metrics.
#' @return Data frame with the columns `lo` and `hi`.
#' @noRd
.rapid_band_bounds <- function(y, metric, epsilon, delta) {
  switch(metric,
    absolute = data.frame(lo = y - epsilon, hi = y + epsilon),
    stabilised_relative = {
      w <- epsilon * (abs(y) + delta)
      data.frame(lo = y - w, hi = y + w)
    },
    symmetric = {
      # Solved for a >= 0; negative values follow by symmetry (v -> -v).
      a <- abs(y)
      if (epsilon >= 2) {
        # The symmetric error never reaches 2, so every value is inside.
        lo <- rep(-Inf, length(a))
        hi <- rep(Inf, length(a))
      } else {
        hi <- (a * (2 + epsilon) + 2 * epsilon * delta) / (2 - epsilon)
        lo <- (a * (2 - epsilon) - 2 * epsilon * delta) / (2 + epsilon)
        # Band crosses zero: below zero the error is 2 (a - v) / (a - v + 2 delta)
        crosses <- lo < 0
        lo[crosses] <- a[crosses] - 2 * epsilon * delta / (2 - epsilon)
      }
      neg <- y < 0
      data.frame(lo = ifelse(neg, -hi, lo), hi = ifelse(neg, -lo, hi))
    },
    .rapid_abort(c("Unknown error metric '", metric, "'. Use one of: ",
                   paste(.rapid_metrics, collapse = ", "), "."), "input")
  )
}

#' Marginal baseline of a numeric target
#'
#' Share of the reference values strictly inside each band, counted on the
#' sorted values.
#' @noRd
.rapid_baseline_num <- function(lo, hi, ref_y) {
  s <- sort(ref_y)
  n_below_hi <- findInterval(hi, s, left.open = TRUE)  # number of s <  hi
  n_upto_lo <- findInterval(lo, s)                      # number of s <= lo
  (n_below_hi - n_upto_lo) / length(s)
}

#' Normalised gain over the baseline
#'
#' The gain is `(g - b) / (1 - b)`, and 0 where `b = 1` (no gain is possible).
#' @noRd
.rapid_gain <- function(g, b) {
  r <- numeric(length(g))
  open <- b < 1
  r[open] <- (g[open] - b[open]) / (1 - b[open])
  r
}

# ---- Predictive distributions --------------------------------------------------
# Task 1 defines the grid type; the exact types follow in rapid_models.R.

#' Predictive distribution given as a grid of values per record
#' @param D Numeric matrix, one row per record.
#' @noRd
.pdist_grid <- function(D) {
  if (!is.matrix(D) || !is.numeric(D)) {
    .rapid_abort("A grid distribution needs a numeric matrix with one row per record.",
                 "input")
  }
  structure(list(D = D, n = nrow(D)), class = c("rapid_pdist_grid", "rapid_pdist"))
}

#' Probability mass of a predictive distribution inside open intervals
#' @param pd A `rapid_pdist` object.
#' @param lo,hi Interval bounds, one per record.
#' @noRd
.pdist_mass <- function(pd, lo, hi) UseMethod(".pdist_mass")

#' @exportS3Method .pdist_mass rapid_pdist_grid
#' @noRd
.pdist_mass.rapid_pdist_grid <- function(pd, lo, hi) {
  rowMeans(pd$D > lo & pd$D < hi)
}

# ---- Rules -------------------------------------------------------------------

#' Scoring rule
#'
#' @param method `"band"`, `"point"`, `"RCS_marginal"`, `"RCS_conditional"` or
#'   `"NCE"`.
#' @param threshold Flag threshold. For `"point"` it is the tolerance.
#' @param metric,epsilon,delta Error metric, tolerance and smoothing constant
#'   (numeric rules only).
#' @return A list of class `c("rapid_rule_<method>", "rapid_rule")`. The
#'   constructor sets `score_var` (the record column that is compared with the
#'   threshold), `direction` (`"above"` or `"below"`) and `gain` (is the score a
#'   normalised gain over a marginal baseline).
#' @noRd
.rapid_rule <- function(method, threshold, metric = NULL, epsilon = NULL, delta = NULL) {
  spec <- switch(as.character(method)[1],
    band = list(score_var = "r", direction = "above", gain = TRUE, numeric = TRUE),
    point = list(score_var = "error", direction = "below", gain = FALSE, numeric = TRUE),
    RCS_marginal = list(score_var = "r", direction = "above", gain = TRUE, numeric = FALSE),
    RCS_conditional = list(score_var = "ratio", direction = "above", gain = FALSE,
                           numeric = FALSE),
    NCE = list(score_var = "nce", direction = "above", gain = FALSE, numeric = FALSE),
    .rapid_abort(c("Unknown evaluation method '", method, "'. Use \"band\" or \"point\" ",
                   "for a numeric target and \"RCS_marginal\", \"RCS_conditional\" or ",
                   "\"NCE\" for a categorical one."), "input")
  )
  if (!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold)) {
    .rapid_abort("The threshold of a rule must be a single number.", "input")
  }
  rule <- list(method = method, threshold = threshold,
               metric = NULL, epsilon = NULL, delta = NULL,
               score_var = spec$score_var, direction = spec$direction, gain = spec$gain)
  if (spec$numeric) {
    if (is.null(metric) || !metric %in% .rapid_metrics) {
      .rapid_abort(c("A numeric rule needs an error metric: ",
                     paste(.rapid_metrics, collapse = ", "), "."), "input")
    }
    rule["metric"] <- list(metric)
    rule["delta"] <- list(if (is.null(delta)) 0.01 else delta)
    # For the point rule the threshold is the tolerance itself.
    rule["epsilon"] <- list(if (spec$direction == "below") threshold else epsilon)
    if (!is.numeric(rule$epsilon) || length(rule$epsilon) != 1 || is.na(rule$epsilon) ||
        rule$epsilon <= 0) {
      .rapid_abort("The tolerance of a numeric rule must be a single positive number.",
                   "tolerance")
    }
  }
  structure(rule, class = c(paste0("rapid_rule_", method), "rapid_rule"))
}

#' Flag records by comparing their score with the threshold of the rule
#' @noRd
.rule_flag <- function(rule, score) {
  flag <- if (rule$direction == "below") score < rule$threshold else score > rule$threshold
  !is.na(flag) & flag
}

#' Replace the threshold of a rule
#'
#' Under the point rule the threshold is the tolerance, so both change.
#' @noRd
.rule_set_threshold <- function(rule, value) {
  rule$threshold <- value
  if (rule$direction == "below") rule$epsilon <- value
  rule
}

.rule_tolerance_label <- function(rule) {
  tol <- if (rule$metric == "absolute") {
    format(rule$epsilon)
  } else {
    paste(format(100 * rule$epsilon), "%")
  }
  paste0(gsub("_", " ", rule$metric, fixed = TRUE), " error < ", tol)
}

#' @export
format.rapid_rule <- function(x, ...) {
  # Fallback for a rule subclass without its own format method
  paste0(x$method, ": ", x$score_var, if (x$direction == "below") " < " else " > ",
         format(x$threshold))
}

#' @export
format.rapid_rule_band <- function(x, ...) {
  paste0("band: ", .rule_tolerance_label(x), ", gain > ", format(x$threshold))
}

#' @export
format.rapid_rule_point <- function(x, ...) {
  paste0("point: ", .rule_tolerance_label(x))
}

#' @export
format.rapid_rule_RCS_marginal <- function(x, ...) {
  paste0("RCS_marginal: gain over the marginal baseline > ", format(x$threshold))
}

#' @export
format.rapid_rule_RCS_conditional <- function(x, ...) {
  paste0("RCS_conditional: ratio to the class-conditional baseline > ",
         format(x$threshold))
}

#' @export
format.rapid_rule_NCE <- function(x, ...) {
  paste0("NCE: normalised cross-entropy score > ", format(x$threshold))
}

#' Score records under a rule
#'
#' @param rule A `rapid_rule`.
#' @param pred The attacker's prediction: a probability matrix with the class
#'   names as column names (categorical target), or
#'   `list(point = <numeric>, dist = <rapid_pdist>)` (numeric target).
#' @param y True values of the scored records.
#' @param ref_y Reference values for the marginal baseline (usually `y`).
#' @return Data frame with the record columns of the rule, one row per record.
#' @noRd
.rule_score <- function(rule, pred, y, ref_y) UseMethod(".rule_score")

.rule_check_numeric <- function(pred, y, need_dist) {
  if (!is.list(pred) || !is.numeric(pred$point) || length(pred$point) != length(y)) {
    .rapid_abort("Numeric predictions must hold one point prediction per record.", "input")
  }
  if (need_dist && (!inherits(pred$dist, "rapid_pdist") || pred$dist$n != length(y))) {
    .rapid_abort(c("The band method needs a predictive distribution with one entry per ",
                   "record."), "input")
  }
  invisible(TRUE)
}

#' @exportS3Method .rule_score rapid_rule_band
#' @noRd
.rule_score.rapid_rule_band <- function(rule, pred, y, ref_y) {
  .rule_check_numeric(pred, y, need_dist = TRUE)
  band <- .rapid_band_bounds(y, rule$metric, rule$epsilon, rule$delta)
  g <- .pdist_mass(pred$dist, band$lo, band$hi)
  b <- .rapid_baseline_num(band$lo, band$hi, ref_y)
  data.frame(predicted_value = as.numeric(pred$point), band_lo = band$lo,
             band_hi = band$hi, g = g, b = b, r = .rapid_gain(g, b))
}

#' @exportS3Method .rule_score rapid_rule_point
#' @noRd
.rule_score.rapid_rule_point <- function(rule, pred, y, ref_y) {
  .rule_check_numeric(pred, y, need_dist = FALSE)
  band <- .rapid_band_bounds(y, rule$metric, rule$epsilon, rule$delta)
  data.frame(predicted_value = as.numeric(pred$point),
             error = .rapid_error(y, as.numeric(pred$point), rule$metric, rule$delta),
             b = .rapid_baseline_num(band$lo, band$hi, ref_y))
}

# Probability of the true class and the predicted class, shared by the
# categorical rules.
.rule_class_probs <- function(pred, y) {
  P <- if (is.data.frame(pred)) as.matrix(pred) else pred
  if (!is.matrix(P) || !is.numeric(P) || is.null(colnames(P)) || nrow(P) != length(y)) {
    .rapid_abort(c("Class predictions must be a probability matrix with one row per ",
                   "record and the class names as column names."), "input")
  }
  col <- match(as.character(y), colnames(P))
  if (anyNA(col)) {
    .rapid_abort(c("No predicted probability for the class(es): ",
                   paste(unique(as.character(y)[is.na(col)]), collapse = ", "), "."),
                 "input")
  }
  list(P = P, g = unname(P[cbind(seq_len(nrow(P)), col)]),
       predicted_class = colnames(P)[max.col(P, ties.method = "first")])
}

#' @exportS3Method .rule_score rapid_rule_RCS_marginal
#' @noRd
.rule_score.rapid_rule_RCS_marginal <- function(rule, pred, y, ref_y) {
  p <- .rule_class_probs(pred, y)
  share <- table(as.character(ref_y)) / length(ref_y)
  b <- as.numeric(share[as.character(y)])
  b[is.na(b)] <- 0                      # class absent from the reference values
  data.frame(predicted_class = p$predicted_class, g = p$g, b = b,
             r = .rapid_gain(p$g, b), stringsAsFactors = FALSE)
}

#' @exportS3Method .rule_score rapid_rule_RCS_conditional
#' @noRd
.rule_score.rapid_rule_RCS_conditional <- function(rule, pred, y, ref_y) {
  p <- .rule_class_probs(pred, y)
  # Baseline: mean probability of the true class among the records of that class
  b <- as.numeric(stats::ave(p$g, as.character(y)))
  data.frame(predicted_class = p$predicted_class, g = p$g, b = b, ratio = p$g / b,
             stringsAsFactors = FALSE)
}

#' @exportS3Method .rule_score rapid_rule_NCE
#' @noRd
.rule_score.rapid_rule_NCE <- function(rule, pred, y, ref_y) {
  p <- .rule_class_probs(pred, y)
  nce <- 1 + log(pmax(p$g, 1e-10)) / log(ncol(p$P))
  data.frame(predicted_class = p$predicted_class, g = p$g, nce = nce,
             stringsAsFactors = FALSE)
}

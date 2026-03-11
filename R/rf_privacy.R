# R/rf_privacy.R

#' RF Proximity Privacy Assessment
#'
#' Detects memorization in synthetic data using Random Forest proximity.
#' Trains a supervised RF to discriminate original from synthetic records,
#' then compares proximity of synthetic records to training vs. holdout
#' subsets of the original data.
#'
#' @param X For \code{rf_privacy.default}: a data.frame of original data.
#'   For \code{rf_privacy.synth_pair}: a \code{synth_pair} object.
#' @param Y data.frame of synthetic data (for default method)
#' @param holdout data.frame or NULL. If NULL, split from X.
#' @param holdout_fraction numeric in (0, 1), fraction of X to use as holdout
#' @param vars character vector of variable names (NULL = all common)
#' @param na.rm logical. If FALSE (default), ranger handles NAs natively
#'   via surrogate splits. If TRUE, remove records with any NA before training.
#' @param seed integer or NULL. Used for holdout splitting (seed) and forest
#'   training (seed + 1) to avoid correlation.
#' @param progress logical, show progress bar during proximity computation
#' @param null_test logical, run permutation null test (default TRUE)
#' @param n_null integer, number of permutations for null test
#' @param n_trees integer (>= 10), number of trees
#' @param mtry integer or NULL, number of variables to consider at each split
#' @param ... additional arguments passed to ranger via modifyList
#'
#' @return An S3 object of class \code{"rf_privacy"} with fields:
#'   \describe{
#'     \item{max_prox_share}{fraction of synthetic records with higher max
#'       proximity to training than holdout (mid-rank tie correction)}
#'     \item{max_prox_ratio}{ratio of mean max-proximities (train/holdout)}
#'     \item{max_prox_train}{per-record max proximity to nearest training record}
#'     \item{max_prox_holdout}{per-record max proximity to nearest holdout record}
#'     \item{prox_share}{fraction with higher mean proximity to training}
#'     \item{prox_ratio}{ratio of mean mean-proximities}
#'     \item{prox_train_mean}{per-record mean proximity to training}
#'     \item{prox_holdout_mean}{per-record mean proximity to holdout}
#'     \item{privacy_pass}{logical, TRUE if no memorization detected}
#'     \item{wilcox_test}{Wilcoxon signed-rank test object (heuristic)}
#'     \item{null_distribution}{list with null stats and p-values (if null_test)}
#'     \item{oob_error}{OOB classification error from the forest}
#'     \item{var_importance}{named numeric vector of variable importances}
#'     \item{n_synthetic, n_train, n_holdout, vars}{dataset metadata}
#'   }
#'
#' @details
#' RF proximity measures how often two records land in the same terminal node
#' across all trees. A proximity of 1 means they always co-terminate; 0 means
#' never. This function trains a supervised RF to discriminate ALL original
#' records from synthetic records, then checks whether synthetic records are
#' more similar (proximate) to training records than to holdout records.
#'
#' Crucially, the forest is trained on the full original dataset (train +
#' holdout combined) so that both subsets have identical in-sample status,
#' eliminating the bias that would arise if holdout records were out-of-sample.
#'
#' If no memorization occurred, synthetic records should have roughly equal
#' proximity to training and holdout (both are real data, both in the forest).
#' If synthetic records consistently land in terminal nodes dominated by
#' training records, this signals memorization.
#'
#' @section When to use this method:
#' Use \code{rf_privacy()} when you have mixed data types with complex
#' interactions (20+ variables), or when you want a data-adaptive alternative
#' to \code{\link{dcr}()}. For simple QIs with interpretable risk, \code{dcr()}
#' with Gower distance is more transparent. For speed on small data (n < 5,000),
#' \code{dcr()} is faster.
#'
#' @section Interpretation:
#' \code{max_prox_share} is the primary metric. Values near 0.5 indicate no
#' memorization; values above 0.5 suggest synthetic records are systematically
#' closer to training than holdout. The \code{max_prox_ratio} captures the
#' magnitude: values near 1 = no memorization, > 1 = memorization signal.
#'
#' Mean-based metrics (\code{prox_share}, \code{prox_ratio}) detect aggregate
#' distributional leakage rather than individual memorized records.
#'
#' @section Comparison with DCR:
#' \code{rf_privacy()} is the RF-proximity analog of \code{dcr()}, using the
#' same holdout design but replacing Gower/Euclidean distance with terminal-node
#' co-occurrence. RF proximity is data-adaptive and handles mixed types
#' natively, but adds forest training overhead.
#'
#' @section Limitations:
#' The supervised forest conflates distributional similarity with memorization.
#' When OOB error is near 0.5 (high-quality synthetic data), memorization
#' detection is most reliable. When OOB error is low (poor utility), asymmetry
#' may reflect distributional differences, not individual copying.
#'
#' The Wilcoxon p-value is anti-conservative because all proximity values share
#' the same forest. Use the permutation null test (\code{null_test = TRUE}) for
#' principled inference.
#'
#' @section Computational considerations:
#' Expected runtimes (n_trees = 500, modern laptop):
#' n = 1,000: ~5 seconds; n = 5,000: ~30 seconds;
#' n = 10,000: 2-5 minutes; n = 50,000: 30+ minutes.
#'
#' @family distance-risk
#' @seealso \code{\link{dcr}}, \code{\link{nndr}}, \code{\link{ims}},
#'   \code{\link{propscore}}, \code{\link{recordLinkage}}
#' @references
#' Breiman, L. (2001). Random Forests. \emph{Machine Learning}, 45(1), 5-32.
#'
#' Lin, Y. & Jeon, Y. (2006). Random forests and adaptive nearest neighbors.
#' \emph{JASA}, 101(474), 578-590.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X_train <- data.frame(a = rnorm(100), b = rnorm(100), c = rnorm(100))
#' X_holdout <- data.frame(a = rnorm(100), b = rnorm(100), c = rnorm(100))
#'
#' # Memorized: synthetic copies training records
#' Y_mem <- X_train[sample(100, 100, replace = TRUE), ]
#' Y_mem <- Y_mem + rnorm(300, 0, 0.01)
#' res_mem <- rf_privacy(X_train, Y_mem, holdout = X_holdout,
#'                       seed = 1, n_trees = 200, null_test = FALSE)
#' print(res_mem)
#'
#' # Random synthetic data (no memorization)
#' Y_rand <- data.frame(a = rnorm(100), b = rnorm(100), c = rnorm(100))
#' res_rand <- rf_privacy(X_train, Y_rand, holdout = X_holdout,
#'                        seed = 1, n_trees = 200, null_test = FALSE)
#' print(res_rand)
#' }
#' @export
rf_privacy <- function(X, ...) UseMethod("rf_privacy")

#' @rdname rf_privacy
#' @export
rf_privacy.synth_pair <- function(X, ...) {
  rf_privacy.default(X = X$original, Y = X$synthetic,
                     vars = X$vars, ...)
}

#' @rdname rf_privacy
#' @export
rf_privacy.default <- function(X, Y,
                               holdout = NULL, holdout_fraction = 0.5,
                               vars = NULL, na.rm = FALSE, seed = NULL,
                               progress = FALSE,
                               null_test = TRUE, n_null = 100L,
                               n_trees = 500L, mtry = NULL, ...) {

  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' required for rf_privacy(). ",
         "Install with install.packages('ranger')", call. = FALSE)
  }

  # Prepare holdout split (splits X into train + holdout)
  prep <- .distance_risk_prepare(X, Y,
                                  holdout = holdout,
                                  holdout_fraction = holdout_fraction,
                                  vars = vars, na.rm = na.rm,
                                  seed = seed, min_holdout = 1L)
  train     <- prep$train
  synthetic <- prep$synthetic
  holdout_d <- prep$holdout
  vars_use  <- prep$vars

  n_syn   <- nrow(synthetic)
  n_train <- nrow(train)
  n_ho    <- nrow(holdout_d)

  # Recombine train + holdout for forest training to avoid in-sample bias.
  # Both subsets are included in the forest so proximity is symmetric.
  all_original <- rbind(train, holdout_d)
  n_orig <- n_train + n_ho

  if (n_orig + n_syn > 10000) {
    message("rf_privacy: n = ", n_orig + n_syn,
            ". Expected runtime: 2-5 minutes with n_trees = ", n_trees, ".")
  }

  # Train RF on ALL original vs synthetic
  rf_seed <- if (!is.null(seed)) seed + 1L else NULL
  rf_res <- .rf_proximity(all_original, synthetic, vars = vars_use,
                          n_trees = n_trees, mtry = mtry,
                          importance = TRUE, seed = rf_seed, ...)

  # Indices in the combined terminal_nodes matrix:
  # rows 1:n_train              = training records
  # rows (n_train+1):n_orig     = holdout records
  # rows (n_orig+1):(n_orig+n_syn) = synthetic records
  idx_train <- seq_len(n_train)
  idx_ho    <- n_train + seq_len(n_ho)
  idx_syn   <- n_orig + seq_len(n_syn)

  # Compute proximity matrices.

  # When null_test = TRUE we need the full syn-vs-all-original matrix anyway

  # for the permutation loop, so compute it once and slice. When null_test =
  # FALSE, compute only the two subsets to save memory/time.
  if (null_test) {
    idx_all_orig <- seq_len(n_orig)
    prox_syn_all <- .proximity_from_nodes(rf_res$terminal_nodes,
                                          idx_all_orig, idx_syn,
                                          progress = progress)
    # prox_syn_all is n_syn x n_orig
    prox_syn_train <- prox_syn_all[, seq_len(n_train), drop = FALSE]
    prox_syn_ho    <- prox_syn_all[, (n_train + 1):n_orig, drop = FALSE]
  } else {
    prox_syn_train <- .proximity_from_nodes(rf_res$terminal_nodes,
                                            idx_train, idx_syn,
                                            progress = progress)
    prox_syn_ho    <- .proximity_from_nodes(rf_res$terminal_nodes,
                                            idx_ho, idx_syn,
                                            progress = FALSE)
  }
  # prox_syn_train is n_syn x n_train
  # prox_syn_ho    is n_syn x n_ho

  # Per-record max and mean proximity
  max_prox_train   <- apply(prox_syn_train, 1, max)
  max_prox_holdout <- apply(prox_syn_ho, 1, max)
  mean_prox_train  <- rowMeans(prox_syn_train)
  mean_prox_ho     <- rowMeans(prox_syn_ho)

  # Max-based metrics (primary)
  # Mid-rank tie correction: mean(a > b) + 0.5 * mean(a == b)
  max_prox_share <- mean(max_prox_train > max_prox_holdout) +
    0.5 * mean(max_prox_train == max_prox_holdout)

  # Zero-guard for ratio
  denom_max <- mean(max_prox_holdout)
  if (denom_max < 1 / n_trees) {
    max_prox_ratio <- NA_real_
    warning("max_prox_holdout mean < 1/n_trees (degenerate). ",
            "max_prox_ratio set to NA.", call. = FALSE)
  } else {
    max_prox_ratio <- mean(max_prox_train) / denom_max
  }

  # Mean-based metrics (supplementary)
  prox_share <- mean(mean_prox_train > mean_prox_ho) +
    0.5 * mean(mean_prox_train == mean_prox_ho)

  denom_mean <- mean(mean_prox_ho)
  if (denom_mean < 1 / n_trees) {
    prox_ratio <- NA_real_
    warning("prox_holdout_mean mean < 1/n_trees (degenerate). ",
            "prox_ratio set to NA.", call. = FALSE)
  } else {
    prox_ratio <- mean(mean_prox_train) / denom_mean
  }

  # Wilcoxon signed-rank test (heuristic -- anti-conservative)
  wilcox_res <- suppressWarnings(
    stats::wilcox.test(max_prox_train, max_prox_holdout,
                       paired = TRUE, alternative = "greater")
  )

  # Null test: permute train/holdout labels
  null_dist <- NULL
  if (null_test) {
    # prox_syn_all (n_syn x n_orig) was already computed above

    null_shares <- numeric(n_null)
    null_ratios <- numeric(n_null)
    null_max_shares <- numeric(n_null)
    null_max_ratios <- numeric(n_null)

    for (p in seq_len(n_null)) {
      perm <- sample(n_orig)
      perm_train <- perm[seq_len(n_train)]
      perm_ho    <- perm[(n_train + 1):n_orig]

      # Mean-based
      perm_train_mean <- rowMeans(prox_syn_all[, perm_train, drop = FALSE])
      perm_ho_mean    <- rowMeans(prox_syn_all[, perm_ho, drop = FALSE])
      null_shares[p] <- mean(perm_train_mean > perm_ho_mean) +
        0.5 * mean(perm_train_mean == perm_ho_mean)

      d_ho <- mean(perm_ho_mean)
      null_ratios[p] <- if (d_ho < 1 / n_trees) NA_real_ else
        mean(perm_train_mean) / d_ho

      # Max-based
      perm_train_max <- apply(prox_syn_all[, perm_train, drop = FALSE], 1, max)
      perm_ho_max    <- apply(prox_syn_all[, perm_ho, drop = FALSE], 1, max)
      null_max_shares[p] <- mean(perm_train_max > perm_ho_max) +
        0.5 * mean(perm_train_max == perm_ho_max)

      d_ho_max <- mean(perm_ho_max)
      null_max_ratios[p] <- if (d_ho_max < 1 / n_trees) NA_real_ else
        mean(perm_train_max) / d_ho_max
    }

    # Permutation p-values: (sum + 1) / (n_null + 1) [Phipson & Smyth]
    pval_share <- (sum(null_max_shares >= max_prox_share, na.rm = TRUE) + 1) /
      (n_null + 1)
    if (is.na(max_prox_ratio)) {
      pval_ratio <- NA_real_
    } else {
      pval_ratio <- (sum(null_max_ratios >= max_prox_ratio,
                         na.rm = TRUE) + 1) / (n_null + 1)
    }

    null_dist <- list(
      null_max_shares = null_max_shares,
      null_max_ratios = null_max_ratios,
      null_shares = null_shares,
      null_ratios = null_ratios,
      null_pvalue = pval_share,
      null_ratio_pvalue = pval_ratio
    )
  }

  # Privacy pass
  if (!is.null(null_dist)) {
    # Null-test-derived: observed within 95th percentile
    q95_share <- stats::quantile(null_dist$null_max_shares, 0.95,
                                 na.rm = TRUE)
    q95_ratio <- stats::quantile(null_dist$null_max_ratios, 0.95,
                                 na.rm = TRUE)
    privacy_pass <- (max_prox_share <= q95_share) &&
      (is.na(max_prox_ratio) || max_prox_ratio <= q95_ratio)
  } else {
    # Heuristic fallback
    privacy_pass <- max_prox_share <= 0.55 &&
      wilcox_res$p.value > 0.05
  }

  result <- list(
    max_prox_share   = max_prox_share,
    max_prox_ratio   = max_prox_ratio,
    max_prox_train   = max_prox_train,
    max_prox_holdout = max_prox_holdout,
    prox_share       = prox_share,
    prox_ratio       = prox_ratio,
    prox_train_mean  = mean_prox_train,
    prox_holdout_mean = mean_prox_ho,
    privacy_pass     = privacy_pass,
    wilcox_test      = wilcox_res,
    null_distribution = null_dist,
    oob_error        = rf_res$oob_error,
    var_importance   = rf_res$importance,
    n_synthetic      = n_syn,
    n_train          = n_train,
    n_holdout        = n_ho,
    vars             = vars_use
  )
  class(result) <- "rf_privacy"
  result
}

#' @export
print.rf_privacy <- function(x, ...) {
  pass_label <- if (x$privacy_pass) "PASS" else "FAIL"
  cat("RF Privacy Assessment (rf_privacy)\n")
  cat("  Max proximity share: ", sprintf("%.2f", x$max_prox_share),
      if (x$max_prox_share <= 0.55) " (training not preferred)"
      else " (training preferred -- memorization signal)", "\n")
  cat("  Max proximity ratio: ",
      if (is.na(x$max_prox_ratio)) "NA (degenerate holdout)"
      else sprintf("%.2f", x$max_prox_ratio), "\n")
  cat("  OOB error:           ", sprintf("%.2f", x$oob_error), "\n")
  if (!is.null(x$null_distribution)) {
    cat("  Null test:           ", pass_label,
        sprintf(" (p = %.2f)", x$null_distribution$null_pvalue), "\n")
  }
  cat("  Privacy:             ", pass_label, "\n")
  invisible(x)
}

#' @export
summary.rf_privacy <- function(object, ...) {
  s <- list(
    max_prox_share   = object$max_prox_share,
    max_prox_ratio   = object$max_prox_ratio,
    prox_share       = object$prox_share,
    prox_ratio       = object$prox_ratio,
    oob_error        = object$oob_error,
    privacy_pass     = object$privacy_pass,
    wilcox_p         = object$wilcox_test$p.value,
    null_distribution = object$null_distribution,
    var_importance   = object$var_importance,
    n_synthetic      = object$n_synthetic,
    n_train          = object$n_train,
    n_holdout        = object$n_holdout,
    # Per-record outliers: top 5 by max_prox_train
    top_outliers     = utils::head(
      order(object$max_prox_train, decreasing = TRUE), 5
    )
  )
  class(s) <- "summary.rf_privacy"
  s
}

#' @export
print.summary.rf_privacy <- function(x, ...) {
  cat("RF Privacy Assessment -- Summary\n")
  cat(strrep("-", 50), "\n")
  cat("Dataset: ", x$n_synthetic, " synthetic, ",
      x$n_train, " train, ", x$n_holdout, " holdout\n\n")

  cat("Max-based metrics (primary):\n")
  cat("  max_prox_share: ", sprintf("%.4f", x$max_prox_share), "\n")
  cat("  max_prox_ratio: ", sprintf("%.4f", x$max_prox_ratio), "\n\n")

  cat("Mean-based metrics (supplementary):\n")
  cat("  prox_share:     ", sprintf("%.4f", x$prox_share), "\n")
  cat("  prox_ratio:     ", sprintf("%.4f", x$prox_ratio), "\n\n")

  cat("Forest diagnostics:\n")
  cat("  OOB error:      ", sprintf("%.4f", x$oob_error), "\n")
  cat("  Wilcoxon p:     ", sprintf("%.4f", x$wilcox_p),
      " (heuristic)\n\n")

  if (!is.null(x$null_distribution)) {
    cat("Null test:\n")
    cat("  p-value (share): ", sprintf("%.4f",
        x$null_distribution$null_pvalue), "\n")
    cat("  p-value (ratio): ", sprintf("%.4f",
        x$null_distribution$null_ratio_pvalue), "\n\n")
  }

  cat("Privacy: ", if (x$privacy_pass) "PASS" else "FAIL", "\n")

  if (!is.null(x$var_importance) && length(x$var_importance) > 0) {
    cat("\nTop variables by importance:\n")
    imp_sorted <- sort(x$var_importance, decreasing = TRUE)
    for (i in seq_len(min(10, length(imp_sorted)))) {
      cat("  ", names(imp_sorted)[i], ": ",
          sprintf("%.4f", imp_sorted[i]), "\n")
    }
  }

  cat("\nTop 5 outlier records (by max proximity to training): ",
      paste(x$top_outliers, collapse = ", "), "\n")
  invisible(x)
}

#' @export
plot.rf_privacy <- function(x, y = NULL, which = 1L, ...) {
  show <- rep(FALSE, 3)
  show[which] <- TRUE

  if (show[1]) {
    # Paired density: max proximity to training vs holdout
    df <- data.frame(
      proximity = c(x$max_prox_train, x$max_prox_holdout),
      group = rep(c("To training", "To holdout"),
                  each = length(x$max_prox_train))
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$proximity,
                                           fill = .data$group)) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::labs(title = "Max Proximity: Training vs Holdout",
                    x = "Max proximity", y = "Density",
                    fill = "Reference") +
      ggplot2::theme_minimal()
    print(p)
  }

  if (show[2]) {
    # Per-record difference histogram
    diffs <- x$max_prox_train - x$max_prox_holdout
    df <- data.frame(diff = diffs)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$diff)) +
      ggplot2::geom_histogram(bins = 30, fill = "steelblue",
                               color = "white") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                           color = "red") +
      ggplot2::labs(
        title = "Per-Record Max Proximity Difference",
        x = "max_prox(train) - max_prox(holdout)",
        y = "Count"
      ) +
      ggplot2::theme_minimal()
    print(p)
  }

  if (show[3]) {
    if (is.null(x$null_distribution)) {
      message("No null distribution available (null_test = FALSE)")
    } else {
      df <- data.frame(null = x$null_distribution$null_max_shares)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$null)) +
        ggplot2::geom_histogram(bins = 20, fill = "grey70",
                                 color = "white") +
        ggplot2::geom_vline(xintercept = x$max_prox_share,
                             color = "red", linewidth = 1.2) +
        ggplot2::labs(
          title = "Null Distribution of Max Proximity Share",
          x = "max_prox_share (permuted)",
          y = "Count",
          caption = sprintf("Observed = %.3f, p = %.3f",
                            x$max_prox_share,
                            x$null_distribution$null_pvalue)
        ) +
        ggplot2::theme_minimal()
      print(p)
    }
  }

  invisible(x)
}

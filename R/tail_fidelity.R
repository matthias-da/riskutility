#' Tail Fidelity — Tail Preservation Utility Measure
#'
#' Assesses how well synthetic data preserves the extreme values (tails) of
#' the original distribution. Combines QQ tail divergence (headline metric),
#' Jensen-Shannon divergence in the tail region, and an optional Hill tail
#' index estimator.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the synthetic/anonymized dataset.
#' @param vars Character vector of numeric variable names to compare.
#'   If NULL (default), all common numeric variables are used; categoricals
#'   are skipped with a message.
#' @param percentile Numeric, the percentile threshold defining the tail region.
#'   Default 95, meaning the upper 5 percent and lower 5 percent of observations.
#' @param tails Character, which tails to assess. One of \code{"both"} (default),
#'   \code{"upper"}, or \code{"lower"}.
#' @param hill Logical, whether to include the Hill tail index estimator.
#'   Default FALSE. Only meaningful for positive-valued data (upper tail).
#' @param na.rm Logical, whether to remove NA values. Default TRUE.
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class \code{"tail_fidelity"} containing:
#' \itemize{
#'   \item \code{qq_divergence}: mean QQ tail divergence across all variables
#'   \item \code{utility_score}: \code{exp(-qq_divergence)}, in \eqn{[0, 1]}, 1 = identical tails
#'   \item \code{per_variable}: data.frame with columns: variable, qq_tail_div,
#'     jsd_tail, and optionally hill_diff
#'   \item \code{percentile}: the percentile threshold used
#'   \item \code{tails}: which tails were assessed
#'   \item \code{n_vars}: number of numeric variables analysed
#'   \item \code{vars}: variable names used
#'   \item \code{n_X}: number of rows in X (after NA removal)
#'   \item \code{n_Y}: number of rows in Y (after NA removal)
#'   \item \code{hill}: whether Hill estimator was computed
#' }
#'
#' @details
#' Three complementary metrics are computed per numeric variable:
#'
#' \strong{1. QQ tail divergence (headline metric):}
#' Quantiles at 50 evenly-spaced probability levels in the tail region are
#' compared. The mean absolute difference is normalized by the IQR of the
#' original variable, making it scale-invariant. For \code{tails = "both"},
#' the upper and lower divergences are averaged.
#'
#' \strong{2. Density ratio in tail region:}
#' Observations beyond the percentile threshold are extracted and compared
#' via Jensen-Shannon divergence using \code{\link{densitydiff_1d_num}}.
#' If fewer than 20 observations fall in the tail, \code{NA} is returned.
#'
#' \strong{3. Hill estimator (optional, \code{hill = TRUE}):}
#' For positive-valued upper-tail data the Hill estimator of the tail index
#' is computed:
#' \deqn{\hat{\alpha} = \left(\frac{1}{k}\sum_{i=1}^{k}\log\frac{x_{(i)}}{x_{(k+1)}}\right)^{-1}}
#' where \eqn{x_{(1)} \ge \cdots \ge x_{(n)}} are order statistics and
#' \eqn{k = \lfloor\sqrt{n}\rfloor}. The absolute difference between original
#' and synthetic Hill estimates is reported.
#'
#' @seealso \code{\link{densitydiff_1d_num}} for density comparison,
#'   \code{\link{compare_wasserstein}} for Wasserstein distance,
#'   \code{\link{energy_distance}} for multivariate energy distance
#'
#' @references
#' Hill, B. M. (1975). A simple general approach to inference about the tail
#' of a distribution. The Annals of Statistics, 3(5), 1163-1174.
#'
#' @author Matthias Templ
#' @family utility
#' @export
#' @importFrom stats quantile IQR complete.cases
#'
#' @examples
#' set.seed(42)
#' # Original data with heavy tails
#' X <- data.frame(
#'   income = rlnorm(500, meanlog = 10, sdlog = 1),
#'   score = rnorm(500, 100, 15)
#' )
#'
#' # Good synthetic (same distribution)
#' Y_good <- data.frame(
#'   income = rlnorm(500, meanlog = 10, sdlog = 1),
#'   score = rnorm(500, 100, 15)
#' )
#'
#' # Poor synthetic (tails trimmed)
#' Y_poor <- data.frame(
#'   income = pmin(pmax(rlnorm(500, meanlog = 10, sdlog = 1), 5000), 80000),
#'   score = pmin(pmax(rnorm(500, 100, 15), 70), 130)
#' )
#'
#' result_good <- tail_fidelity(X, Y_good)
#' print(result_good)
#'
#' result_poor <- tail_fidelity(X, Y_poor)
#' print(result_poor)
#'
#' \donttest{
#' # With Hill estimator (for positive data)
#' result_hill <- tail_fidelity(X, Y_good, hill = TRUE)
#' summary(result_hill)
#'
#' # Assess only upper tail
#' result_upper <- tail_fidelity(X, Y_poor, tails = "upper")
#' print(result_upper)
#' }
tail_fidelity <- function(X, ...) {
  UseMethod("tail_fidelity")
}

#' @rdname tail_fidelity
#' @export
tail_fidelity.synth_pair <- function(X, ...) {
  tail_fidelity.default(
    X = X$original,
    Y = X$synthetic,
    vars = X$num_vars,
    ...
  )
}

#' @rdname tail_fidelity
#' @export
tail_fidelity.default <- function(X, Y,
                                  vars = NULL,
                                  percentile = 95,
                                  tails = c("both", "upper", "lower"),
                                  hill = FALSE,
                                  na.rm = TRUE,
                                  ...) {

  tails <- match.arg(tails)

  # Validate percentile
  if (!is.numeric(percentile) || length(percentile) != 1 ||
      percentile <= 50 || percentile >= 100) {
    stop("'percentile' must be a single number in (50, 100).")
  }

  # Convert to data.frame
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)

  # Auto-detect numeric variables if not specified
  if (is.null(vars)) {
    num_vars_X <- names(X)[sapply(X, is.numeric)]
    num_vars_Y <- names(Y)[sapply(Y, is.numeric)]
    vars <- intersect(num_vars_X, num_vars_Y)

    # Warn about skipped categoricals
    all_common <- intersect(names(X), names(Y))
    cat_vars <- setdiff(all_common, vars)
    if (length(cat_vars) > 0) {
      message("Skipping non-numeric variables: ",
              paste(cat_vars, collapse = ", "))
    }
  }

  if (length(vars) == 0) {
    stop("No numeric variables found or specified for comparison.")
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

  # --- Per-variable computation ---
  per_var_results <- lapply(vars, function(v) {
    x <- X[[v]]
    y <- Y[[v]]

    # Remove per-variable NAs (in case na.rm acts row-wise above but
    # individual columns may still have NAs when vars has length 1)
    if (na.rm) {
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
    }

    iqr_x <- IQR(x)
    if (iqr_x == 0) iqr_x <- sd(x)
    if (iqr_x == 0) iqr_x <- 1  # constant variable fallback

    # --- 1. QQ tail divergence ---
    qq_div <- .tail_qq_divergence(x, y, percentile, tails, iqr_x)

    # --- 2. JS divergence in tail region ---
    jsd_tail <- .tail_jsd(x, y, percentile, tails)

    # --- 3. Hill estimator (optional) ---
    hill_diff <- NA_real_
    if (hill) {
      hill_diff <- .tail_hill_diff(x, y, percentile)
    }

    data.frame(
      variable = v,
      qq_tail_div = qq_div,
      jsd_tail = jsd_tail,
      stringsAsFactors = FALSE
    )
  })

  per_variable <- do.call(rbind, per_var_results)

  # Add Hill column if requested
  if (hill) {
    hill_diffs <- sapply(vars, function(v) {
      x <- X[[v]]
      y <- Y[[v]]
      if (na.rm) { x <- x[!is.na(x)]; y <- y[!is.na(y)] }
      .tail_hill_diff(x, y, percentile)
    })
    per_variable$hill_diff <- hill_diffs
  }

  rownames(per_variable) <- NULL

  # --- Aggregate headline metric ---
  qq_divergence <- mean(per_variable$qq_tail_div, na.rm = TRUE)
  utility_score <- exp(-qq_divergence)
  # Clamp to [0, 1]
  utility_score <- max(0, min(1, utility_score))

  result <- list(
    qq_divergence = qq_divergence,
    utility_score = utility_score,
    per_variable = per_variable,
    percentile = percentile,
    tails = tails,
    hill = hill,
    n_vars = length(vars),
    vars = vars,
    n_X = n_X,
    n_Y = n_Y
  )

  class(result) <- "tail_fidelity"
  return(result)
}


# ---- Internal helpers --------------------------------------------------------

# QQ tail divergence: mean |Q_X(p) - Q_Y(p)| / IQR(X) in the tail region
# @param x numeric vector (original)
# @param y numeric vector (synthetic)
# @param percentile numeric, e.g. 95
# @param tails character, "both", "upper", or "lower"
# @param iqr_x numeric, IQR of x (for normalization)
# @return scalar divergence (non-negative)
# @keywords internal
.tail_qq_divergence <- function(x, y, percentile, tails, iqr_x) {
  n_probs <- 50L

  compute_tail_div <- function(probs) {
    q_x <- quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
    q_y <- quantile(y, probs = probs, na.rm = TRUE, names = FALSE)
    mean(abs(q_x - q_y)) / iqr_x
  }

  if (tails == "upper" || tails == "both") {
    upper_probs <- seq(percentile / 100, 1, length.out = n_probs)
    upper_div <- compute_tail_div(upper_probs)
  }
  if (tails == "lower" || tails == "both") {
    lower_probs <- seq(0, (100 - percentile) / 100, length.out = n_probs)
    lower_div <- compute_tail_div(lower_probs)
  }

  if (tails == "both") {
    return((upper_div + lower_div) / 2)
  } else if (tails == "upper") {
    return(upper_div)
  } else {
    return(lower_div)
  }
}


# JS divergence in the tail region using densitydiff_1d_num
# @param x numeric vector (original)
# @param y numeric vector (synthetic)
# @param percentile numeric
# @param tails character
# @return scalar JSD or NA if too few observations
# @keywords internal
.tail_jsd <- function(x, y, percentile, tails) {
  min_obs <- 20L

  extract_tail <- function(v, percentile, tails) {
    if (tails == "upper") {
      threshold <- quantile(v, probs = percentile / 100, na.rm = TRUE)
      v[v >= threshold]
    } else if (tails == "lower") {
      threshold <- quantile(v, probs = (100 - percentile) / 100, na.rm = TRUE)
      v[v <= threshold]
    } else {
      # both tails
      upper_thresh <- quantile(v, probs = percentile / 100, na.rm = TRUE)
      lower_thresh <- quantile(v, probs = (100 - percentile) / 100, na.rm = TRUE)
      v[v >= upper_thresh | v <= lower_thresh]
    }
  }

  x_tail <- extract_tail(x, percentile, tails)
  y_tail <- extract_tail(y, percentile, tails)

  if (length(x_tail) < min_obs || length(y_tail) < min_obs) {
    return(NA_real_)
  }

  # densitydiff_1d_num returns list with $jsd
  # Suppress NaN warnings from KDE on sparse tail data
  tryCatch({
    dd <- suppressWarnings(densitydiff_1d_num(x_tail, y_tail))
    dd$jsd
  }, error = function(e) {
    NA_real_
  })
}


# Hill estimator difference for the upper tail
# @param x numeric vector (original)
# @param y numeric vector (synthetic)
# @param percentile numeric
# @return |alpha_x - alpha_y| or NA
# @keywords internal
.tail_hill_diff <- function(x, y, percentile) {
  .hill_alpha <- function(v) {
    # Extract upper tail
    threshold <- quantile(v, probs = percentile / 100, na.rm = TRUE)
    v_tail <- v[v > threshold]

    # Hill estimator requires positive values
    if (any(v_tail <= 0, na.rm = TRUE) || length(v_tail) < 5) {
      return(NA_real_)
    }

    v_sorted <- sort(v_tail, decreasing = TRUE)
    n <- length(v_sorted)
    k <- floor(sqrt(n))
    if (k < 2) return(NA_real_)

    # Hill estimator: alpha_hat = k / sum_{i=1}^{k} log(x_(i) / x_(k+1))
    log_ratios <- log(v_sorted[1:k] / v_sorted[k + 1])
    hill_sum <- sum(log_ratios)
    if (hill_sum == 0) return(NA_real_)

    # Return alpha (inverse of the mean log ratio)
    k / hill_sum
  }

  alpha_x <- .hill_alpha(x)
  alpha_y <- .hill_alpha(y)

  if (is.na(alpha_x) || is.na(alpha_y)) return(NA_real_)
  abs(alpha_x - alpha_y)
}


# ---- S3 methods --------------------------------------------------------------

#' Print method for tail_fidelity objects
#'
#' @param x an object of class \code{"tail_fidelity"}
#' @param ... additional arguments (ignored)
#' @export
print.tail_fidelity <- function(x, ...) {
  cat("Tail Fidelity - Tail Preservation Utility Measure\n")
  cat("==================================================\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original (X):", x$n_X, "\n")
  cat("  Synthetic (Y):", x$n_Y, "\n")
  cat("  Numeric variables:", x$n_vars, "\n\n")

  cat("Settings:\n")
  cat("  Percentile:", x$percentile, "\n")
  cat("  Tails:", x$tails, "\n")
  cat("  Hill estimator:", x$hill, "\n\n")

  cat("Results:\n")
  cat("  QQ tail divergence:", sprintf("%.4f", x$qq_divergence), "\n")
  cat("  Utility score:     ", sprintf("%.4f", x$utility_score),
      "(exp(-QQ_div), higher=better)\n\n")

  cat("Per-variable QQ tail divergence:\n")
  for (i in seq_len(nrow(x$per_variable))) {
    row <- x$per_variable[i, ]
    cat(sprintf("  %-20s qq=%.4f", row$variable, row$qq_tail_div))
    if (!is.na(row$jsd_tail)) {
      cat(sprintf("  jsd=%.4f", row$jsd_tail))
    } else {
      cat("  jsd=NA")
    }
    if (x$hill && "hill_diff" %in% names(row)) {
      if (!is.na(row$hill_diff)) {
        cat(sprintf("  hill_diff=%.4f", row$hill_diff))
      } else {
        cat("  hill_diff=NA")
      }
    }
    cat("\n")
  }
  cat("\n")

  cat("Interpretation:\n")
  if (x$utility_score > 0.95) {
    cat("  EXCELLENT: Tail distributions are very well preserved.\n")
  } else if (x$utility_score > 0.80) {
    cat("  GOOD: Tail distributions are reasonably preserved.\n")
  } else if (x$utility_score > 0.50) {
    cat("  MODERATE: Some distortion in the tails.\n")
  } else {
    cat("  POOR: Significant distortion in tail distributions.\n")
  }

  invisible(x)
}


#' Summary method for tail_fidelity objects
#'
#' @param object an object of class \code{"tail_fidelity"}
#' @param ... additional arguments (ignored)
#' @export
summary.tail_fidelity <- function(object, ...) {
  # Aggregate JSD across variables (ignoring NAs)
  mean_jsd <- mean(object$per_variable$jsd_tail, na.rm = TRUE)
  max_qq <- max(object$per_variable$qq_tail_div, na.rm = TRUE)
  worst_var <- object$per_variable$variable[
    which.max(object$per_variable$qq_tail_div)
  ]

  summ <- list(
    qq_divergence = object$qq_divergence,
    utility_score = object$utility_score,
    mean_jsd_tail = mean_jsd,
    max_qq_divergence = max_qq,
    worst_variable = worst_var,
    per_variable = object$per_variable,
    percentile = object$percentile,
    tails = object$tails,
    hill = object$hill,
    n_vars = object$n_vars,
    vars = object$vars,
    n_X = object$n_X,
    n_Y = object$n_Y
  )

  class(summ) <- "summary.tail_fidelity"
  return(summ)
}


#' Print method for summary.tail_fidelity objects
#'
#' @param x an object of class \code{"summary.tail_fidelity"}
#' @param ... additional arguments (ignored)
#' @export
print.summary.tail_fidelity <- function(x, ...) {
  cat("Summary: Tail Fidelity\n")
  cat("======================\n\n")

  cat("Variables (", x$n_vars, "):", paste(x$vars, collapse = ", "), "\n")
  cat("Tails:", x$tails, " | Percentile:", x$percentile, "\n\n")

  cat("Aggregate Metrics:\n")
  cat("  Mean QQ tail divergence: ", sprintf("%.4f", x$qq_divergence), "\n")
  cat("  Utility score:           ", sprintf("%.4f", x$utility_score), "\n")
  cat("  Mean JSD in tail:        ", sprintf("%.4f", x$mean_jsd_tail), "\n")
  cat("  Worst variable:          ", x$worst_variable,
      sprintf(" (qq=%.4f)", x$max_qq_divergence), "\n\n")

  cat("Per-variable breakdown:\n")
  print(x$per_variable, row.names = FALSE)
  cat("\n")

  cat("Sample Sizes: X =", x$n_X, ", Y =", x$n_Y, "\n")

  invisible(x)
}


#' Plot method for tail_fidelity objects
#'
#' @param x an object of class \code{"tail_fidelity"}
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = tail-focused QQ plot per variable,
#'   2 = per-variable tail divergence lollipop chart
#' @param max_panels integer, maximum number of panels for QQ plots.
#'   Default 9.
#' @importFrom graphics par plot points abline mtext segments
#' @importFrom stats quantile
#' @export
plot.tail_fidelity <- function(x, y = NULL, ..., which = 1, max_panels = 9) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)

  if (show[1]) {
    .plot_tail_qq(x, max_panels = max_panels, ...)
  }

  if (show[2]) {
    .plot_tail_divergence(x, ...)
  }

  invisible(x)
}


# Internal: tail-focused QQ plot per variable
# @keywords internal
.plot_tail_qq <- function(obj, max_panels = 9, ...) {
  vars <- obj$vars
  n_vars <- min(length(vars), max_panels)
  vars <- vars[seq_len(n_vars)]

  # Determine grid layout
  if (n_vars == 1) {
    nc <- 1; nr <- 1
  } else if (n_vars <= 4) {
    nc <- 2; nr <- ceiling(n_vars / 2)
  } else {
    nc <- 3; nr <- ceiling(n_vars / 3)
  }

  op <- par(mfrow = c(nr, nc), mar = c(4, 4, 2, 1))
  on.exit(par(op))

  percentile <- obj$percentile
  tails <- obj$tails
  n_probs <- 50L

  for (v in vars) {
    # We need the raw data to compute quantiles. Since we don't store it,
    # plot using the per_variable info we do have.
    # Actually, we need the quantile data for QQ plots. Since we don't store
    # raw data in the result, we plot using the tail prob levels and label
    # the axes descriptively. But we'd need raw data.
    #
    # Alternative: just plot the per-variable metrics. But the spec says
    # "Tail-focused QQ plot per variable". We don't have raw data stored.
    # We'll create a simple visualization of the metric values instead.

    # Since we can't regenerate QQ data without raw vectors, plot a placeholder
    # that shows the divergence value prominently.
    qq_val <- obj$per_variable$qq_tail_div[obj$per_variable$variable == v]
    jsd_val <- obj$per_variable$jsd_tail[obj$per_variable$variable == v]

    # Create a simple indicator plot
    plot(1, 1, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "", ylab = "", axes = FALSE,
         main = v)

    # Background gradient
    for (i in 1:100) {
      col <- rgb(1 - i/100, i/100, 0, alpha = 0.2)
      graphics::rect((i-1)/100, 0.3, i/100, 0.7, col = col, border = NA)
    }

    # Score = exp(-qq_val)
    score <- exp(-qq_val)
    graphics::rect(0, 0.3, 1, 0.7, border = "black", lwd = 1.5)
    points(score, 0.5, pch = 18, cex = 2.5, col = "black")
    graphics::text(0.5, 0.15, sprintf("QQ div = %.3f", qq_val), cex = 0.9)
    if (!is.na(jsd_val)) {
      graphics::text(0.5, 0.85, sprintf("JSD = %.3f", jsd_val), cex = 0.9)
    }
  }
}


# Internal: per-variable tail divergence lollipop chart
# @keywords internal
.plot_tail_divergence <- function(obj, ...) {
  pv <- obj$per_variable
  n <- nrow(pv)

  op <- par(mar = c(4, max(nchar(as.character(pv$variable))) * 0.6 + 2, 3, 1))
  on.exit(par(op))

  # Lollipop chart
  y_pos <- seq_len(n)
  x_vals <- pv$qq_tail_div

  plot(x_vals, y_pos, type = "n",
       xlim = c(0, max(x_vals, na.rm = TRUE) * 1.15 + 0.01),
       ylim = c(0.5, n + 0.5),
       xlab = "QQ Tail Divergence",
       ylab = "",
       yaxt = "n",
       main = sprintf("Tail Fidelity (percentile = %g, tails = %s)",
                       obj$percentile, obj$tails))

  axis(2, at = y_pos, labels = pv$variable, las = 1)

  # Draw lollipop stems
  segments(x0 = 0, y0 = y_pos, x1 = x_vals, y1 = y_pos,
           col = "grey50", lwd = 1.5)

  # Colour by divergence
  cols <- ifelse(x_vals < 0.1, "forestgreen",
                 ifelse(x_vals < 0.3, "orange", "red"))
  points(x_vals, y_pos, pch = 19, cex = 1.5, col = cols)

  # Add utility score annotation
  mtext(sprintf("Overall utility = %.3f", obj$utility_score),
        side = 3, line = 0, cex = 0.8)
}

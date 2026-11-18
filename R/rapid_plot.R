# =============================================================================
# Internal plot helpers
# =============================================================================

#' Plot QI Variable Importance Barplot
#'
#' @param x A rapid object (must have model stored via store_model=TRUE)
#' @param top_n Maximum number of variables to display
#' @return Invisible NULL
#' @noRd
.plot_qi_importance <- function(x, top_n = 15) {
  .rapid_require(x, c("version", "model"))
  learner <- if (!is.null(x$attack)) x$attack$learner else x$predictor$learner
  imp <- .learner_importance(learner, x$model)

  if (is.null(imp) || length(imp) == 0) {
    message("Variable importance not available for model type: ", x$model_type)
    return(invisible(NULL))
  }

  # Sort and limit
  imp <- sort(imp, decreasing = TRUE)
  if (length(imp) > top_n) {
    imp <- imp[1:top_n]
  }

  # Normalize to 0-1 scale
  imp <- imp / max(imp)

  # Reverse for bottom-to-top display (largest at top)
  imp <- rev(imp)
  n <- length(imp)
  y_pos <- seq_len(n)

  # Adjust margins for long labels
  max_label_len <- max(nchar(names(imp)))
  left_margin <- min(max(8, max_label_len * 0.6), 20)

  old_par <- par(mar = c(5, left_margin, 4, 2))
  on.exit(par(old_par))

  # Lollipop chart: segments from 0 + points
  plot(imp, y_pos, type = "n",
       xlim = c(0, 1.05),
       ylim = c(0.5, n + 0.5),
       yaxt = "n",
       xlab = "Relative Importance",
       ylab = "",
       main = paste0("QI Variable Importance (", x$model_type, ")"))
  segments(0, y_pos, imp, y_pos, col = "gray60", lwd = 1.5)
  points(imp, y_pos, pch = 19, col = "steelblue", cex = 1.3)
  axis(2, at = y_pos, labels = names(imp), las = 1, cex.axis = 0.85)

  invisible(NULL)
}

#' Discretize QI data for attribution GLM
#'
#' Converts numeric and high-cardinality factor variables to quantile bins
#' to prevent quasi-complete separation in logistic regression.
#' @param qi_data data frame of QI columns
#' @param max_levels integer, factors with more levels are rebinned (default 10)
#' @param n_bins integer, number of quantile bins for numeric variables (default 4)
#' @return list with qi_data (discretized), level_map, contrast_list
#' @noRd
.prepare_attribution_data <- function(qi_data, at_risk = NULL, max_levels = 10,
                                      n_bins = 4) {
  # Adaptive binning: reduce bins when at-risk count is low to prevent separation
  # Aim for >= 5 expected at-risk records per bin
  if (!is.null(at_risk)) {
    n_at_risk <- sum(at_risk)
    n_minority <- min(n_at_risk, length(at_risk) - n_at_risk)
    n_bins <- max(2L, min(n_bins, as.integer(floor(n_minority / 5))))
  }

  level_map <- list()
  contrast_list <- list()

  for (qi in names(qi_data)) {
    col <- qi_data[[qi]]

    if (is.numeric(col)) {
      # Discretize numeric into quantile bins
      breaks <- unique(stats::quantile(col, probs = seq(0, 1, length.out = n_bins + 1),
                                       na.rm = TRUE))
      if (length(breaks) <= 2) {
        qi_data[[qi]] <- as.factor(col)
      } else {
        qi_data[[qi]] <- cut(col, breaks = breaks, include.lowest = TRUE,
                             dig.lab = 3)
      }
    } else {
      qi_data[[qi]] <- as.factor(qi_data[[qi]])
    }

    # Rebin factors with too many levels
    if (is.factor(qi_data[[qi]]) && nlevels(qi_data[[qi]]) > max_levels) {
      col_num <- as.numeric(qi_data[[qi]])
      reb <- max(2L, min(n_bins, nlevels(qi_data[[qi]])))
      breaks <- unique(stats::quantile(col_num,
                                       probs = seq(0, 1, length.out = reb + 1),
                                       na.rm = TRUE))
      qi_data[[qi]] <- cut(col_num, breaks = breaks, include.lowest = TRUE,
                           labels = paste0("G", seq_len(length(breaks) - 1)))
    }

    # Set up sum-to-zero contrasts and label mapping
    lvls <- levels(qi_data[[qi]])
    n_levels <- length(lvls)
    if (n_levels > 1) {
      contrast_list[[qi]] <- stats::contr.sum(n_levels)
      for (i in seq_len(n_levels - 1)) {
        coef_name <- paste0(qi, i)
        level_map[[coef_name]] <- paste0(qi, ": ", lvls[i])
      }
      level_map[[paste0(qi, "_ref")]] <- paste0(qi, ": ", lvls[n_levels])
    }
  }

  list(qi_data = qi_data, level_map = level_map, contrast_list = contrast_list)
}

#' Fit attribution GLM and extract coefficients with reference levels
#' @noRd
.fit_attribution_glm <- function(formula, qi_data, at_risk, contrast_list, level_map) {
  model_data <- cbind(at_risk = as.numeric(at_risk), qi_data)

  fit <- tryCatch({
    suppressWarnings(
      stats::glm(formula, data = model_data, family = stats::binomial(),
                 contrasts = contrast_list)
    )
  }, error = function(e) {
    .rapid_warn(c("Could not fit the attribution model: ", conditionMessage(e)),
                "attribution_model")
    return(NULL)
  })

  if (is.null(fit)) return(NULL)

  coefs <- summary(fit)$coefficients
  V <- stats::vcov(fit)
  coef_df <- data.frame(
    term = rownames(coefs),
    estimate = coefs[, "Estimate"],
    se = coefs[, "Std. Error"]
  )

  # Remove intercept
  coef_df <- coef_df[coef_df$term != "(Intercept)", , drop = FALSE]
  if (nrow(coef_df) == 0) return(NULL)

  # Add implied reference-level coefficients (sum-to-zero: beta_ref = -sum(others))
  for (qi in names(contrast_list)) {
    lvls <- levels(qi_data[[qi]])
    n_levels <- length(lvls)
    if (n_levels < 2) next

    est_terms <- paste0(qi, seq_len(n_levels - 1))
    est_terms <- est_terms[est_terms %in% coef_df$term]
    if (length(est_terms) == 0) next

    ref_est <- -sum(coef_df$estimate[coef_df$term %in% est_terms])

    V_sub <- V[est_terms, est_terms, drop = FALSE]
    ref_var <- sum(V_sub)
    ref_se <- if (ref_var > 0) sqrt(ref_var) else NA_real_

    ref_label <- level_map[[paste0(qi, "_ref")]]
    if (is.null(ref_label)) ref_label <- paste0(qi, ": ", lvls[n_levels])

    coef_df <- rbind(coef_df, data.frame(
      term = paste0(qi, "_ref"),
      estimate = ref_est,
      se = ref_se
    ))
  }

  # Improve term labels
  coef_df$label <- vapply(coef_df$term, function(t) {
    if (t %in% names(level_map)) return(level_map[[t]])
    t_clean <- t
    for (qi in names(qi_data)) {
      if (is.factor(qi_data[[qi]])) {
        lvls <- levels(qi_data[[qi]])
        for (i in seq_along(lvls)) {
          t_clean <- gsub(paste0("\\b", qi, i, "\\b"), lvls[i], t_clean)
        }
      }
    }
    gsub(":", " x ", t_clean)
  }, character(1))

  coef_df
}

#' Draw forest-plot style dot chart from coefficient data frame
#' @noRd
.draw_forest_plot <- function(coef_df, main = "QI Attribution (Effect Coding)",
                              xlab = "Log-odds contribution to at-risk") {
  n <- nrow(coef_df)
  y_pos <- seq_len(n)

  colors <- ifelse(coef_df$estimate > 0, "firebrick", "steelblue")

  valid_se <- is.finite(coef_df$se)
  ci_lo <- coef_df$estimate - 1.96 * coef_df$se
  ci_hi <- coef_df$estimate + 1.96 * coef_df$se

  # Base x-axis range on point estimates; include CIs only if they are reasonable
  # (SE < 10 is a practical threshold; larger indicates quasi-separation)
  reasonable <- valid_se & coef_df$se < 10
  pts_range <- range(coef_df$estimate, na.rm = TRUE)
  if (any(reasonable)) {
    x_range <- range(c(pts_range, ci_lo[reasonable], ci_hi[reasonable]))
  } else {
    x_range <- pts_range
  }
  x_pad <- max(diff(x_range) * 0.15, 0.05)
  x_lim <- c(x_range[1] - x_pad, x_range[2] + x_pad)

  max_label_len <- max(nchar(coef_df$label))
  left_margin <- min(max(8, max_label_len * 0.5), 20)

  old_par <- par(mar = c(5, left_margin, 4, 2))
  on.exit(par(old_par))

  plot(coef_df$estimate, y_pos, type = "n",
       xlim = x_lim, ylim = c(0.5, n + 0.5),
       yaxt = "n", xlab = xlab, ylab = "", main = main)
  abline(v = 0, lty = 2, col = "gray40")

  # Draw CIs: only for reasonable SEs, truncated to plot limits
  if (any(reasonable)) {
    ci_lo_draw <- pmax(ci_lo[reasonable], x_lim[1])
    ci_hi_draw <- pmin(ci_hi[reasonable], x_lim[2])
    ci_width <- abs(ci_hi_draw - ci_lo_draw)
    draw <- ci_width > .Machine$double.eps * 100
    if (any(draw)) {
      suppressWarnings(
        arrows(ci_lo_draw[draw], y_pos[reasonable][draw],
               ci_hi_draw[draw], y_pos[reasonable][draw],
               angle = 90, code = 3, length = 0.03, col = "gray50")
      )
    }
  }

  points(coef_df$estimate, y_pos, pch = 19, col = colors, cex = 1.2)
  axis(2, at = y_pos, labels = coef_df$label, las = 1, cex.axis = 0.8)
}

#' Plot QI Attribution Analysis using Effect Coding
#'
#' @param x A rapid object
#' @param formula Optional formula for the attribution model
#' @param top_n Maximum number of terms to display
#' @return Invisible NULL
#' @noRd
.plot_qi_attribution <- function(x, formula = NULL, top_n = 20) {
  records <- x$records

  if (is.null(records) || !("at_risk" %in% names(records))) {
    .rapid_abort("Records with at_risk column required.", "object")
  }

  if (all(records$at_risk) || !any(records$at_risk)) {
    message("Cannot fit attribution model: all records have the same at_risk status.")
    return(invisible(NULL))
  }

  key_vars <- x$key_vars

  # Discretize numeric / high-cardinality variables
  prep <- .prepare_attribution_data(records[, key_vars, drop = FALSE],
                                     at_risk = records$at_risk)

  # Build formula: main effects only
  if (is.null(formula)) {
    formula <- stats::as.formula(paste("at_risk ~", paste(key_vars, collapse = " + ")))
  }

  coef_df <- .fit_attribution_glm(formula, prep$qi_data, records$at_risk,
                                  prep$contrast_list, prep$level_map)
  if (is.null(coef_df) || nrow(coef_df) == 0) {
    message("No coefficients to display.")
    return(invisible(NULL))
  }

  # Limit and sort
  if (nrow(coef_df) > top_n) {
    coef_df <- coef_df[order(abs(coef_df$estimate), decreasing = TRUE), ]
    coef_df <- coef_df[1:top_n, ]
  }
  coef_df <- coef_df[order(coef_df$estimate), ]

  .draw_forest_plot(coef_df, main = "QI Attribution (Effect Coding)")

  invisible(NULL)
}

#' Plot QI Interaction Effects
#'
#' Fits a logistic regression with all pairwise interactions between QI variables
#' and shows only the interaction coefficients.
#' @param x A rapid object
#' @param top_n Maximum number of interaction terms to display
#' @return Invisible NULL
#' @noRd
.plot_qi_interactions <- function(x, top_n = 20) {
  records <- x$records

  if (is.null(records) || !("at_risk" %in% names(records))) {
    .rapid_abort("Records with at_risk column required.", "object")
  }

  if (all(records$at_risk) || !any(records$at_risk)) {
    message("Cannot fit interaction model: all records have the same at_risk status.")
    return(invisible(NULL))
  }

  key_vars <- x$key_vars
  if (length(key_vars) < 2) {
    message("Interaction plot requires at least 2 key variables.")
    return(invisible(NULL))
  }

  # Discretize
  prep <- .prepare_attribution_data(records[, key_vars, drop = FALSE],
                                     at_risk = records$at_risk)

  # Formula with all pairwise interactions
  pairs <- utils::combn(key_vars, 2, simplify = FALSE)
  interaction_terms <- vapply(pairs, function(p) paste(p, collapse = ":"), character(1))
  formula <- stats::as.formula(
    paste("at_risk ~", paste(c(key_vars, interaction_terms), collapse = " + "))
  )

  coef_df <- .fit_attribution_glm(formula, prep$qi_data, records$at_risk,
                                  prep$contrast_list, prep$level_map)
  if (is.null(coef_df) || nrow(coef_df) == 0) {
    message("No coefficients to display.")
    return(invisible(NULL))
  }

  # Keep only interaction terms (those containing " x " in the label, or ":" in the term)
  is_interaction <- grepl(":", coef_df$term)
  coef_df <- coef_df[is_interaction, , drop = FALSE]

  if (nrow(coef_df) == 0) {
    message("No interaction coefficients to display.")
    return(invisible(NULL))
  }

  if (nrow(coef_df) > top_n) {
    coef_df <- coef_df[order(abs(coef_df$estimate), decreasing = TRUE), ]
    coef_df <- coef_df[1:top_n, ]
  }
  coef_df <- coef_df[order(coef_df$estimate), ]

  .draw_forest_plot(coef_df, main = "QI Interaction Effects (Effect Coding)")

  invisible(NULL)
}

#' Plot Threshold Sensitivity Curve
#'
#' @param x A rapid object
#' @param tau_range Thresholds for rules that flag scores above the threshold
#' @param epsilon_range Tolerances for the point rule
#' @return Invisible NULL
#' @noRd
.plot_threshold_sensitivity <- function(x, tau_range = seq(0, 1, by = 0.05),
                                        epsilon_range = seq(0, 1, by = 0.05)) {
  scores <- x$score[is.finite(x$score)]

  if (length(scores) == 0) {
    message("No records available for threshold sensitivity plot.")
    return(invisible(NULL))
  }

  if (x$rule$direction == "below") {
    # The score is an error; a record is at risk below the tolerance
    thresholds <- epsilon_range
    if (identical(epsilon_range, seq(0, 1, by = 0.05)) && max(scores) > 1) {
      # Errors in the units of the target (absolute metric)
      thresholds <- seq(0, max(scores), length.out = 21)
    }
    rapid_values <- vapply(thresholds, function(eps) mean(scores < eps), numeric(1))
    xlab_text <- "Tolerance (epsilon)"
    label <- "epsilon"
  } else {
    thresholds <- tau_range
    rapid_values <- vapply(thresholds, function(tau) mean(scores > tau), numeric(1))
    xlab_text <- "Threshold (tau)"
    label <- "tau"
  }

  plot(thresholds, rapid_values,
       type = "l", lwd = 2, col = "blue",
       xlab = xlab_text,
       ylab = "RAPID (Proportion at Risk)",
       main = "Threshold Sensitivity Curve",
       ylim = c(0, 1))

  points(x$threshold, x$rapid, pch = 19, col = "red", cex = 1.5)
  text(x$threshold, x$rapid,
       paste0(label, " = ", round(x$threshold, 3), "\nRAPID = ", round(x$rapid, 3)),
       pos = 4, col = "red")
  graphics::grid()

  invisible(NULL)
}

#' Plot method for rapid objects
#'
#' Draws one or several of six plot types for the result of
#' \code{\link{rapid}()}. The other methods for \code{rapid} objects are
#' documented in \code{\link{rapid-methods}}.
#'
#' @param x an object of class \code{"rapid"}
#' @param y not used
#' @param ... additional arguments passed to the plotting functions
#' @param which integer vector of plot types, drawn in the order given
#'   (default 1):
#'   \describe{
#'     \item{1}{Distribution of the scores with the threshold.}
#'     \item{2}{For the band and \code{RCS_marginal} rules the confidence
#'       \eqn{g} against the baseline \eqn{b}. For the other rules of a
#'       continuous target the prediction against the true value.}
#'     \item{3}{RAPID against the threshold.}
#'     \item{4}{Importance of the quasi-identifiers in the attacker model (needs
#'       \code{store_model = TRUE}).}
#'     \item{5}{Main effects of the quasi-identifiers on the flag.}
#'     \item{6}{Pairwise interactions of the quasi-identifiers on the flag (needs
#'       at least two quasi-identifiers).}
#'   }
#' @param xlim x-axis limits of plots 1 and 2: \code{"full"} (default, the
#'   range of the data), \code{"zero_to_one"}, \code{"threshold_to_one"} (the
#'   region of the records at risk), or a numeric vector of length 2
#' @param annotate logical, if \code{TRUE} (default) add the percentage at risk
#' @param facet_by character vector of variables of the record table. Plot 1
#'   is then drawn with one panel per subgroup (ggplot2).
#' @param bins integer, number of histogram bins (default 30)
#' @param tau_range numeric vector of thresholds for plot 3, for rules that flag
#'   scores above the threshold (default \code{seq(0, 1, by = 0.05)})
#' @param epsilon_range numeric vector of tolerances for plot 3 under the point
#'   rule (default \code{seq(0, 1, by = 0.05)}). For the absolute metric the
#'   range of the errors is used.
#' @param formula optional formula of the attribution model of plot 5, for
#'   example \code{at_risk ~ age * sex}. The default uses the main effects of
#'   all quasi-identifiers.
#' @param top_n integer, the largest number of variables or terms shown in plots
#'   4, 5 and 6 (default 15)
#'
#' @details
#' \strong{Gain rules.} For the band rule and \code{RCS_marginal}, types 1 to 3
#' are drawn with ggplot2. Type 1 is the histogram of the gain \eqn{r} with the
#' threshold. When the predictions are stored, a grey outline adds the
#' permutation null, the average histogram of one permutation. Type 2 shows the
#' confidence \eqn{g} against the baseline \eqn{b} with the lines \eqn{g = b}
#' and \eqn{g = b + \tau (1 - b)}{g = b + tau (1 - b)}. Records above the dashed
#' line are at risk, records near the diagonal gain nothing over the marginal
#' distribution. Type 3 shows RAPID against the threshold together with the
#' chance curve, the share of records flagged under permuted true values. The
#' permutations are drawn once per call from the stored predictions and the
#' \code{n_perm} of the fit, starting from the seed of the fit when there is
#' one. The chance
#' curve then passes through \code{x$chance} at the object's threshold, as
#' long as the caller has not changed \code{RNGkind()} since the fit. The random
#' stream of the caller is restored afterwards. Without stored predictions or
#' with \code{n_perm = 0}, the outline and the chance curve are left out.
#'
#' \strong{Other rules.} The point rule, \code{RCS_conditional} and
#' \code{NCE} draw types 1 to 3 with base graphics. Type 1 is the histogram of
#' the score with the threshold as a red dashed line. Under the point rule the
#' score is the error and low values are at risk. Type 2 (continuous target
#' only, nothing is drawn for a categorical one) shows the prediction against
#' the true value, coloured by the flag.
#' Type 3 marks the threshold of the object with a red point. \code{which =
#' 1:2} draws both panels side by side.
#'
#' \strong{Quasi-identifiers.} Type 4 is a lollipop chart of the variable
#' importance of the attacker (impurity for rf, the importance of rpart for
#' cart, gain for gbm, and the absolute coefficients scaled to sum 1 for lm and
#' logit). Type 5
#' fits a logistic regression of the flag on the quasi-identifiers with
#' sum-to-zero contrasts and shows every level, including the reference level,
#' with 95 percent intervals. Numeric quasi-identifiers are cut into quartile
#' bins and factors with many levels are grouped, which avoids separation.
#' Positive coefficients (red) raise the chance of being flagged, negative ones
#' (blue) lower it. Type 6 fits the same model with all pairwise interactions
#' and shows the interaction terms only.
#'
#' @return Invisibly, the ggplot2 object for one ggplot2 type, a list named by the type
#'   numbers (\code{"1"}, \code{"3"}, ...) for several types, with \code{NULL} entries for
#'   the types drawn with base graphics.
#'
#' @seealso \code{\link{rapid}}, \code{\link{rapid-methods}}
#' @family rapid
#' @examples
#' set.seed(42)
#' n <- 200
#' X <- data.frame(
#'   age = sample(20:60, n, replace = TRUE),
#'   gender = factor(sample(c("M", "F"), n, replace = TRUE)),
#'   education = factor(sample(c("High", "Medium", "Low"), n, replace = TRUE))
#' )
#' X$income <- round(exp(9 + 0.03 * X$age + 0.3 * (X$gender == "M") + rnorm(n, 0, 0.1)))
#' X$health <- factor(ifelse(
#'   X$age > 40 & X$gender == "M",
#'   sample(c("Poor", "Fair", "Good"), n, replace = TRUE, prob = c(0.6, 0.3, 0.1)),
#'   sample(c("Poor", "Fair", "Good"), n, replace = TRUE, prob = c(0.1, 0.3, 0.6))
#' ))
#' Y <- X
#' Y$income <- round(Y$income * exp(rnorm(n, 0, 0.05)))
#' idx <- sample(n, n * 0.3)
#' Y$health[idx] <- sample(levels(X$health), length(idx), replace = TRUE)
#'
#' # Continuous target, band rule: score histogram, g against b, chance curve
#' r <- rapid(X, Y, key_vars = c("age", "gender", "education"),
#'            target_var = "income", model_type = "lm", store_model = TRUE, seed = 1)
#' plot(r, which = 1:3)
#' plot(r, which = 4)
#'
#' # The point rule draws with base graphics
#' plot(update(r, num_eval_method = "point"), which = 1:2)
#'
#' # Categorical target: the quasi-identifiers behind the flags
#' if (requireNamespace("rpart", quietly = TRUE)) {
#'   rc <- rapid(X, Y, key_vars = c("age", "gender", "education"),
#'               target_var = "health", model_type = "cart", store_model = TRUE,
#'               seed = 1)
#'   plot(rc, which = 4:6)
#'   plot(rc, facet_by = "gender")
#' }
#'
#' @export
#' @importFrom graphics hist abline par legend points text mtext arrows grid segments axis
#' @importFrom stats glm binomial contr.sum as.formula coef vcov quantile
#' @importFrom utils combn modifyList
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline annotate facet_wrap
#'   labs theme_minimal theme element_text scale_fill_manual after_stat
plot.rapid <- function(x, y = NULL, ...,
                       which = 1,
                       xlim = "full",
                       annotate = TRUE,
                       facet_by = NULL,
                       bins = 30,
                       tau_range = seq(0, 1, by = 0.05),
                       epsilon_range = seq(0, 1, by = 0.05),
                       formula = NULL,
                       top_n = 15) {

  .rapid_require(x, "version")
  which <- .rapid_plot_which(which)
  gain <- isTRUE(x$rule$gain)

  # The null scores are drawn once per call and serve the outline of type 1 and
  # the chance curve of type 3
  null <- NULL
  drawn <- FALSE
  get_null <- function() {
    if (!drawn) {
      null <<- .rapid_plot_null(x)
      drawn <<- TRUE
    }
    null
  }

  # Two base-graphics panels side by side, as before
  base12 <- !gain && is.null(facet_by) && identical(which, 1:2)
  if (base12) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  res <- vector("list", length(which))
  names(res) <- which
  for (k in seq_along(which)) {
    w <- which[k]
    out <- switch(
      w,
      if (!is.null(facet_by)) {
        .plot_score_facets(x, facet_by, xlim, bins)
      } else if (gain) {
        .plot_gain_histogram(x, get_null(), xlim, annotate, bins)
      } else {
        .plot_score_histogram(x, xlim, annotate, bins)
      },
      if (gain) .plot_gain_scatter(x, annotate) else .plot_prediction_scatter(x, annotate),
      if (gain) {
        .plot_gain_sensitivity(x, get_null(), tau_range)
      } else {
        .plot_threshold_sensitivity(x, tau_range = tau_range, epsilon_range = epsilon_range)
      },
      .plot_qi_importance(x, top_n = top_n),
      .plot_qi_attribution(x, formula = formula, top_n = top_n),
      .plot_qi_interactions(x, top_n = top_n)
    )
    if (inherits(out, "ggplot")) {
      print(out)
      res[[k]] <- out
    }
  }
  invisible(if (length(which) == 1) res[[1]] else res)
}

#' The plot types asked for, checked
#' @noRd
.rapid_plot_which <- function(which) {
  if (!is.numeric(which) || length(which) == 0 || anyNA(which) ||
      any(which != round(which)) || any(which < 1 | which > 6)) {
    .rapid_abort("'which' must be one or more of the plot types 1 to 6.", "input")
  }
  unique(as.integer(which))
}

#' Null scores for the plots, drawn from the seed of the fit
#'
#' Needs the stored predictions and `n_perm > 0`; `NULL` otherwise. The draw
#' starts from the seed of the fit (when there is one), so that the chance
#' curve at the threshold of the object equals `x$chance`. The random stream of
#' the caller is restored.
#' @return Matrix from `.rapid_null()` or `NULL`.
#' @noRd
.rapid_plot_null <- function(x) {
  if (is.null(x$attack) || !isTRUE(x$n_perm > 0)) return(NULL)
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = env, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  })
  if (!is.null(x$seed)) set.seed(x$seed)
  .rapid_null(x$attack, x$rule, x$n_perm)
}

#' Finite scores, records and x-axis limits shared by the type 1 plots
#' @return A list, or `NULL` (with a message) without finite scores.
#' @noRd
.rapid_plot_prep <- function(x, xlim) {
  records <- x$records
  keep <- is.finite(x$score)
  if (!any(keep)) {
    message("No finite scores to plot.")
    return(NULL)
  }
  if (!all(keep)) records <- records[keep, , drop = FALSE]
  scores <- x$score[keep]
  below <- x$rule$direction == "below"
  score_range <- range(scores)
  full_range <- if (below) c(0, score_range[2]) else score_range
  xlim_values <- full_range
  if (is.character(xlim)) {
    xlim_values <- switch(xlim,
      "full" = full_range,
      "zero_to_one" = c(0, 1),
      "threshold_to_one" = c(x$threshold, 1),
      full_range
    )
  } else if (is.numeric(xlim) && length(xlim) == 2) {
    xlim_values <- xlim
  }
  list(records = records, keep = keep, scores = scores, below = below,
       score_range = score_range, xlim = xlim_values,
       label = x$rule$score_label)
}

# Colours of the ggplot2 histograms, as in the faceted plot
.rapid_plot_cols <- c("Not at risk" = "steelblue", "At risk" = "coral")

#' Histogram of the gain with the permutation null as an outline
#' @param null Matrix of null scores or `NULL`.
#' @return A ggplot object (or `NULL` without finite scores).
#' @noRd
.plot_gain_histogram <- function(x, null, xlim, annotate, bins) {
  pr <- .rapid_plot_prep(x, xlim)
  if (is.null(pr)) return(NULL)
  vals <- pr$scores
  if (!is.null(null)) null <- null[is.finite(null)]
  lo <- min(pr$xlim[1], pr$score_range[1], null)
  hi <- max(pr$xlim[2], pr$score_range[2], null)
  if (hi <= lo) hi <- lo + 1
  brk <- seq(lo, hi, length.out = bins + 1)

  d <- data.frame(
    r = vals,
    status = factor(ifelse(x$at_risk[pr$keep], "At risk", "Not at risk"),
                    levels = c("Not at risk", "At risk"))
  )
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$r)) +
    ggplot2::geom_histogram(ggplot2::aes(fill = .data$status), breaks = brk,
                            colour = "white", alpha = 0.8) +
    ggplot2::geom_vline(xintercept = x$threshold, colour = "red", linetype = "dashed",
                        linewidth = 1)
  if (length(null) > 0) {
    cnt <- tabulate(findInterval(null, brk, rightmost.closed = TRUE), nbins = bins) /
      x$n_perm
    null_label <- paste0("Permutation null (mean of ", x$n_perm, " permutations)")
    steps <- data.frame(r = brk, n = c(cnt, cnt[bins]), line = null_label)
    p <- p +
      ggplot2::geom_step(data = steps, ggplot2::aes(x = .data$r, y = .data$n,
                                                    colour = .data$line),
                         inherit.aes = FALSE, linewidth = 0.7) +
      ggplot2::scale_colour_manual(values = stats::setNames("grey30", null_label),
                                   name = NULL)
  }
  subtitle <- if (annotate) {
    paste0(sprintf("%.1f%% at risk (threshold = %s)", x$pct_at_risk, format(x$threshold)),
           if (is.finite(x$chance)) sprintf(", chance level %.1f%%", 100 * x$chance))
  }
  p +
    ggplot2::scale_fill_manual(values = .rapid_plot_cols, name = "Status") +
    ggplot2::coord_cartesian(xlim = pr$xlim) +
    ggplot2::labs(title = paste(pr$label, "Distribution"), subtitle = subtitle,
                  x = pr$label, y = "Count") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#' Confidence against baseline with the line of the threshold
#' @return A ggplot object (or `NULL` without finite values).
#' @noRd
.plot_gain_scatter <- function(x, annotate) {
  d <- data.frame(
    b = .rapid_col(x, "b"), g = .rapid_col(x, "g"),
    status = factor(ifelse(x$at_risk, "At risk", "Not at risk"),
                    levels = c("Not at risk", "At risk"))
  )
  d <- d[is.finite(d$b) & is.finite(d$g), , drop = FALSE]
  if (nrow(d) == 0) {
    message("No finite values to plot.")
    return(NULL)
  }
  tau <- x$threshold
  lines <- data.frame(
    intercept = c(0, tau), slope = c(1, 1 - tau),
    line = c("g = b", paste0("g = b + ", format(tau), " (1 - b)"))
  )
  lines$line <- factor(lines$line, levels = lines$line)
  # Plot the records at risk last, so that they are not hidden
  d <- d[order(d$status), , drop = FALSE]
  ggplot2::ggplot(d, ggplot2::aes(x = .data$b, y = .data$g)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$status), alpha = 0.6, size = 1.6) +
    ggplot2::geom_abline(data = lines,
                         ggplot2::aes(intercept = .data$intercept, slope = .data$slope,
                                      linetype = .data$line),
                         colour = "grey30", linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = .rapid_plot_cols, name = "Status") +
    ggplot2::scale_linetype_manual(values = c("dotted", "dashed"), name = NULL) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = "Confidence against baseline",
      subtitle = if (annotate) {
        sprintf("%.1f%% at risk (n = %d / %d)", x$pct_at_risk, x$n_at_risk, x$n_total)
      },
      x = if (x$is_categorical) "Baseline b (share of the true class)" else
        "Baseline b (chance of the band)",
      y = "Confidence g (attacker)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#' RAPID over a threshold grid with the chance curve
#' @param null Matrix of null scores or `NULL`.
#' @return A ggplot object (or `NULL` without records).
#' @noRd
.plot_gain_sensitivity <- function(x, null, tau_range) {
  if (!is.numeric(tau_range) || length(tau_range) == 0 || anyNA(tau_range)) {
    .rapid_abort("'tau_range' must be a numeric vector.", "input")
  }
  tau <- sort(unique(tau_range))
  rate <- function(tt) {
    rule <- .rule_set_threshold(x$rule, tt)
    c(mean(.rule_flag(rule, x$score)),
      if (!is.null(null)) mean(.rapid_chance_rates(rule, null)))
  }
  vals <- vapply(tau, function(tt) rate(tt), numeric(1 + !is.null(null)))
  vals <- matrix(vals, nrow = 1 + !is.null(null))
  d <- data.frame(tau = tau, value = vals[1, ], series = "RAPID")
  cols <- c(RAPID = "firebrick")
  if (!is.null(null)) {
    d <- rbind(d, data.frame(tau = tau, value = vals[2, ], series = "Chance level"))
    cols <- c(cols, "Chance level" = "grey40")
  }
  pt <- data.frame(tau = x$threshold, value = x$rapid)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$tau, y = .data$value, colour = .data$series)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(data = pt, ggplot2::aes(x = .data$tau, y = .data$value),
                        inherit.aes = FALSE, colour = "firebrick", size = 2.5) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Threshold Sensitivity Curve", x = "Threshold (tau)",
      y = "Share of records at risk",
      subtitle = sprintf("Point: tau = %s, RAPID = %.3f", format(x$threshold), x$rapid)
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#' Faceted histogram of the score (ggplot2, all rules)
#' @noRd
.plot_score_facets <- function(x, facet_by, xlim, bins) {
  pr <- .rapid_plot_prep(x, xlim)
  if (is.null(pr)) return(NULL)
  records <- pr$records
  score_var <- x$record_cols[[x$rule$score_var]]
  missing_vars <- setdiff(facet_by, names(records))
  if (length(missing_vars) > 0) {
    .rapid_abort(c("Facet variables not found in records: ",
                   paste(missing_vars, collapse = ", "), "."), "input")
  }
  if (length(facet_by) == 1) {
    records$facet_var <- records[[facet_by]]
  } else {
    records$facet_var <- interaction(records[, facet_by, drop = FALSE], sep = " / ")
  }
  records$at_risk <- x$at_risk[pr$keep]

  ggplot2::ggplot(records, ggplot2::aes(x = .data[[score_var]])) +
    ggplot2::geom_histogram(
      ggplot2::aes(fill = at_risk),
      bins = bins,
      color = "white",
      alpha = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = x$threshold,
      color = "red",
      linetype = "dashed",
      linewidth = 1
    ) +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "coral"),
      labels = c("FALSE" = "Not at Risk", "TRUE" = "At Risk"),
      name = "Status"
    ) +
    ggplot2::facet_wrap(~ facet_var, scales = "free_y") +
    ggplot2::coord_cartesian(xlim = pr$xlim) +
    ggplot2::labs(
      title = paste(pr$label, "Distribution"),
      subtitle = sprintf("Overall: %.1f%% at risk (threshold = %.2f)", x$pct_at_risk,
                         x$threshold),
      x = pr$label,
      y = "Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

#' Histogram of the score in base graphics (point rule and the other rules)
#' @noRd
.plot_score_histogram <- function(x, xlim, annotate, bins) {
  pr <- .rapid_plot_prep(x, xlim)
  if (is.null(pr)) return(NULL)
  scores <- pr$scores
  below <- pr$below
  xlim_values <- pr$xlim
  score_range <- pr$score_range
  x_label <- pr$label
  main_title <- paste(x_label, "Distribution")
  risk_labels <- paste0(c("At Risk (", "Not at Risk ("),
                        if (below) c("< ", ">= ") else c("> ", "<= "),
                        format(x$threshold), ")")

  # Create histogram breaks that span the xlim range
  break_min <- min(xlim_values[1], score_range[1])
  break_max <- max(xlim_values[2], score_range[2])
  if (break_max <= break_min) break_max <- break_min + 1
  breaks_seq <- seq(break_min, break_max, length.out = bins + 1)

  h <- hist(scores, breaks = breaks_seq, plot = FALSE)

  # Color bars based on whether they're in the at-risk region
  bar_colors <- if (below) {
    # Low error means an accurate prediction, hence disclosure risk
    ifelse(h$mids < x$threshold, "coral", "steelblue")
  } else {
    ifelse(h$mids >= x$threshold, "coral", "steelblue")
  }

  # Plot histogram with explicit xlim (don't pass ... to avoid xlim conflict)
  plot(h, main = main_title,
       xlab = x_label, col = bar_colors, border = "white",
       xlim = xlim_values, ylab = "Frequency")
  abline(v = x$threshold, col = "red", lwd = 2, lty = 2)

  if (annotate) {
    usr <- par("usr")
    text_x <- usr[2] - 0.05 * (usr[2] - usr[1])
    text_y <- usr[4] - 0.05 * (usr[4] - usr[3])
    text(text_x, text_y,
         labels = sprintf("%.1f%% at risk\n(n = %d / %d)",
                          x$pct_at_risk, x$n_at_risk, x$n_total),
         adj = c(1, 1), cex = 0.9, font = 2)
  }

  legend(if (below) "topright" else "topleft",
         legend = c(risk_labels, "Threshold"),
         fill = c("coral", "steelblue", NA),
         border = c("white", "white", NA),
         col = c(NA, NA, "red"),
         lty = c(NA, NA, 2),
         lwd = c(NA, NA, 2),
         cex = 0.8)
  invisible(NULL)
}

#' Predicted against actual values in base graphics (numeric targets)
#'
#' Draws nothing for a categorical target.
#' @noRd
.plot_prediction_scatter <- function(x, annotate) {
  if (x$is_categorical) return(invisible(NULL))
  # Plot not-at-risk points first (background), then at-risk points on top
  col_at_risk <- adjustcolor("coral", alpha.f = 0.6)
  col_not_risk <- adjustcolor("steelblue", alpha.f = 0.6)
  actual <- .rapid_col(x, "original_value")
  predicted <- .rapid_col(x, "predicted_value")

  plot(range(actual), range(predicted),
       type = "n",
       main = "Predicted vs Actual Values",
       xlab = "Actual Value", ylab = "Predicted Value")

  if (any(!x$at_risk)) {
    points(actual[!x$at_risk], predicted[!x$at_risk],
           col = col_not_risk, pch = 19, cex = 0.7)
  }
  if (any(x$at_risk)) {
    points(actual[x$at_risk], predicted[x$at_risk],
           col = col_at_risk, pch = 19, cex = 0.7)
  }

  abline(0, 1, col = "darkgray", lwd = 2, lty = 2)

  if (annotate) {
    usr <- par("usr")
    text_x <- usr[1] + 0.05 * (usr[2] - usr[1])
    text_y <- usr[4] - 0.05 * (usr[4] - usr[3])
    text(text_x, text_y,
         labels = sprintf("%.1f%% at risk\n(n = %d / %d)",
                          x$pct_at_risk, x$n_at_risk, x$n_total),
         adj = c(0, 1), cex = 0.9, font = 2)
  }

  legend("bottomright",
         legend = c("At Risk", "Not at Risk", "Perfect Prediction"),
         col = c(col_at_risk, col_not_risk, "darkgray"),
         pch = c(19, 19, NA), lty = c(NA, NA, 2), lwd = c(NA, NA, 2),
         cex = 0.8)
  invisible(NULL)
}

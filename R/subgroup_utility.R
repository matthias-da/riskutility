#' Stratified Utility Assessment Across Subgroups
#'
#' Evaluates whether synthetic data utility holds across subgroups defined by
#' a grouping variable. Synthetic data can look good overall but perform poorly
#' for minority subgroups. This function applies a utility measure to each
#' subgroup and identifies groups with low utility.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the synthetic/anonymized dataset.
#' @param group_var Character string naming the column to stratify by. Must be
#'   present in both datasets. If the column is numeric, it is converted to a
#'   factor with a message.
#' @param utility_fun Function to apply per subgroup. Must accept two data.frames
#'   as its first two arguments and return an object with a \code{$utility_score}
#'   element. Default \code{\link{energy_distance}}.
#' @param threshold Numeric, subgroups with utility below this value are flagged.
#'   Default 0.5.
#' @param na.rm Logical, whether to remove rows where \code{group_var} is NA.
#'   Default TRUE.
#' @param ... Additional arguments passed to \code{utility_fun}.
#'
#' @return An object of class \code{"subgroup_utility"} containing:
#' \itemize{
#'   \item \code{overall_score}: utility computed on the full (ungrouped) data
#'   \item \code{utility_score}: worst subgroup score (conservative measure)
#'   \item \code{per_group}: data.frame with columns \code{group}, \code{n_orig},
#'     \code{n_synth}, \code{utility_score}, \code{flagged}
#'   \item \code{worst_group}: name of the worst-performing subgroup
#'   \item \code{ratio}: worst-subgroup score / overall score
#'   \item \code{group_var}: grouping variable name
#'   \item \code{threshold}: threshold used for flagging
#'   \item \code{n_groups}: number of subgroups evaluated (excluding skipped)
#' }
#'
#' @details
#' The function splits both original and synthetic data by the levels of
#' \code{group_var}, then applies \code{utility_fun} to each subgroup pair.
#' Subgroups with fewer than 5 observations in either dataset are skipped
#' (their utility is set to NA) with a warning.
#'
#' The conservative \code{utility_score} returned is the minimum across all
#' evaluated subgroups, reflecting a worst-case perspective. The \code{ratio}
#' of worst to overall utility helps identify whether particular subgroups
#' are disproportionately affected.
#'
#' @seealso \code{\link{energy_distance}}, \code{\link{mmd}},
#'   \code{\link{hellinger}} for utility functions that can be passed as
#'   \code{utility_fun}
#'
#' @author Matthias Templ
#' @family utility
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' X <- data.frame(
#'   group = sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.6, 0.3, 0.1)),
#'   x1 = rnorm(n),
#'   x2 = rnorm(n)
#' )
#' # Good synthetic for groups A and B, poor for C
#' Y <- X
#' Y$x1 <- Y$x1 + ifelse(Y$group == "C", 3, rnorm(n, 0, 0.1))
#' Y$x2 <- Y$x2 + ifelse(Y$group == "C", 3, rnorm(n, 0, 0.1))
#'
#' result <- subgroup_utility(X, Y, group_var = "group", seed = 42)
#' print(result)
#' summary(result)
#'
#' \donttest{
#' plot(result)
#' }
subgroup_utility <- function(X, ...) {
  UseMethod("subgroup_utility")
}

#' @rdname subgroup_utility
#' @export
subgroup_utility.synth_pair <- function(X, ...) {
  subgroup_utility.default(
    X = X$original,
    Y = X$synthetic,
    ...
  )
}

#' @rdname subgroup_utility
#' @export
subgroup_utility.default <- function(X, Y,
                                     group_var,
                                     utility_fun = energy_distance,
                                     threshold = 0.5,
                                     na.rm = TRUE,
                                     ...) {

  # Convert to data.frame if needed
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)

  # Validate group_var
  if (missing(group_var) || !is.character(group_var) || length(group_var) != 1) {
    stop("'group_var' must be a single character string naming a column.")
  }

  if (!group_var %in% names(X)) {
    stop(paste0("'group_var' (\"", group_var, "\") not found in X."))
  }
  if (!group_var %in% names(Y)) {
    stop(paste0("'group_var' (\"", group_var, "\") not found in Y."))
  }

  # Convert numeric group_var to factor
  if (is.numeric(X[[group_var]])) {
    message(paste0("Converting numeric '", group_var, "' to factor."))
    X[[group_var]] <- as.factor(X[[group_var]])
    Y[[group_var]] <- as.factor(Y[[group_var]])
  }

  # Ensure factor/character
  X[[group_var]] <- as.factor(X[[group_var]])
  Y[[group_var]] <- as.factor(Y[[group_var]])

  # Handle NAs in group_var
  if (na.rm) {
    X <- X[!is.na(X[[group_var]]), , drop = FALSE]
    Y <- Y[!is.na(Y[[group_var]]), , drop = FALSE]
  }

  # Compute overall utility (on full data, excluding group_var from utility computation)
  analysis_vars <- setdiff(names(X), group_var)
  analysis_vars <- intersect(analysis_vars, names(Y))
  X_full <- X[, analysis_vars, drop = FALSE]
  Y_full <- Y[, analysis_vars, drop = FALSE]
  overall_result <- utility_fun(X_full, Y_full, ...)
  overall_score <- overall_result$utility_score

  # Get all levels present in either dataset
  all_levels <- sort(unique(c(levels(X[[group_var]]), levels(Y[[group_var]]))))
  # Keep only levels with data in at least one dataset
  all_levels <- all_levels[all_levels %in% c(as.character(X[[group_var]]),
                                              as.character(Y[[group_var]]))]

  # Per-group computation
  group_names <- character(0)
  n_orig_vec <- integer(0)
  n_synth_vec <- integer(0)
  utility_vec <- numeric(0)

  for (lev in all_levels) {
    X_sub <- X[X[[group_var]] == lev, analysis_vars, drop = FALSE]
    Y_sub <- Y[Y[[group_var]] == lev, analysis_vars, drop = FALSE]
    n_x <- nrow(X_sub)
    n_y <- nrow(Y_sub)

    group_names <- c(group_names, lev)
    n_orig_vec <- c(n_orig_vec, n_x)
    n_synth_vec <- c(n_synth_vec, n_y)

    if (n_x < 5 || n_y < 5) {
      warning(paste0("Subgroup '", lev, "' has fewer than 5 observations ",
                     "(n_orig=", n_x, ", n_synth=", n_y, "). Skipping."))
      utility_vec <- c(utility_vec, NA_real_)
    } else {
      sub_result <- utility_fun(X_sub, Y_sub, ...)
      utility_vec <- c(utility_vec, sub_result$utility_score)
    }
  }

  per_group <- data.frame(
    group = group_names,
    n_orig = n_orig_vec,
    n_synth = n_synth_vec,
    utility_score = utility_vec,
    flagged = ifelse(is.na(utility_vec), NA, utility_vec < threshold),
    stringsAsFactors = FALSE
  )

  # Worst subgroup (among non-NA)
  valid_scores <- utility_vec[!is.na(utility_vec)]
  n_groups <- length(valid_scores)
  if (n_groups > 0) {
    worst_idx <- which.min(utility_vec)
    worst_group <- group_names[worst_idx]
    worst_score <- utility_vec[worst_idx]
    ratio <- if (overall_score > 0) worst_score / overall_score else NA_real_
  } else {
    worst_group <- NA_character_
    worst_score <- NA_real_
    ratio <- NA_real_
  }

  result <- list(
    overall_score = overall_score,
    utility_score = worst_score,
    per_group = per_group,
    worst_group = worst_group,
    ratio = ratio,
    group_var = group_var,
    threshold = threshold,
    n_groups = n_groups
  )

  class(result) <- "subgroup_utility"
  return(result)
}


#' Print method for subgroup_utility objects
#'
#' @param x an object of class \code{"subgroup_utility"}
#' @param ... additional arguments (ignored)
#' @export
print.subgroup_utility <- function(x, ...) {
  cat("Subgroup Utility Assessment\n")
  cat("===========================\n\n")

  cat("Group variable:", x$group_var, "\n")
  cat("Number of subgroups:", x$n_groups, "\n")
  cat("Threshold:", x$threshold, "\n\n")

  cat("Overall utility score:", sprintf("%.4f", x$overall_score), "\n")
  cat("Worst subgroup score: ", sprintf("%.4f", x$utility_score),
      " (", x$worst_group, ")\n", sep = "")
  cat("Worst / Overall ratio:", sprintf("%.4f", x$ratio), "\n\n")

  n_flagged <- sum(x$per_group$flagged, na.rm = TRUE)
  n_skipped <- sum(is.na(x$per_group$utility_score))
  if (n_flagged > 0) {
    flagged_groups <- x$per_group$group[which(x$per_group$flagged)]
    cat("Flagged subgroups (utility <", x$threshold, "):",
        paste(flagged_groups, collapse = ", "), "\n")
  } else {
    cat("No subgroups flagged (all above threshold", x$threshold, ").\n")
  }
  if (n_skipped > 0) {
    skipped_groups <- x$per_group$group[is.na(x$per_group$utility_score)]
    cat("Skipped subgroups (< 5 obs):",
        paste(skipped_groups, collapse = ", "), "\n")
  }

  invisible(x)
}


#' Summary method for subgroup_utility objects
#'
#' @param object an object of class \code{"subgroup_utility"}
#' @param ... additional arguments (ignored)
#' @export
summary.subgroup_utility <- function(object, ...) {
  summ <- list(
    overall_score = object$overall_score,
    utility_score = object$utility_score,
    per_group = object$per_group,
    worst_group = object$worst_group,
    ratio = object$ratio,
    group_var = object$group_var,
    threshold = object$threshold,
    n_groups = object$n_groups,
    mean_utility = mean(object$per_group$utility_score, na.rm = TRUE),
    sd_utility = if (sum(!is.na(object$per_group$utility_score)) > 1)
      sd(object$per_group$utility_score, na.rm = TRUE) else NA_real_
  )

  class(summ) <- "summary.subgroup_utility"
  return(summ)
}


#' Print method for summary.subgroup_utility objects
#'
#' @param x an object of class \code{"summary.subgroup_utility"}
#' @param ... additional arguments (ignored)
#' @export
print.summary.subgroup_utility <- function(x, ...) {
  cat("Summary: Subgroup Utility Assessment\n")
  cat("=====================================\n\n")

  cat("Group variable:", x$group_var, "\n")
  cat("Threshold:", x$threshold, "\n\n")

  cat("Aggregate Scores:\n")
  cat("  Overall:  ", sprintf("%.4f", x$overall_score), "\n")
  cat("  Mean:     ", sprintf("%.4f", x$mean_utility), "\n")
  if (!is.na(x$sd_utility)) {
    cat("  SD:       ", sprintf("%.4f", x$sd_utility), "\n")
  }
  cat("  Worst:    ", sprintf("%.4f", x$utility_score),
      " (", x$worst_group, ")\n", sep = "")
  cat("  Ratio:    ", sprintf("%.4f", x$ratio), "\n\n")

  cat("Per-Group Results:\n")
  pg <- x$per_group
  pg$utility_score <- sprintf("%.4f", pg$utility_score)
  pg$flagged <- ifelse(is.na(x$per_group$flagged), "skipped",
                       ifelse(x$per_group$flagged, "YES", "no"))
  print(pg, row.names = FALSE, right = FALSE)

  invisible(x)
}


#' Plot method for subgroup_utility objects
#'
#' @param x an object of class \code{"subgroup_utility"}
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = bar chart of utility scores per subgroup
#' @importFrom graphics barplot abline legend axis par mtext
#' @export
plot.subgroup_utility <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 1)
  show[which] <- TRUE

  if (show[1]) {
    pg <- x$per_group

    # Prepare values; replace NA with 0 for plotting
    vals <- pg$utility_score
    bar_names <- pg$group
    is_na <- is.na(vals)
    vals[is_na] <- 0

    # Color: flagged = coral, OK = steelblue, skipped = grey
    cols <- ifelse(is_na, "grey70",
                   ifelse(!is.na(pg$flagged) & pg$flagged, "coral", "steelblue"))

    bp <- barplot(vals,
                  names.arg = bar_names,
                  main = paste("Subgroup Utility by", x$group_var),
                  ylab = "Utility Score",
                  xlab = x$group_var,
                  col = cols,
                  ylim = c(0, max(1, max(vals, na.rm = TRUE) * 1.1)),
                  las = 1,
                  ...)

    # Threshold line
    abline(h = x$threshold, lty = 2, col = "red", lwd = 1.5)

    # Overall reference line
    abline(h = x$overall_score, lty = 3, col = "darkgreen", lwd = 1.5)

    legend("topright",
           legend = c(
             sprintf("Threshold = %.2f", x$threshold),
             sprintf("Overall = %.2f", x$overall_score)
           ),
           col = c("red", "darkgreen"),
           lty = c(2, 3),
           lwd = c(1.5, 1.5),
           bty = "n",
           cex = 0.8)
  }

  invisible(x)
}

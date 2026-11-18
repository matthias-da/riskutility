#' Cross-Validation of Synthesizer Disclosure Risk
#'
#' Evaluates the RAPID disclosure risk of a synthesis method using
#' k-fold cross-validation.
#' In each fold the synthesizer is trained on the training split and
#' RAPID is evaluated on the held-out test split, giving a
#' distribution of risk scores that accounts for variability in both
#' the synthesis and the evaluation.  Every fold is one call of
#' \code{\link{rapid}}, so the fold values come with their chance level
#' and excess.
#'
#' @param X A data frame of original (confidential) records.
#' @param synthesizer A function that takes a data frame (and optionally
#'   a \code{seed} argument) and returns a synthetic data frame.  If
#'   \code{NULL}, \code{synthpop::syn()} is used as a default.
#' @param key_vars Character vector of quasi-identifier column names.
#' @param target_var Character string naming the sensitive (target)
#'   variable.
#' @param k Number of folds (default 5).
#' @param stratified Logical; use stratified folds for categorical
#'   targets (default \code{TRUE})?
#' @param seed Optional random seed; there is no default.  With a seed, the
#'   folds are drawn from it and the synthesizer of fold \eqn{i} receives
#'   \code{seed + i} when it has a \code{seed} argument (or \code{...}).
#'   A synthesizer without one draws from the random stream the seed started.
#' @param verbose Logical; print progress messages?
#' @param ... Additional arguments passed to \code{\link{rapid}}
#'   (e.g. \code{model_type}, \code{cat_tau}, \code{num_tau},
#'   \code{n_perm}).  The old argument names \code{original_data},
#'   \code{quasi_identifiers} and \code{sensitive_attribute} are accepted
#'   as aliases of \code{X}, \code{key_vars} and \code{target_var}.
#' @param return_details,return_all_records Deprecated and ignored: the
#'   object always holds the per-fold table \code{folds} and the fold
#'   objects \code{fits}.
#'
#' @return An object of class \code{"rapid_cv"} with components:
#' \describe{
#'   \item{folds}{Data frame with one row per fold: \code{fold},
#'     \code{n} (records of the test split), \code{rapid},
#'     \code{chance}, \code{excess}, the model metrics and the
#'     \code{threshold} of the fold (it differs between folds under
#'     \code{threshold_rule = "permutation"}).}
#'   \item{fits}{List of the \code{rapid} objects of the folds.}
#'   \item{cv_summary}{List with \code{mean}, \code{sd}, \code{se},
#'     \code{median}, \code{min}, \code{max}, \code{ci_lower},
#'     \code{ci_upper} of the fold-level RAPID values, their
#'     \code{mean_chance} and \code{mean_excess}, plus the mean model
#'     metrics.}
#'   \item{rapid, chance, excess}{Means over the folds.}
#'   \item{rule}{The scoring rule, taken from the first fold. Under a
#'     permutation threshold the print adds the range of the fold thresholds.}
#'   \item{method, threshold, model_type, attacker}{Method and threshold of
#'     the rule, the attacker model and its settings line, all taken from the
#'     first fold.}
#'   \item{settings}{List of the settings: \code{k}, \code{cv_type},
#'     \code{stratified}, \code{is_categorical}, \code{n_original},
#'     \code{target_var}, \code{key_vars}.}
#' }
#' The interval \code{ci_lower}, \code{ci_upper} is the t interval
#' \eqn{\bar R \pm t_{0.975, k-1} s / \sqrt{k}}{mean +/- qt(0.975, k - 1) sd /
#' sqrt(k)} of the mean of the \eqn{k} fold values, where \eqn{s} is their
#' standard deviation. It is a rough guide: the folds share the original data,
#' so the fold values are not independent.
#' \code{format()} returns the printed lines and \code{print()} the object,
#' invisibly. \code{summary()} returns an object of class
#' \code{"summary.rapid_cv"} with the statistics over the folds and the fold
#' table. \code{plot()} draws RAPID and the chance level by fold with the mean
#' and its interval and returns the ggplot object, invisibly.
#' \code{as.data.frame()} returns \code{folds}.
#'
#' @references
#' Templ, M., Thees, O. and Müller, R. (2026). RAPID: Risk of Attribute
#' Prediction-Induced Disclosure. Manuscript under review at EPJ Data Science.
#'
#' @seealso \code{\link{rapid}}, \code{\link{rapid_test}}
#'
#' @family rapid
#' @author Oscar Thees, Matthias Templ
#' @examples
#' set.seed(1)
#' df <- data.frame(age = sample(20:80, 150, replace = TRUE),
#'                  sex = factor(sample(c("m", "f"), 150, replace = TRUE)))
#' df$income <- 20000 + 500 * df$age + rnorm(150, 0, 5000)
#' noisy <- function(data, seed = NULL) {
#'   if (!is.null(seed)) set.seed(seed)
#'   data$income <- data$income + rnorm(nrow(data), 0, 3000)
#'   data[sample(nrow(data)), ]
#' }
#' cv <- rapid_synthesizer_cv(df, noisy, key_vars = c("age", "sex"),
#'                            target_var = "income", model_type = "lm",
#'                            k = 3, seed = 1)
#' cv
#' summary(cv)
#' as.data.frame(cv)
#' plot(cv)
#' \donttest{
#' # With synthpop as the synthesizer
#' if (requireNamespace("synthpop", quietly = TRUE)) {
#'   cv <- rapid_synthesizer_cv(
#'     df,
#'     synthesizer = function(data, seed = NULL) {
#'       synthpop::syn(data, m = 1, seed = seed, print.flag = FALSE)$syn
#'     },
#'     key_vars = c("age", "sex"), target_var = "income",
#'     model_type = "rf", k = 5, seed = 1
#'   )
#'   print(cv)
#' }
#' }
#'
#' @export
rapid_synthesizer_cv <- function(X,
                                 synthesizer = NULL,
                                 key_vars,
                                 target_var,
                                 k = 5,
                                 stratified = TRUE,
                                 seed = NULL,
                                 verbose = FALSE,
                                 ...,
                                 return_details,
                                 return_all_records) {
  deprecated <- c(return_details = !missing(return_details),
                  return_all_records = !missing(return_all_records))
  for (arg in names(deprecated)[deprecated]) {
    .rapid_warn(c("'", arg, "' is deprecated and ignored: the object always keeps the ",
                  "per-fold table 'folds' and the fold objects 'fits'."), "deprecated")
  }
  resolved <- .rapid_resolve_aliases(
    list(X = if (missing(X)) NULL else X,
         key_vars = if (missing(key_vars)) NULL else key_vars,
         target_var = if (missing(target_var)) NULL else target_var),
    list(...))
  X <- resolved$arg$X
  key_vars <- resolved$arg$key_vars
  target_var <- resolved$arg$target_var
  dots <- resolved$dots
  .rapid_check_reserved(dots, c(
    Y = "the synthetic data of each fold come from 'synthesizer'",
    synthetic_data = "the synthetic data of each fold come from 'synthesizer'"),
    "rapid_synthesizer_cv")
  if (!is.data.frame(X)) {
    .rapid_abort("The original data are missing: pass 'X' (a data frame).", "input")
  }
  if (!is.character(key_vars) || length(key_vars) == 0 || !is.character(target_var) ||
      length(target_var) != 1) {
    .rapid_abort("'key_vars' and 'target_var' must be given as variable names.", "input")
  }
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k < 2 || nrow(X) < k) {
    .rapid_abort("'k' must be a single number of at least 2 and at most the number of rows.",
                 "input")
  }
  if (!is.null(synthesizer) && !is.function(synthesizer)) {
    .rapid_abort("'synthesizer' must be a function.", "input")
  }
  if (!is.null(seed)) set.seed(seed)

  # Default synthesizer: synthpop
  if (is.null(synthesizer)) {
    if (!requireNamespace("synthpop", quietly = TRUE)) {
      .rapid_abort(c("Package 'synthpop' is required for the default synthesizer. ",
                     "Install it or provide a custom synthesizer function."), "dependency")
    }
    synthesizer <- function(data, seed = NULL) {
      synthpop::syn(data, m = 1, seed = seed, print.flag = FALSE)$syn
    }
  }

  # A synthesizer without a seed argument draws from the stream seeded above
  syn_args <- names(formals(args(synthesizer)))
  takes_seed <- any(c("seed", "...") %in% syn_args)

  target_vec <- X[[target_var]]
  is_categorical <- is.factor(target_vec) || is.character(target_vec) || is.logical(target_vec)

  # ---- Create folds ----
  if (stratified && is_categorical) {
    folds <- .create_stratified_folds(target_vec, k)
  } else {
    folds <- split(seq_len(nrow(X)), sample(rep(seq_len(k), length.out = nrow(X))))
  }

  if (verbose) {
    msg <- sprintf("  RAPID Synthesizer CV (%d folds", k)
    if (stratified && is_categorical) msg <- paste0(msg, ", stratified")
    message(msg, ")")
  }

  # ---- Run CV ----
  fits <- vector("list", k)
  for (fold_id in seq_len(k)) {
    if (verbose) message(sprintf("  Fold %d/%d: ", fold_id, k), appendLF = FALSE)

    test_idx  <- folds[[fold_id]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)

    # Synthesize from training split
    if (verbose) message(sprintf("synthesizing %d  ", length(train_idx)), appendLF = FALSE)
    synth <- if (!is.null(seed) && takes_seed) {
      synthesizer(X[train_idx, ], seed = seed + fold_id)
    } else {
      synthesizer(X[train_idx, ])
    }

    # Evaluate RAPID on test split
    if (verbose) message(sprintf("evaluating %d... ", length(test_idx)), appendLF = FALSE)
    fits[[fold_id]] <- do.call(rapid, c(
      list(X = X[test_idx, ], Y = synth, key_vars = key_vars, target_var = target_var,
           verbose = FALSE),
      dots))
    if (verbose) message(sprintf("done (RAPID: %.3f)", fits[[fold_id]]$rapid))
  }

  # ---- Aggregate ----
  fold_df <- do.call(rbind, lapply(seq_len(k), function(i) {
    f <- fits[[i]]
    data.frame(fold = i, n = f$n_total, rapid = f$rapid, chance = f$chance,
               excess = f$excess, as.data.frame(f$model_metrics), threshold = f$threshold)
  }))
  rapids <- fold_df$rapid
  half <- stats::qt(0.975, k - 1) * stats::sd(rapids) / sqrt(k)

  cv_summary <- list(
    mean   = mean(rapids),
    sd     = stats::sd(rapids),
    se     = stats::sd(rapids) / sqrt(k),
    median = stats::median(rapids),
    min    = min(rapids),
    max    = max(rapids),
    ci_lower = mean(rapids) - half,
    ci_upper = mean(rapids) + half,
    mean_chance = mean(fold_df$chance),
    mean_excess = mean(fold_df$excess)
  )
  metrics <- setdiff(names(fold_df), c("fold", "n", "rapid", "chance", "excess", "threshold"))
  for (m in metrics) {
    cv_summary[[paste0("mean_", m)]] <- mean(fold_df[[m]], na.rm = TRUE)
  }
  if (is_categorical) cv_summary$sd_accuracy <- stats::sd(fold_df$accuracy, na.rm = TRUE)

  out <- structure(list(
    folds       = fold_df,
    fits        = fits,
    cv_summary  = cv_summary,
    rapid       = cv_summary$mean,
    chance      = cv_summary$mean_chance,
    excess      = cv_summary$mean_excess,
    rule        = fits[[1]]$rule,
    method      = fits[[1]]$method,
    threshold   = fits[[1]]$threshold,
    model_type  = fits[[1]]$model_type,
    attacker    = fits[[1]]$attacker,
    settings    = list(
      k = k,
      cv_type = "synthesizer",
      stratified = stratified,
      is_categorical = is_categorical,
      n_original = nrow(X),
      target_var = target_var,
      key_vars = key_vars
    )
  ), class = "rapid_cv")

  if (verbose) {
    message(sprintf("  Mean RAPID: %.4f [%.4f, %.4f]",
                    cv_summary$mean, cv_summary$ci_lower, cv_summary$ci_upper))
  }
  out
}


# ---- Internal helper: stratified folds ----

#' Create stratified k-fold indices
#' @noRd
.create_stratified_folds <- function(y, k) {
  folds <- vector("list", k)
  for (cls in unique(y)) {
    class_indices <- which(y == cls)
    shuffled <- sample(class_indices)
    fold_assignment <- rep(seq_len(k), length.out = length(class_indices))
    for (fold_id in seq_len(k)) {
      folds[[fold_id]] <- c(folds[[fold_id]],
                            shuffled[fold_assignment == fold_id])
    }
  }
  lapply(folds, sample)
}


# A permutation threshold is calibrated in each fold: the rule line shows the
# first, this line the range
.rapid_cv_threshold_line <- function(x, f) {
  th <- x$folds$threshold
  if (is.null(th) || length(unique(th)) < 2) return(NULL)
  paste0("  The threshold differs by fold: ", f(min(th)), " to ", f(max(th)),
         " (folds$threshold)")
}


#' @rdname rapid_synthesizer_cv
#' @param x,object an object of class \code{"rapid_cv"}
#' @param digits number of decimal places
#' @export
format.rapid_cv <- function(x, digits = 3, ...) {
  f <- function(v) formatC(v, digits = digits, format = "f")
  s <- x$cv_summary
  c(
    paste0("RAPID cross-validation of a synthesizer (", x$settings$k, " folds",
           if (x$settings$stratified && x$settings$is_categorical) ", stratified", ")"),
    paste0("Target '", x$settings$target_var, "', attacker ", x$attacker),
    paste0("Rule: ", format(x$rule)),
    .rapid_cv_threshold_line(x, f),
    paste0("RAPID: mean ", f(s$mean), ", SD ", f(s$sd), ", 95% CI [", f(s$ci_lower), ", ",
           f(s$ci_upper), "]"),
    if (!is.na(s$mean_chance)) {
      paste0("Chance level: mean ", f(s$mean_chance), "; excess: mean ", f(s$mean_excess))
    },
    "Next: summary(), plot(), as.data.frame()"
  )
}

#' @rdname rapid_synthesizer_cv
#' @export
print.rapid_cv <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname rapid_synthesizer_cv
#' @export
summary.rapid_cv <- function(object, ...) {
  structure(c(
    list(rule = format(object$rule),
         method = object$method,
         threshold = object$threshold,
         model_type = object$model_type,
         attacker = object$attacker,
         k = object$settings$k,
         stratified = object$settings$stratified,
         is_categorical = object$settings$is_categorical,
         n_original = object$settings$n_original,
         target_var = object$settings$target_var,
         key_vars = object$settings$key_vars),
    object$cv_summary[c("mean", "sd", "se", "median", "min", "max", "ci_lower", "ci_upper",
                        "mean_chance", "mean_excess")],
    list(cv_summary = object$cv_summary, folds = object$folds)
  ), class = "summary.rapid_cv")
}

#' @rdname rapid_synthesizer_cv
#' @export
format.summary.rapid_cv <- function(x, ...) {
  f <- function(v) formatC(v, digits = 4, format = "f")
  s <- x$cv_summary
  perf <- if (x$is_categorical) {
    if (!is.null(s$mean_accuracy)) paste0("  Mean accuracy: ", f(s$mean_accuracy))
  } else if (!is.null(s$mean_mae)) {
    c(paste0("  Mean MAE:  ", f(s$mean_mae)), paste0("  Mean RMSE: ", f(s$mean_rmse)))
  }
  c(
    "Summary: RAPID Cross-Validation",
    "================================",
    "Configuration:",
    paste0("  Rule:                ", x$rule),
    paste0("  Attacker:            ", x$attacker),
    paste0("  Target variable:     ", x$target_var),
    paste0("  K-folds:             ", x$k, if (x$stratified && x$is_categorical) " (stratified)"),
    paste0("  Original data size:  ", x$n_original),
    "RAPID over the folds:",
    paste0("  Mean:   ", f(x$mean)),
    paste0("  Median: ", f(x$median)),
    paste0("  SD:     ", f(x$sd)),
    paste0("  SE:     ", f(x$se)),
    paste0("  Range:  [", f(x$min), ", ", f(x$max), "]"),
    paste0("  95% CI: [", f(x$ci_lower), ", ", f(x$ci_upper), "]"),
    if (!is.na(x$mean_chance)) {
      c(paste0("  Mean chance level: ", f(x$mean_chance)),
        paste0("  Mean excess:       ", f(x$mean_excess)))
    },
    if (length(perf) > 0) c("Model performance:", perf),
    "Folds:",
    paste0("  ", utils::capture.output(print(x$folds, row.names = FALSE, digits = 4)))
  )
}

#' @rdname rapid_synthesizer_cv
#' @export
print.summary.rapid_cv <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname rapid_synthesizer_cv
#' @export
plot.rapid_cv <- function(x, ...) {
  s <- x$cv_summary
  df <- rbind(
    data.frame(fold = x$folds$fold, value = x$folds$rapid, series = "RAPID"),
    if (!anyNA(x$folds$chance)) {
      data.frame(fold = x$folds$fold, value = x$folds$chance, series = "Chance level")
    }
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$fold, y = .data$value,
                                        colour = .data$series, shape = .data$series)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = s$ci_lower, ymax = s$ci_upper,
                      fill = "steelblue", alpha = 0.15) +
    ggplot2::geom_hline(yintercept = s$mean, colour = "steelblue", linetype = "dashed") +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_x_continuous(breaks = x$folds$fold) +
    ggplot2::scale_colour_manual(values = c(RAPID = "steelblue", `Chance level` = "grey50")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Fold", y = "Share of records at risk", colour = NULL, shape = NULL,
      title = "RAPID by cross-validation fold",
      subtitle = sprintf("Dashed = mean RAPID %.3f, band = 95%% interval of the mean", s$mean))
  print(p)
  invisible(p)
}

#' @rdname rapid_synthesizer_cv
#' @param row.names,optional not used
#' @export
as.data.frame.rapid_cv <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$folds
}

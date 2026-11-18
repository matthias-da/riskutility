# Methods for objects of class "rapid"

#' Basic methods for rapid objects
#'
#' \code{format()} returns the lines that \code{print()} shows: the target, the
#' quasi-identifiers and the attacker with the settings that drive its result,
#' the scoring rule, the share of records at risk, and where to go next.
#' \code{as.data.frame()} returns the record table and \code{nobs()} the number
#' of scored records.
#'
#' @param x,object an object of class "rapid"
#' @param digits number of decimal places of the RAPID value
#' @param row.names \code{NULL} or a character vector of row names
#' @param optional not used
#' @param at_risk_only logical, if TRUE return only the records at risk
#' @param ... additional arguments (ignored)
#' @return \code{format()}: a character vector. \code{print()}: the input
#'   object, invisibly. \code{as.data.frame()}: a data frame with one row per
#'   record (see the \code{records} component of \code{\link{rapid}}).
#'   \code{nobs()}: the number of scored records.
#' @examples
#' set.seed(1)
#' X <- data.frame(age = sample(20:60, 100, replace = TRUE))
#' X$income <- 20000 + 600 * X$age + rnorm(100, 0, 3000)
#' Y <- X
#' Y$income <- Y$income + rnorm(100, 0, 2000)
#' r <- rapid(X, Y, key_vars = "age", target_var = "income", model_type = "lm")
#' format(r)
#' print(r)
#' nobs(r)
#' head(as.data.frame(r, at_risk_only = TRUE))
#' @name rapid_methods
#' @family rapid
NULL

#' @rdname rapid_methods
#' @export
format.rapid <- function(x, digits = 3, ...) {
  .rapid_require(x, "version")
  qi <- x$key_vars
  shown <- if (length(qi) > 6) c(qi[1:5], "...") else qi
  chance <- if (!is.na(x$chance)) {
    c(paste0("Chance level = ", formatC(x$chance, digits = digits, format = "f"), " (",
             x$n_perm, " permutation", if (x$n_perm != 1) "s", ")"),
      paste0("Excess = ", formatC(x$excess, digits = digits, format = "f")))
  }
  spread <- if (length(x$chance_rates) > 1) stats::sd(x$chance_rates) else NA_real_
  indistinct <- if (!is.na(spread) && x$excess <= 2 * spread) {
    "RAPID is not distinguishable from the chance level; see rapid_test()."
  }
  c(
    paste0("RAPID of '", x$target_var, "' given ", length(qi), " quasi-identifier",
           if (length(qi) != 1) "s", " (", paste(shown, collapse = ", "), "), attacker ",
           x$attacker),
    paste0("Rule: ", format(x$rule)),
    paste0("RAPID = ", formatC(x$rapid, digits = digits, format = "f"), ": ", x$n_at_risk,
           " of ", x$n_total, " records at risk"),
    chance, indistinct,
    "Next: summary(), plot(), confint(), rapid_test()"
  )
}

#' @rdname rapid_methods
#' @export
print.rapid <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname rapid_methods
#' @export
as.data.frame.rapid <- function(x, row.names = NULL, optional = FALSE, ...,
                                at_risk_only = FALSE) {
  .rapid_require(x, "version")
  records <- x$records
  if (isTRUE(at_risk_only)) records <- records[x$at_risk, , drop = FALSE]
  if (!is.null(row.names)) rownames(records) <- row.names
  records
}

#' @rdname rapid_methods
#' @importFrom stats nobs
#' @export
nobs.rapid <- function(object, ...) {
  .rapid_require(object, "version")
  object$n_total
}

#' Summary method for rapid objects
#'
#' @param object an object of class "rapid"
#' @param ... additional arguments (ignored)
#' @return A list of summary statistics for the corresponding object.
#' @export
summary.rapid <- function(object, ...) {
  .rapid_require(object, "version")
  summ <- list(
    rapid = object$rapid,
    n_at_risk = object$n_at_risk,
    pct_at_risk = object$pct_at_risk,
    n_total = object$n_total,
    method = object$method,
    rule = format(object$rule),
    threshold = object$threshold,
    model_type = object$model_type,
    attacker = object$attacker,
    model_metrics = object$model_metrics,
    is_categorical = object$is_categorical,
    key_vars = object$key_vars,
    target_var = object$target_var,
    risk_level = if (object$rapid < 0.05) "LOW" else if (object$rapid < 0.15) "MEDIUM" else "HIGH"
  )
  class(summ) <- "summary.rapid"
  summ
}

#' Print method for summary.rapid objects
#'
#' @param x an object of class "summary.rapid"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.rapid <- function(x, ...) {
  cat("Summary: RAPID Disclosure Risk Assessment\n")
  cat("==========================================\n\n")

  cat("Configuration:\n")
  cat("  Attacker:", x$attacker, "\n")
  cat("  Rule:", x$rule, "\n")
  cat("  Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("  Target variable:", x$target_var, "\n\n")

  cat("Risk Assessment:\n")
  cat("  RAPID score:", sprintf("%.4f", x$rapid), "\n")
  cat("  Records at risk:", x$n_at_risk, "/", x$n_total,
      sprintf("(%.2f%%)", x$pct_at_risk), "\n")
  cat("  Risk level:", x$risk_level, "\n\n")

  cat("Model Performance:\n")
  if (x$is_categorical) {
    cat("  Accuracy:", sprintf("%.4f (%.1f%%)", x$model_metrics$accuracy,
                              100 * x$model_metrics$accuracy), "\n")
  } else {
    cat("  MAE:", sprintf("%.4f", x$model_metrics$mae), "\n")
    cat("  RMSE:", sprintf("%.4f", x$model_metrics$rmse), "\n")
    cat("  Relative MAE:", sprintf("%.4f", x$model_metrics$rmae), "\n")
    cat("  Relative RMSE:", sprintf("%.4f", x$model_metrics$rrmse), "\n")
  }

  cat("\nInterpretation Guidelines:\n")
  if (x$risk_level == "LOW") {
    cat("  Low risk: Synthetic data provides good privacy protection.\n")
  } else if (x$risk_level == "MEDIUM") {
    cat("  Medium risk: Some records may be at risk. Review the at-risk records.\n")
  } else {
    cat("  High risk: Significant disclosure risk detected. Consider additional protection.\n")
  }

  invisible(x)
}

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

#' Update a rapid object
#'
#' Changes arguments of the call that created a \code{rapid} object. When all
#' changed arguments only affect how the attack is scored (\code{num_eval_method},
#' \code{num_error_metric}, \code{num_epsilon}, \code{num_delta}, \code{num_tau},
#' \code{cat_eval_method}, \code{cat_tau}, \code{n_perm}, \code{threshold_rule},
#' \code{perm_level}, \code{seed}) and the object holds the attacker's
#' predictions, the stored attack is scored again under the new rule and the
#' model is not fitted again. Any other argument, among them the model
#' arguments, \code{num_scale}, \code{num_dist} and \code{prune}, refits: the
#' changed call is evaluated in the caller's frame, as
#' \code{\link[stats]{update}} does.
#'
#' A rescored object equals the result of the changed call, except where the
#' fit depends on a changed argument: under \code{num_scale = "auto"} the
#' working scale of the linear model depends on the error metric, which a
#' rescored object keeps from the first fit.
#'
#' @param object an object of class "rapid"
#' @param ... arguments of \code{\link{rapid}} to change, as \code{name = value}
#' @param evaluate logical, if FALSE return the changed call instead of evaluating it
#' @return A \code{rapid} object whose \code{call} holds the changed arguments,
#'   or the changed call if \code{evaluate = FALSE}.
#' @examples
#' set.seed(1)
#' X <- data.frame(age = sample(20:60, 100, replace = TRUE))
#' X$income <- 20000 + 600 * X$age + rnorm(100, 0, 3000)
#' Y <- X
#' Y$income <- Y$income + rnorm(100, 0, 2000)
#' r <- rapid(X, Y, key_vars = "age", target_var = "income", model_type = "lm")
#' update(r, num_tau = 0.5)                  # rescored, no refit
#' update(r, num_eval_method = "point")      # rescored, no refit
#' update(r, num_tau = 0.5, evaluate = FALSE)
#' @importFrom stats update
#' @family rapid
#' @export
update.rapid <- function(object, ..., evaluate = TRUE) {
  .rapid_require(object, "version")
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras) > 0 && (is.null(names(extras)) || any(!nzchar(names(extras))))) {
    .rapid_abort("Every argument of update() must be named.", "input")
  }

  # Change the stored call (as stats::update.default does)
  call <- object$call
  if (length(extras) > 0) {
    existing <- names(extras) %in% names(call)
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) call <- as.call(c(as.list(call), extras[!existing]))
  }
  if (!isTRUE(evaluate)) return(call)

  env <- parent.frame()
  if (!all(names(extras) %in% .rapid_scoring_args) || is.null(object$attack)) {
    return(eval(call, env))
  }

  # Rescore: the scoring arguments of the changed call, the defaults of rapid() otherwise
  defaults <- formals(rapid.default)
  value <- function(name) {
    if (name %in% names(call)) eval(call[[name]], env)
    else if (name %in% names(defaults)) eval(defaults[[name]])
  }
  scoring <- .rapid_scoring(
    object$is_categorical, names(call)[-1L],
    num_eval_method = value("num_eval_method"), num_error_metric = value("num_error_metric"),
    num_epsilon = value("num_epsilon"),
    num_epsilon_type = if ("num_epsilon_type" %in% names(call)) value("num_epsilon_type"),
    num_delta = value("num_delta"), num_tau = value("num_tau"),
    cat_eval_method = value("cat_eval_method"), cat_tau = value("cat_tau"),
    n_perm = value("n_perm"), threshold_rule = value("threshold_rule"),
    perm_level = value("perm_level")
  )
  seed <- value("seed")
  if (!is.null(seed)) set.seed(seed)
  .rapid_rescore(object, scoring$rule, scoring$n_perm, scoring$threshold_rule,
                 scoring$perm_level, call)
}

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
  f <- function(v) formatC(v, digits = digits, format = "f")
  holdout <- if (!is.null(x$holdout)) {
    h <- x$holdout
    c(paste0("Holdout: RAPID = ", f(h$rapid), " on ", h$n, " record", if (h$n != 1) "s",
             if (!is.na(h$chance)) {
               paste0(", chance level = ", f(h$chance), ", excess = ", f(h$excess))
             }),
      paste0("Release-specific = ", f(x$release_specific), " (RAPID - holdout RAPID)"))
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
    chance, indistinct, holdout, x$note,
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

#' Predictions of the attacker for new records
#'
#' Applies the attacker of a \code{rapid} object to new records. The records
#' are prepared like the original data (column types, \code{na_strategy},
#' factor levels). When they hold the target, they are scored under the rule of
#' the object as calibrated, with baselines taken from the original values.
#'
#' @param object an object of class "rapid", created with \code{store_model = TRUE}
#' @param newdata data frame with the quasi-identifiers and, optionally, the target
#' @param ... not used
#' @return The prepared records. With the target, followed by the record
#'   columns of the rule and \code{at_risk} (as in \code{as.data.frame()});
#'   without it, followed by \code{predicted_class} or \code{predicted_value}.
#' @examples
#' set.seed(1)
#' X <- data.frame(age = sample(20:60, 100, replace = TRUE))
#' X$income <- 20000 + 600 * X$age + rnorm(100, 0, 3000)
#' Y <- X
#' Y$income <- Y$income + rnorm(100, 0, 2000)
#' r <- rapid(X, Y, key_vars = "age", target_var = "income", model_type = "lm",
#'            store_model = TRUE)
#' new <- data.frame(age = c(25, 50), income = c(36000, 50000))
#' predict(r, new)
#' predict(r, new["age"])
#' @family rapid
#' @export
predict.rapid <- function(object, newdata, ...) {
  .rapid_require(object, c("version", "model"))
  .rapid_check_empty_dots(list(...), "predict() takes 'newdata' only.")
  if (missing(newdata)) .rapid_abort("'newdata' is missing: pass the records to predict.",
                                     "input")
  p <- object$predictor
  target <- object$target_var
  new <- .rapid_prepare_new(newdata, p$train, object$key_vars, target, object$na,
                            p$learner, "newdata", need_target = FALSE)
  pred <- .rapid_predict(p$learner, object$model, new$data, p$train, target,
                         object$is_categorical, p$num_dist, p$n_quantiles,
                         point_only = !new$has_target)
  if (new$has_target) {
    ref_y <- .rapid_data(object)[[target]]
    return(.rapid_score_records(object$rule, pred, new$data[[target]], ref_y, new$data,
                                object$is_categorical)$rec$records)
  }
  out <- new$data
  name <- if (object$is_categorical) "predicted_class" else "predicted_value"
  if (name %in% names(out)) name <- paste0(name, ".rapid")
  out[[name]] <- if (object$is_categorical) .rapid_predicted_class(pred) else pred
  out
}

#' Summary method for rapid objects
#'
#' \code{summary()} reports RAPID with its Wilson confidence interval, the
#' chance level and the excess, and how the flagged records are spread over
#' common and rare values of the target. \code{as.data.frame()} returns the
#' one-row table, so that the summaries of several releases can be stacked with
#' \code{rbind()}.
#'
#' A value is rare when the baseline of its record (the share of original
#' records that fall in its tolerance band, or whose class it shares) is below
#' the \code{1 - common_share} quantile of the baselines. Close predictions of
#' common values carry little information: an attacker who always guesses the
#' typical value is often right. \code{rare_hit_rate} is the share of flagged
#' records among the rare ones and \code{rare_rate} the share of all records
#' that are flagged and rare. Both are \code{NA} for a rule without baseline.
#'
#' @param object an object of class "rapid"
#' @param common_share share of records treated as common (between 0 and 1)
#' @param level confidence level of the Wilson interval, as in
#'   \code{\link{confint.rapid}}
#' @param x an object of class "summary.rapid"
#' @param digits number of decimal places
#' @param row.names,optional not used
#' @param ... additional arguments (ignored)
#' @return \code{summary()}: an object of class "summary.rapid", a list with the
#'   one-row data frame \code{table} (\code{target}, \code{model}, \code{rule},
#'   \code{threshold}, \code{n}, \code{rapid}, \code{lower}, \code{upper},
#'   \code{chance}, \code{excess}), \code{rare_share}, \code{n_rare},
#'   \code{rare_hit_rate}, \code{rare_rate} and the settings of the fit.
#'   \code{print()}: the input, invisibly. \code{as.data.frame()}: \code{table}.
#' @examples
#' set.seed(1)
#' X <- data.frame(age = sample(20:60, 100, replace = TRUE))
#' X$income <- 20000 + 600 * X$age + rnorm(100, 0, 3000)
#' Y <- X
#' Y$income <- Y$income + rnorm(100, 0, 2000)
#' r <- rapid(X, Y, key_vars = "age", target_var = "income", model_type = "lm")
#' summary(r)
#' # Summaries of several releases in one table
#' r2 <- rapid(X, Y, key_vars = "age", target_var = "income", model_type = "lm",
#'             num_tau = 0.5)
#' rbind(as.data.frame(summary(r)), as.data.frame(summary(r2)))
#' @family rapid
#' @export
summary.rapid <- function(object, common_share = 0.7, level = 0.95, ...) {
  .rapid_require(object, "version")
  if (!is.numeric(common_share) || length(common_share) != 1 || is.na(common_share) ||
      common_share <= 0 || common_share >= 1) {
    .rapid_abort("'common_share' must be a single number between 0 and 1.", "input")
  }
  ci <- stats::confint(object, level = level)

  # Rare values: baseline below the 1 - common_share quantile
  b <- object$baseline
  if (is.null(b)) {
    n_rare <- NA_integer_
    rare_hit_rate <- rare_rate <- NA_real_
  } else {
    rare <- b < stats::quantile(b, 1 - common_share, names = FALSE)
    n_rare <- sum(rare)
    rare_hit_rate <- if (n_rare > 0) mean(object$at_risk[rare]) else NA_real_
    rare_rate <- mean(object$at_risk & rare)
  }

  structure(list(
    table = data.frame(
      target = object$target_var, model = object$model_type,
      rule = format(object$rule), threshold = object$threshold, n = object$n_total,
      rapid = object$rapid, lower = ci[[1]], upper = ci[[2]],
      chance = object$chance, excess = object$excess, stringsAsFactors = FALSE
    ),
    rare_share = 1 - common_share,
    n_rare = n_rare,
    rare_hit_rate = rare_hit_rate,
    rare_rate = rare_rate,
    level = level,
    n_at_risk = object$n_at_risk,
    pct_at_risk = object$pct_at_risk,
    n_perm = object$n_perm,
    method = object$method,
    attacker = object$attacker,
    model_metrics = object$model_metrics,
    is_categorical = object$is_categorical,
    key_vars = object$key_vars
  ), class = "summary.rapid")
}

#' @rdname summary.rapid
#' @export
format.summary.rapid <- function(x, digits = 4, ...) {
  f <- function(v) formatC(v, digits = digits, format = "f")
  tab <- x$table
  m <- x$model_metrics
  perf <- if (x$is_categorical) {
    paste0("  Accuracy: ", f(m$accuracy))
  } else {
    c(paste0("  MAE: ", f(m$mae)), paste0("  RMSE: ", f(m$rmse)),
      paste0("  Relative MAE: ", f(m$rmae)), paste0("  Relative RMSE: ", f(m$rrmse)))
  }
  c(
    "Summary: RAPID Disclosure Risk Assessment",
    "==========================================",
    "Configuration:",
    paste0("  Attacker: ", x$attacker),
    paste0("  Rule: ", tab$rule),
    paste0("  Key variables: ", paste(x$key_vars, collapse = ", ")),
    paste0("  Target variable: ", tab$target),
    "Risk assessment:",
    paste0("  RAPID: ", f(tab$rapid), " [", f(tab$lower), ", ", f(tab$upper), "] (",
           format(100 * x$level), "% Wilson interval)"),
    paste0("  Records at risk: ", x$n_at_risk, " / ", tab$n, sprintf(" (%.2f%%)", x$pct_at_risk)),
    if (!is.na(tab$chance)) {
      paste0("  Chance level: ", f(tab$chance), " (", x$n_perm, " permutations), excess: ",
             f(tab$excess))
    },
    if (!is.na(x$n_rare)) {
      c(paste0("Rare values (baseline below the ", format(100 * x$rare_share),
               "% quantile; ", x$n_rare, " records):"),
        paste0("  Flagged share among rare records: ", f(x$rare_hit_rate)),
        paste0("  Share of all records flagged and rare: ", f(x$rare_rate)))
    },
    "Model performance:", perf
  )
}

#' @rdname summary.rapid
#' @export
print.summary.rapid <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname summary.rapid
#' @export
as.data.frame.summary.rapid <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$table
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
#' A changed \code{num_error_metric} refits when the model has a working scale
#' chosen under \code{num_scale = "auto"} (the linear model), because that
#' choice depends on the error metric.
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
  # The arguments of the changed call, the defaults of rapid() otherwise
  defaults <- formals(rapid.default)
  value <- function(name) {
    if (name %in% names(call)) eval(call[[name]], env)
    else if (name %in% names(defaults)) eval(defaults[[name]])
  }
  # Under num_scale = "auto" a learner with a working scale chose it by the error metric
  metric_drives_fit <- "num_error_metric" %in% names(extras) && !object$is_categorical &&
    !is.null(object$attack) && "num_scale" %in% object$attack$learner$options &&
    identical(.rapid_choice(value("num_scale"), c("auto", "identity", "log"), "num_scale"),
              "auto")
  if (!all(names(extras) %in% .rapid_scoring_args) || is.null(object$attack) ||
      metric_drives_fit) {
    return(eval(call, env))
  }

  # Rescore the stored attack
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
                 scoring$perm_level, call, seed = seed)
}

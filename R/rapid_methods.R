# Methods for objects of class "rapid"

#' Methods for rapid objects
#'
#' The methods that work on the result of \code{\link{rapid}()}: printing, the
#' summary, the record table, the number of records, changing arguments,
#' applying the attacker to new records and confidence intervals. The plot
#' method has its own page, \code{\link{plot.rapid}}, and the functions that
#' test, calibrate and combine fitted objects are \code{\link{rapid_test}()},
#' \code{\link{rapid_threshold_select}()} and \code{\link{rapid_suite}()}.
#'
#' @param x for \code{format()}, \code{print()} and \code{as.data.frame()}: an
#'   object of class \code{"rapid"} or \code{"summary.rapid"}
#' @param object an object of class \code{"rapid"}
#' @param digits number of decimal places (default 3 for a \code{rapid}
#'   object, 4 for its summary)
#' @param row.names \code{NULL} or a character vector of row names for the
#'   record table (ignored for a summary)
#' @param optional not used
#' @param at_risk_only logical, if \code{TRUE} return only the records at risk
#' @param common_share share of records treated as common, between 0 and 1
#'   (default 0.7). See the section "Summary".
#' @param level confidence level (default 0.95)
#' @param newdata data frame with the quasi-identifiers and, optionally, the
#'   target
#' @param evaluate logical, if \code{FALSE} return the changed call instead of
#'   evaluating it
#' @param parm not used, for compatibility with the generic
#' @param method interval of \code{confint()}: \code{"wilson"} (default),
#'   \code{"wald"}, \code{"clopper-pearson"} or \code{"bootstrap"}
#' @param n_bootstrap number of bootstrap replicates (default 1000), only for
#'   \code{method = "bootstrap"}
#' @param seed random seed of the bootstrap (default \code{NULL})
#' @param ... for \code{update()}: the arguments of \code{\link{rapid}()} to
#'   change, as \code{name = value}. \code{predict()} takes none and stops with
#'   an error of class \code{"rapid_error_input"} when given one. Ignored by the
#'   other methods.
#'
#' @section Printing:
#' \code{format()} returns the lines that \code{print()} shows: the target, the
#' quasi-identifiers and the attacker with the settings that drive its result
#' (field \code{attacker}), the rule, RAPID with the number of records at risk,
#' the chance level and the excess, the result on the holdout and the
#' release-specific RAPID when there is a holdout, and the note of the lookup
#' attacker. A line says that RAPID is not distinguishable from the chance
#' level when the excess is at most twice the standard deviation of the shares
#' flagged in the permutations.
#'
#' @section Summary:
#' \code{summary()} reports RAPID with its Wilson interval at \code{level}, the
#' chance level and the excess, how the flagged records are spread over rare
#' and common values of the target, and the accuracy of the attacker's point
#' prediction or predicted class. \code{as.data.frame()} of the summary returns
#' its one-row table, so that the summaries of several releases or rules can be
#' stacked with \code{rbind()}.
#'
#' A record counts as rare when its baseline \eqn{b_i} is below the \code{1 -
#' common_share} quantile of the baselines of all records. The baseline is that
#' of the rule. Under the band and the point rule it is the share of original
#' values inside the record's tolerance band, and under \code{RCS_marginal} the
#' share of the record's class in the original data. A low baseline then marks a
#' rare value, where a close guess carries information. Under
#' \code{RCS_conditional} the baseline is the mean predicted probability of the
#' true class among the records of that class. A low baseline there marks a
#' class to which the attacker gives little probability on average, which need
#' not be a rare class. \code{rare_hit_rate} is the share of flagged records
#' among the rare ones and \code{rare_rate} the share of all records that are
#' both flagged and rare. Under \code{NCE}, which has no baseline, both are
#' \code{NA}.
#'
#' @section Record table:
#' \code{as.data.frame()} returns the field \code{records}: one row per scored
#' record with the prepared target and quasi-identifiers, the true value of a
#' continuous target (\code{original_value}), the columns of the rule (for the
#' band rule \code{predicted_value, band_lo, band_hi, g, b, r}) and
#' \code{at_risk}. The columns of each rule are listed under the value of
#' \code{\link{rapid}}. \code{nobs()} returns the number of scored records, which
#' is smaller than \code{nrow(X)} when \code{na_strategy = "drop"} removed
#' records.
#'
#' @section Update:
#' \code{update()} changes arguments of the call that created the object. When
#' every changed argument only affects how the attack is scored
#' (\code{num_eval_method}, \code{num_error_metric}, \code{num_epsilon},
#' \code{num_delta}, \code{num_tau}, \code{cat_eval_method}, \code{cat_tau},
#' \code{n_perm}, \code{threshold_rule}, \code{perm_level}, \code{seed}) and the
#' object holds the attacker's predictions, the stored attack is scored again
#' under the new rule and the model is not fitted again. The permutations then
#' start from the seed, as in a fresh call with the same seed. Any other
#' argument, among them the model arguments, \code{model_type},
#' \code{num_scale}, \code{num_dist} and \code{prune}, refits: the changed call
#' is evaluated in the caller's frame, as \code{\link[stats]{update}()} does, so
#' the data must still be there under their names. A changed
#' \code{num_error_metric} also refits when the linear model chose its scale
#' under \code{num_scale = "auto"}, because that choice depends on the metric.
#'
#' @section Predict:
#' \code{predict()} applies the attacker to new records. It needs an object
#' created with \code{store_model = TRUE}, which keeps the model and, in the
#' field \code{predictor}, the prepared synthetic data. The records are prepared
#' like the original data: the column types of the synthetic data, the factor
#' levels, and \code{na_strategy} for the target. A missing target value gets
#' the value a missing value of the original data got (\code{na_constant_value},
#' or the median of the original values), or the record is dropped under
#' \code{na_strategy = "drop"}. Factor levels that occur neither in the original
#' nor in the synthetic data are refused with an error of class
#' \code{"rapid_error_levels"}, and so are, for the linear and logistic model
#' and user-defined attackers, levels of the quasi-identifiers that the
#' synthetic data do not contain. When \code{newdata} holds the target, the
#' records are scored under the rule of the object as calibrated (a permutation
#' threshold stays as found), with baselines from the original values.
#'
#' @section Confidence intervals:
#' RAPID is a share of \eqn{n} flags, \eqn{\hat p = k / n}{p = k / n}. Given the
#' fitted attacker, \code{confint()} treats the flags as independent Bernoulli
#' trials. The interval is for RAPID itself. It does not include the
#' uncertainty of the chance level, for which see \code{\link{rapid_test}()}.
#' \describe{
#'   \item{wilson}{Score interval (Wilson, 1927), the default. It stays inside
#'     \eqn{[0, 1]} and keeps its coverage near 0 and 1.}
#'   \item{wald}{\eqn{\hat p \pm z_{1-\alpha/2} \sqrt{\hat p (1 - \hat p) /
#'     n}}{p +/- z * sqrt(p (1 - p) / n)}, which can leave \eqn{[0, 1]} and has
#'     width 0 when \eqn{\hat p}{p} is 0 or 1.}
#'   \item{clopper-pearson}{Exact interval from beta quantiles (Clopper and
#'     Pearson, 1934), equal to that of \code{binom.test()}, with the bounds 0
#'     at \eqn{k = 0} and 1 at \eqn{k = n}. It is conservative.}
#'   \item{bootstrap}{Percentile interval of \code{n_bootstrap} binomial
#'     counts with share \eqn{\hat p}{p}, which is the same as resampling the
#'     \eqn{n} flags with replacement. The model is not fitted again.}
#' }
#'
#' @return
#' \describe{
#'   \item{\code{format()}}{A character vector of lines.}
#'   \item{\code{print()}}{The input object, invisibly.}
#'   \item{\code{summary()}}{An object of class \code{"summary.rapid"}, a list
#'     with the one-row data frame \code{table} (\code{target}, \code{model},
#'     \code{rule}, \code{threshold}, \code{n}, \code{rapid}, \code{lower},
#'     \code{upper}, \code{chance}, \code{excess}), \code{rare_share},
#'     \code{n_rare}, \code{rare_hit_rate}, \code{rare_rate}, \code{level}, and
#'     \code{n_at_risk}, \code{pct_at_risk}, \code{n_perm}, \code{method},
#'     \code{attacker}, \code{model_metrics}, \code{is_categorical},
#'     \code{key_vars} of the object.}
#'   \item{\code{as.data.frame()}}{For a \code{rapid} object the record table,
#'     for a summary its one-row \code{table}.}
#'   \item{\code{nobs()}}{The number of scored records.}
#'   \item{\code{update()}}{A \code{rapid} object whose \code{call} holds the
#'     changed arguments, or the changed call if \code{evaluate = FALSE}.}
#'   \item{\code{predict()}}{The prepared records. With the target they are
#'     followed by the columns of the rule and \code{at_risk}, as in
#'     \code{as.data.frame()}. Without it they are followed by
#'     \code{predicted_value} (for the linear model on the log scale the
#'     predictive median) or \code{predicted_class}.}
#'   \item{\code{confint()}}{A 1 x 2 matrix with the lower and upper limit,
#'     following the convention of \code{\link[stats]{confint}()}.}
#' }
#'
#' @references
#' Wilson, E.B. (1927). Probable Inference, the Law of Succession, and
#' Statistical Inference. \emph{Journal of the American Statistical Association},
#' 22(158), 209--212.
#'
#' Clopper, C.J. and Pearson, E.S. (1934). The use of confidence or fiducial limits
#' illustrated in the case of the binomial. \emph{Biometrika}, 26(4), 404--413.
#'
#' @seealso \code{\link{rapid}}, \code{\link{plot.rapid}},
#'   \code{\link{rapid_test}}, \code{\link{rapid_threshold_select}},
#'   \code{\link{rapid_suite}}
#' @author Oscar Thees, Matthias Templ
#' @examples
#' set.seed(1)
#' n <- 200
#' X <- data.frame(age = sample(20:70, n, replace = TRUE),
#'                 sex = factor(sample(c("F", "M"), n, replace = TRUE)))
#' X$income <- round(exp(9 + 0.03 * X$age + 0.3 * (X$sex == "M") + rnorm(n, 0, 0.1)))
#' Y <- X
#' Y$income <- round(Y$income * exp(rnorm(n, 0, 0.05)))
#' r <- rapid(X, Y, key_vars = c("age", "sex"), target_var = "income",
#'            model_type = "lm", store_model = TRUE, seed = 1)
#'
#' # Printing and the record table
#' format(r)
#' r
#' nobs(r)
#' head(as.data.frame(r))
#' head(as.data.frame(r, at_risk_only = TRUE))
#'
#' # Summary, and the summaries of two rules in one table
#' s <- summary(r)
#' s
#' r_point <- update(r, num_eval_method = "point")   # rescored, no refit
#' rbind(as.data.frame(s), as.data.frame(summary(r_point)))
#'
#' # update(): the changed call, a rescore, a refit
#' update(r, num_tau = 0.5, evaluate = FALSE)
#' update(r, num_tau = 0.5)
#' update(r, num_scale = "identity")                   # refits the linear model
#'
#' # predict(): new records with and without the target
#' new <- data.frame(age = c(25, 45, 65), sex = c("F", "M", "M"),
#'                   income = c(20000, 40000, 60000))
#' predict(r, new)
#' predict(r, new[c("age", "sex")])
#'
#' # Confidence intervals
#' confint(r)
#' confint(r, method = "clopper-pearson")
#' confint(r, method = "bootstrap", n_bootstrap = 500, seed = 1)
#'
#' # Plots, see ?plot.rapid
#' plot(r, which = 1:2)
#' @name rapid-methods
#' @aliases rapid_methods
#' @family rapid
NULL

#' @rdname rapid-methods
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

#' @rdname rapid-methods
#' @export
print.rapid <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname rapid-methods
#' @export
as.data.frame.rapid <- function(x, row.names = NULL, optional = FALSE, ...,
                                at_risk_only = FALSE) {
  .rapid_require(x, "version")
  records <- x$records
  if (isTRUE(at_risk_only)) records <- records[x$at_risk, , drop = FALSE]
  if (!is.null(row.names)) rownames(records) <- row.names
  records
}

#' @rdname rapid-methods
#' @importFrom stats nobs
#' @export
nobs.rapid <- function(object, ...) {
  .rapid_require(object, "version")
  object$n_total
}

#' @rdname rapid-methods
#' @export
predict.rapid <- function(object, newdata, ...) {
  .rapid_require(object, c("version", "model", "predictor"))
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

#' @rdname rapid-methods
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

#' @rdname rapid-methods
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

#' @rdname rapid-methods
#' @export
print.summary.rapid <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname rapid-methods
#' @export
as.data.frame.summary.rapid <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$table
}

#' @rdname rapid-methods
#' @importFrom stats update
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

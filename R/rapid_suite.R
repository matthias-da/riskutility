#' RAPID over a Suite of Attacker Models
#'
#' Fits one \code{\link{rapid}} object per attacker model on the same data and
#' collects their RAPID, chance level and excess. The RAPID of a suite is
#' judged by its worst case, the strongest attacker, not by one model.
#'
#' Every model is fitted with the arguments in \code{...} and the same
#' \code{seed}, so that \code{suite[["rf"]]} equals
#' \code{rapid(X, Y, key_vars, target_var, model_type = "rf", ..., seed = seed)}.
#' A model whose backend package is not installed, or that does not support the
#' type of the target, is skipped with a warning of class
#' \code{"rapid_warning_skipped_model"}; if no model is left the call stops with
#' an error of class \code{"rapid_error_model"}. Every other error of
#' \code{rapid()} is passed on.
#'
#' The lookup attacker predicts a record from the synthetic records with the
#' same quasi-identifiers. Its RAPID is readable only against its chance level,
#' so compare the models by \code{excess_max} and \code{worst_model_excess},
#' not by \code{rapid_max}, when one is in the suite; the print says so.
#'
#' \code{rapid_suite()} on \code{rapid} objects combines objects that were
#' fitted already. They must agree in \code{target_var}, \code{key_vars},
#' \code{n_total} and scoring rule (else an error of class
#' \code{"rapid_error_object"} names the field). The results are named by the
#' names given in \code{...}, or by the \code{model_type} of the object; two
#' objects with the same name are an error.
#'
#' @param X the original data, a \code{\link{synth_pair}}, or a \code{rapid}
#'   object
#' @param Y data frame of synthetic data (not needed if X is a synth_pair)
#' @param key_vars character vector of quasi-identifier variable names
#' @param target_var character, name of the sensitive target variable
#' @param model_types character vector of attacker models, any of
#'   \code{"rf"}, \code{"lm"}, \code{"cart"}, \code{"gbm"}, \code{"logit"},
#'   \code{"lookup"} (default \code{c("rf", "cart", "gbm", "lookup")})
#' @param ... for the data methods, arguments of \code{\link{rapid}} that go to
#'   every model (for example \code{num_tau}, \code{n_perm}); the old argument
#'   names are accepted. For the \code{rapid} method, further \code{rapid}
#'   objects.
#' @param seed random seed, passed to every \code{rapid()} call (default NULL:
#'   the random stream is left alone)
#'
#' @return An object of class \code{"rapid_suite"} with components:
#' \describe{
#'   \item{results}{Data frame with \code{model_type, rapid, chance, excess,
#'     n_at_risk}, one row per model.}
#'   \item{fits}{Named list of the \code{rapid} objects.}
#'   \item{rapid_mean, rapid_max}{Mean and maximum RAPID over the models.}
#'   \item{worst_model}{Name of the model with the largest RAPID (the first one
#'     in case of ties).}
#'   \item{excess_max, worst_model_excess}{Largest excess over the chance level
#'     and the model that has it (\code{NA} without a chance level, that is
#'     with \code{n_perm = 0}).}
#'   \item{n_flagged_any}{Number of records flagged by at least one model.}
#'   \item{n_total, key_vars, target_var, rule, notes}{Settings shared by the
#'     fits, and the notes of the models that carry one.}
#' }
#' \code{suite[["rf"]]} returns the \code{rapid} object of a model.
#'
#' @seealso \code{\link{rapid}}, \code{\link{rapid_test}}
#'
#' @family rapid
#' @author Oscar Thees, Matthias Templ
#' @examples
#' set.seed(42)
#' X <- data.frame(
#'   age = sample(20:60, 80, replace = TRUE),
#'   sex = sample(c("M", "F"), 80, replace = TRUE),
#'   income = rnorm(80, 50000, 10000)
#' )
#' X$income <- X$income + 500 * X$age
#' Y <- X
#' Y$income <- Y$income + rnorm(80, 0, 5000)
#'
#' # The lm attacker only needs base R
#' s <- rapid_suite(X, Y, key_vars = c("age", "sex"), target_var = "income",
#'                  model_types = "lm", n_perm = 3, seed = 1)
#' s
#'
#' \donttest{
#' if (requireNamespace("ranger", quietly = TRUE) &&
#'     requireNamespace("rpart", quietly = TRUE)) {
#'   s <- rapid_suite(X, Y, key_vars = c("age", "sex"), target_var = "income",
#'                    model_types = c("rf", "cart", "lookup"), n_perm = 3, seed = 1)
#'   summary(s)
#'   s[["rf"]]
#'   plot(s)
#' }
#' }
#'
#' @export
rapid_suite <- function(X, ...) {
  if (missing(X)) {
    # Named rapid objects, or the data under their alias names
    objects <- list(...)
    if (length(objects) > 0 && all(vapply(objects, inherits, logical(1), what = "rapid"))) {
      return(do.call(rapid_suite.rapid, objects))
    }
    cl <- sys.call()
    cl[[1L]] <- rapid_suite.default
    return(eval(cl, parent.frame()))
  }
  UseMethod("rapid_suite")
}

#' @rdname rapid_suite
#' @export
rapid_suite.synth_pair <- function(X, ...) {
  # The default method unpacks the pair through .rapid_inputs()
  NextMethod()
}

#' @rdname rapid_suite
#' @export
rapid_suite.default <- function(X, Y, key_vars, target_var,
                                model_types = c("rf", "cart", "gbm", "lookup"), ...,
                                seed = NULL) {
  inputs <- .rapid_inputs(
    if (missing(X)) NULL else X, if (missing(Y)) NULL else Y,
    if (missing(key_vars)) NULL else key_vars,
    if (missing(target_var)) NULL else target_var,
    list(...)
  )
  if (!is.character(model_types) || length(model_types) == 0 || anyNA(model_types) ||
      anyDuplicated(model_types) > 0 || !all(model_types %in% .rapid_learner_names)) {
    .rapid_abort(c("'model_types' must be distinct names out of: ",
                   paste0("\"", .rapid_learner_names, "\"", collapse = ", "), "."), "input")
  }

  # Decide before fitting which models can run
  is_categorical <- .rapid_var_kind(inputs$X[[inputs$target_var]], inputs$target_var,
                                    "X") == "categorical"
  skipped <- character(0)
  for (m in model_types) {
    why <- tryCatch({
      .learner_check(.rapid_learner(m), is_categorical, "auto")
      NULL
    }, rapid_error_dependency = conditionMessage, rapid_error_model = conditionMessage)
    if (!is.null(why)) skipped[m] <- why
  }
  usable <- setdiff(model_types, names(skipped))
  if (length(skipped) > 0) {
    msg <- paste0("Skipped model(s) ", paste(names(skipped), collapse = ", "), ": ",
                  paste(skipped, collapse = " "))
    if (length(usable) == 0) .rapid_abort(c("No usable model. ", msg), "model")
    .rapid_warn(msg, "skipped_model")
  }

  fits <- lapply(usable, function(m) {
    do.call(rapid, c(list(inputs$X, inputs$Y, key_vars = inputs$key_vars,
                          target_var = inputs$target_var, model_type = m, seed = seed),
                     inputs$dots))
  })
  names(fits) <- usable
  .rapid_suite_new(fits)
}

#' @rdname rapid_suite
#' @export
rapid_suite.rapid <- function(X, ...) {
  fits <- c(if (!missing(X)) list(X), list(...))
  fit_names <- names(fits)
  if (is.null(fit_names)) fit_names <- rep("", length(fits))
  if (!all(vapply(fits, inherits, logical(1), what = "rapid"))) {
    .rapid_abort("Every element of 'X' and '...' must be a rapid object.", "object")
  }
  for (f in fits) .rapid_require(f, "version")
  # The first object has no name in the call: name it by its model
  fit_names[!nzchar(fit_names)] <- vapply(fits[!nzchar(fit_names)], function(f) f$model_type,
                                          character(1))
  if (anyDuplicated(fit_names) > 0) {
    .rapid_abort(c("Two objects have the name '", fit_names[anyDuplicated(fit_names)],
                   "'. Name them in the call, e.g. rapid_suite(a = fit1, b = fit2)."),
                 "object")
  }
  names(fits) <- fit_names
  first <- fits[[1L]]
  for (nm in names(fits)[-1L]) {
    f <- fits[[nm]]
    differ <- c(
      target_var = !identical(f$target_var, first$target_var),
      key_vars = !setequal(f$key_vars, first$key_vars),
      n_total = !identical(as.integer(f$n_total), as.integer(first$n_total)),
      rule = !identical(format(f$rule), format(first$rule))
    )
    if (any(differ)) {
      .rapid_abort(c("The objects '", names(fits)[1L], "' and '", nm, "' differ in ",
                     paste(names(differ)[differ], collapse = ", "), "."), "object")
    }
  }
  .rapid_suite_new(fits)
}

#' Assemble the suite from named, compatible rapid objects
#' @noRd
.rapid_suite_new <- function(fits) {
  num <- function(field) vapply(fits, function(f) as.numeric(f[[field]]), numeric(1))
  results <- data.frame(model_type = names(fits), rapid = unname(num("rapid")),
                        chance = unname(num("chance")), excess = unname(num("excess")),
                        n_at_risk = unname(num("n_at_risk")), stringsAsFactors = FALSE)
  first <- fits[[1L]]
  has_excess <- !all(is.na(results$excess))
  flagged <- Reduce(`|`, lapply(fits, function(f) f$at_risk))
  notes <- lapply(fits, function(f) f$note)
  notes <- notes[!vapply(notes, is.null, logical(1))]
  structure(list(
    results = results,
    fits = fits,
    rapid_mean = mean(results$rapid),
    rapid_max = max(results$rapid),
    worst_model = results$model_type[which.max(results$rapid)],
    excess_max = if (has_excess) max(results$excess, na.rm = TRUE) else NA_real_,
    worst_model_excess = if (has_excess) results$model_type[which.max(results$excess)]
                         else NA_character_,
    n_flagged_any = sum(flagged),
    n_total = first$n_total,
    key_vars = first$key_vars,
    target_var = first$target_var,
    rule = first$rule,
    notes = notes
  ), class = "rapid_suite")
}

#' @rdname rapid_suite
#' @param x,object an object of class \code{"rapid_suite"}
#' @param digits number of decimal places
#' @export
format.rapid_suite <- function(x, digits = 3, ...) {
  f <- function(v) formatC(v, digits = digits, format = "f")
  r <- x$results
  tab <- data.frame(model = r$model_type, RAPID = f(r$rapid), chance = f(r$chance),
                    excess = f(r$excess), flagged = r$n_at_risk, stringsAsFactors = FALSE)
  c(
    paste0("RAPID model suite (", nrow(r), if (nrow(r) == 1) " model, " else " models, ",
           "target '", x$target_var, "', ", x$n_total, " records)"),
    paste0("Rule: ", format(x$rule)),
    utils::capture.output(print(tab, row.names = FALSE, right = FALSE)),
    paste0("Mean RAPID = ", f(x$rapid_mean), ", maximum = ", f(x$rapid_max),
           " (", x$worst_model, ")"),
    if (!is.na(x$excess_max)) {
      paste0("Largest excess over the chance level = ", f(x$excess_max), " (",
             x$worst_model_excess, ")")
    },
    paste0("Records flagged by at least one model: ", x$n_flagged_any, " / ", x$n_total),
    if (length(x$notes) > 0) {
      paste0("Note (", paste(names(x$notes), collapse = ", "), "): ",
             paste(unique(unlist(x$notes)), collapse = " "),
             " Compare the models by their excess.")
    }
  )
}

#' @rdname rapid_suite
#' @export
print.rapid_suite <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname rapid_suite
#' @export
summary.rapid_suite <- function(object, ...) {
  flags <- lapply(object$fits, function(f) f$at_risk)
  all_flagged <- Reduce(`&`, flags)
  structure(list(
    results = object$results,
    rapid_mean = object$rapid_mean,
    rapid_max = object$rapid_max,
    worst_model = object$worst_model,
    excess_max = object$excess_max,
    worst_model_excess = object$worst_model_excess,
    n_flagged_any = object$n_flagged_any,
    n_flagged_all = sum(all_flagged),
    n_total = object$n_total,
    key_vars = object$key_vars,
    target_var = object$target_var,
    rule = format(object$rule),
    notes = object$notes
  ), class = "summary.rapid_suite")
}

#' @rdname rapid_suite
#' @export
format.summary.rapid_suite <- function(x, ...) {
  f <- function(v) formatC(v, digits = 4, format = "f")
  shown <- x$results
  shown[c("rapid", "chance", "excess")] <- lapply(shown[c("rapid", "chance", "excess")], f)
  c(
    "Summary: RAPID Model Suite",
    "==========================",
    paste0("Target:         ", x$target_var),
    paste0("Quasi-identifiers: ", paste(x$key_vars, collapse = ", ")),
    paste0("Rule:           ", x$rule),
    paste0("Records:        ", x$n_total),
    "",
    utils::capture.output(print(shown, row.names = FALSE)),
    "",
    paste0("Mean RAPID:     ", f(x$rapid_mean)),
    paste0("Maximum RAPID:  ", f(x$rapid_max), " (", x$worst_model, ")"),
    if (!is.na(x$excess_max)) {
      paste0("Largest excess: ", f(x$excess_max), " (", x$worst_model_excess, ")")
    },
    paste0("Flagged by at least one model: ", x$n_flagged_any, " / ", x$n_total),
    paste0("Flagged by every model:        ", x$n_flagged_all, " / ", x$n_total),
    if (length(x$notes) > 0) {
      paste0("Note (", paste(names(x$notes), collapse = ", "), "): ",
             paste(unique(unlist(x$notes)), collapse = " "))
    }
  )
}

#' @rdname rapid_suite
#' @export
print.summary.rapid_suite <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname rapid_suite
#' @export
plot.rapid_suite <- function(x, ...) {
  r <- x$results
  d <- data.frame(
    model = factor(rep(r$model_type, 2), levels = r$model_type[order(r$rapid)]),
    measure = factor(rep(c("RAPID", "Chance level"), each = nrow(r)),
                     levels = c("RAPID", "Chance level")),
    value = c(r$rapid, r$chance)
  )
  d <- d[!is.na(d$value), , drop = FALSE]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$model,
                                        colour = .data$measure, shape = .data$measure)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(values = c(RAPID = "firebrick", "Chance level" = "grey40")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Share of records flagged", y = NULL, colour = NULL, shape = NULL,
                  title = "RAPID by attacker model",
                  subtitle = paste0("Target: ", x$target_var))
  print(p)
  invisible(p)
}

#' @rdname rapid_suite
#' @param row.names,optional not used
#' @export
as.data.frame.rapid_suite <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$results
}

#' @rdname rapid_suite
#' @param i name or position of the model
#' @export
`[[.rapid_suite` <- function(x, i, ...) {
  unclass(x)$fits[[i]]
}

#' Risk Contribution of Quasi-Identifiers to RAPID
#'
#' Quantifies how much each quasi-identifier contributes to the RAPID
#' disclosure risk under a declared attack scenario.  Starting from the
#' full set of quasi-identifiers \eqn{Q}, each quasi-identifier \eqn{j}
#' (or a group of them) is removed from the attacker's knowledge, the
#' attacker model is refitted on the synthetic data without it, and the
#' drop in risk is recorded:
#' \deqn{\Delta_j = \mathrm{RAPID}(Q) - \mathrm{RAPID}(Q \setminus \{j\}).}
#'
#' Because part of every RAPID value is present even when the predictions
#' carry no information about the true values, every scenario also reports
#' the chance level of \code{\link{rapid}} (the share of records flagged
#' when the true values are permuted) and the function reports the change
#' in the excess over it,
#' \eqn{\Delta^{exc}_j = [\mathrm{RAPID}(Q) - c(Q)] -
#' [\mathrm{RAPID}(Q \setminus \{j\}) - c(Q \setminus \{j\})]}.
#' This is the quantity to interpret.
#'
#' @details
#' Every row of the result comes from one call of \code{\link{rapid}} with
#' \code{n_perm = n_permutations}, so that the arguments in \code{...}
#' (model, scoring rule, thresholds) define the attack in the same way as
#' there.
#'
#' \strong{Interpretation.} \eqn{\Delta_j} answers "how much risk would
#' disappear if the attacker did not know \eqn{j}?", which is also the
#' first-order effect of withholding or strongly coarsening \eqn{j}.  It is
#' a drop-one (drop-column) quantity: the attacker model is refitted, so it
#' measures the value of the \emph{knowledge}, not how much a fixed model
#' relies on a variable (permutation importance).
#'
#' \strong{Substitutes.} Quasi-identifiers that carry the same information
#' mask each other: dropping either alone changes little, while dropping
#' both changes a lot.  Use \code{groups} to drop such sets jointly
#' (e.g. relationship and marital status in the Adult data).  A small
#' \eqn{\Delta_j} therefore means "not needed given the others", not
#' "harmless".
#'
#' \strong{Uncertainty.} A single synthetic release gives one value per
#' quasi-identifier. Differences of a few thousandths are within model
#' noise.  Repeat over several releases (or seeds) and report the spread.
#' If \code{seed} is given, every attacker model and every set of
#' permutations starts from the same seed (common random numbers), which
#' reduces this noise.
#'
#' @param X The original (confidential) data: a data frame, or a
#'   \code{\link{synth_pair}}.
#' @param Y The synthetic data used to train the attacker model.
#' @param key_vars Character vector: the full set of quasi-identifiers
#'   known to the attacker in the reference scenario (at least two).
#' @param target_var Character string naming the sensitive (target)
#'   variable.
#' @param groups Optional named list of character vectors; each element is
#'   a set of quasi-identifiers that is dropped jointly, in addition to the
#'   single drops.  Names are used as labels.
#' @param n_permutations Number of permutations of the true values for the
#'   chance level of each attacker model (default 99).  Set to 0 to skip
#'   the chance level. The excess columns are then \code{NA}.
#' @param seed Optional random seed (see Details).
#' @param verbose Logical; print progress?
#' @param ... Further arguments passed to \code{\link{rapid}}: the model
#'   (\code{model_type}), the scoring rule and its thresholds.  The old
#'   argument names \code{original_data}, \code{synthetic_data},
#'   \code{quasi_identifiers} and \code{sensitive_attribute} are accepted
#'   as aliases of \code{X}, \code{Y}, \code{key_vars} and
#'   \code{target_var}.
#'
#' @return An object of class \code{"rapid_contribution"}, a list with
#' \describe{
#'   \item{results}{Data frame with one row for the full scenario
#'     (\code{dropped = "(none)"}) and one per dropped quasi-identifier or
#'     group: \code{dropped}, \code{n_qi}, \code{rapid}, \code{chance},
#'     \code{excess}, \code{delta_rapid}, \code{delta_excess}.  Rows are
#'     sorted by decreasing contribution.}
#'   \item{reference}{The row of the full scenario.}
#'   \item{rapid, chance, excess}{Values of the full scenario.}
#'   \item{rule}{The scoring rule of the full scenario.}
#'   \item{key_vars, target_var, groups, is_categorical,
#'     n_permutations}{Settings used.}
#' }
#' \code{format()} returns the printed lines and \code{print()} the object,
#' invisibly. \code{summary()} returns an object of class
#' \code{"summary.rapid_contribution"} with the full table and the
#' quasi-identifier of the largest contribution, printed by its own
#' \code{print()} method. \code{as.data.frame()} returns \code{results}.
#' \code{plot()} draws the drop for each scenario as a lollipop chart (groups
#' as triangles) and returns the ggplot object, or a list of two for both
#' measures, invisibly.
#'
#' @seealso \code{\link{rapid}}, \code{\link{rapid_test}}
#' @family rapid
#' @author Matthias Templ
#' @examples
#' set.seed(1)
#' n <- 300
#' X <- data.frame(age = sample(20:70, n, replace = TRUE),
#'                 sex = factor(sample(c("F", "M"), n, replace = TRUE)),
#'                 region = factor(sample(letters[1:4], n, replace = TRUE)))
#' X$smoker <- factor(ifelse(X$age < 40 & X$sex == "M", "yes", "no"))
#' Y <- X
#' Y$age <- Y$age + sample(-2:2, n, replace = TRUE)
#' if (requireNamespace("rpart", quietly = TRUE)) {
#'   res <- rapid_qi_contribution(X, Y,
#'                                key_vars = c("age", "sex", "region"),
#'                                target_var = "smoker",
#'                                model_type = "cart", n_permutations = 19,
#'                                seed = 1)
#'   print(res)
#'   print(summary(res))
#'   print(as.data.frame(res))
#'   plot(res)
#' }
#'
#' # The linear model drops a group of quasi-identifiers as well
#' X$income <- round(exp(9 + 0.03 * X$age + 0.3 * (X$sex == "M") + rnorm(n, 0, 0.1)))
#' Y$income <- round(X$income * exp(rnorm(n, 0, 0.05)))
#' res_g <- rapid_qi_contribution(X, Y,
#'                                key_vars = c("age", "sex", "region"),
#'                                target_var = "income",
#'                                groups = list(demography = c("age", "sex")),
#'                                model_type = "lm", n_permutations = 19,
#'                                seed = 1)
#' res_g
#' plot(res_g, measure = c("excess", "rapid"))
#' @export
rapid_qi_contribution <- function(X, ...) {
  if (missing(X)) {
    # Without X the data arrive under their alias names
    cl <- sys.call()
    cl[[1L]] <- rapid_qi_contribution.default
    return(eval(cl, parent.frame()))
  }
  UseMethod("rapid_qi_contribution")
}

#' @rdname rapid_qi_contribution
#' @export
rapid_qi_contribution.synth_pair <- function(X, ...) {
  # The default method unpacks the pair
  NextMethod()
}

#' @rdname rapid_qi_contribution
#' @export
rapid_qi_contribution.default <- function(X, Y, key_vars, target_var, groups = NULL,
                                          n_permutations = 99, seed = NULL,
                                          verbose = FALSE, ...) {
  inputs <- .rapid_inputs(
    if (missing(X)) NULL else X, if (missing(Y)) NULL else Y,
    if (missing(key_vars)) NULL else key_vars,
    if (missing(target_var)) NULL else target_var,
    list(...)
  )
  .rapid_check_reserved(inputs$dots, c(n_perm = paste(
    "use 'n_permutations' for the permutations of each scenario")),
    "rapid_qi_contribution")
  key_vars <- inputs$key_vars
  if (length(key_vars) < 2) {
    .rapid_abort("At least two quasi-identifiers are needed to drop one.", "input")
  }
  if (!is.numeric(n_permutations) || length(n_permutations) != 1 ||
      !is.finite(n_permutations) ||
      n_permutations < 0 || n_permutations != round(n_permutations)) {
    .rapid_abort("'n_permutations' must be a single whole number (0 skips the chance level).",
                 "input")
  }
  if (!is.null(groups)) {
    if (!is.list(groups) || is.null(names(groups)) || any(names(groups) == "")) {
      .rapid_abort("'groups' must be a named list of character vectors.", "input")
    }
    bad <- setdiff(unlist(groups), key_vars)
    if (length(bad)) {
      .rapid_abort(c("'groups' contains variables not in 'key_vars': ",
                     paste(bad, collapse = ", ")), "input")
    }
    if (any(vapply(groups, function(g) length(setdiff(key_vars, g)) == 0, logical(1)))) {
      .rapid_abort("A group may not contain all quasi-identifiers.", "input")
    }
  }

  drops <- c(list("(none)" = character(0)),
             stats::setNames(as.list(key_vars), key_vars),
             groups)

  fits <- list()
  one <- function(label, drop) {
    qi <- setdiff(key_vars, drop)
    if (verbose) message("  dropped: ", label)
    fit <- do.call(rapid, c(
      list(X = inputs$X, Y = inputs$Y, key_vars = qi, target_var = inputs$target_var,
           n_perm = n_permutations, seed = seed, verbose = FALSE),
      inputs$dots))
    if (length(fits) == 0) fits <<- list(fit)
    data.frame(dropped = label, n_qi = length(qi), rapid = fit$rapid,
               chance = fit$chance, excess = fit$excess, stringsAsFactors = FALSE)
  }

  res <- do.call(rbind, Map(one, names(drops), drops))
  rownames(res) <- NULL
  ref <- res[1, ]
  res$delta_rapid <- ref$rapid - res$rapid
  res$delta_excess <- ref$excess - res$excess
  key <- if (n_permutations > 0) res$delta_excess else res$delta_rapid
  res <- res[c(1, 1 + order(key[-1], decreasing = TRUE)), ]
  rownames(res) <- NULL
  ref <- res[1, ]

  structure(list(results = res, reference = ref,
                 rapid = ref$rapid, chance = ref$chance, excess = ref$excess,
                 rule = fits[[1]]$rule,
                 key_vars = key_vars, target_var = inputs$target_var,
                 groups = groups, is_categorical = fits[[1]]$is_categorical,
                 n_permutations = n_permutations),
            class = "rapid_contribution")
}


#' @rdname rapid_qi_contribution
#' @param x,object An object of class \code{"rapid_contribution"}.
#' @param digits Number of digits to print.
#' @export
format.rapid_contribution <- function(x, digits = 3, ...) {
  f <- function(v) formatC(v, digits = digits, format = "f")
  tab <- x$results[-1, c("dropped", "n_qi", "rapid", "delta_rapid",
                         if (x$n_permutations > 0) c("excess", "delta_excess"))]
  num <- vapply(tab, is.double, logical(1))
  tab[num] <- lapply(tab[num], round, digits)
  c(
    "RAPID risk contribution of quasi-identifiers (drop-one, refitted)",
    paste0("Target: ", x$target_var, "; rule: ", format(x$rule)),
    paste0("Full scenario: ", length(x$key_vars), " quasi-identifiers, RAPID = ", f(x$rapid),
           if (x$n_permutations > 0) {
             paste0(", chance level = ", f(x$chance), ", excess = ", f(x$excess))
           }),
    utils::capture.output(print(tab, row.names = FALSE)),
    if (x$n_permutations > 0) {
      "delta_excess = drop in the excess over the chance level."
    },
    "Next: summary(), plot()"
  )
}

#' @rdname rapid_qi_contribution
#' @export
print.rapid_contribution <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname rapid_qi_contribution
#' @export
summary.rapid_contribution <- function(object, ...) {
  structure(list(
    target_var = object$target_var,
    key_vars = object$key_vars,
    rule = format(object$rule),
    n_permutations = object$n_permutations,
    rapid = object$rapid,
    chance = object$chance,
    excess = object$excess,
    top = if (nrow(object$results) > 1) object$results$dropped[2],
    table = object$results
  ), class = "summary.rapid_contribution")
}

#' @rdname rapid_qi_contribution
#' @export
format.summary.rapid_contribution <- function(x, digits = 4, ...) {
  tab <- x$table
  num <- vapply(tab, is.double, logical(1))
  tab[num] <- lapply(tab[num], round, digits)
  c(
    "Summary: RAPID Risk Contribution of Quasi-Identifiers",
    "=====================================================",
    paste0("Target variable: ", x$target_var),
    paste0("Rule: ", x$rule),
    paste0("Quasi-identifiers: ", paste(x$key_vars, collapse = ", ")),
    paste0("Permutations for the chance level: ", x$n_permutations),
    "Scenarios (first row: all quasi-identifiers known):",
    utils::capture.output(print(tab, row.names = FALSE)),
    if (!is.null(x$top)) paste0("Largest contribution: ", x$top)
  )
}

#' @rdname rapid_qi_contribution
#' @export
print.summary.rapid_contribution <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname rapid_qi_contribution
#' @param row.names,optional Not used.
#' @export
as.data.frame.rapid_contribution <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$results
}

#' @rdname rapid_qi_contribution
#' @param measure \code{"excess"} (default; change in excess over the
#'   chance level) or \code{"rapid"} (change in raw RAPID). Both may be
#'   given and are drawn in this order.
#' @export
plot.rapid_contribution <- function(x, measure = "excess", ...) {
  measure <- match.arg(measure, c("excess", "rapid"), several.ok = TRUE)
  if (x$n_permutations == 0) measure <- unique(sub("excess", "rapid", measure))
  plots <- lapply(measure, function(m) {
    df <- x$results[-1, ]
    df$value <- if (m == "excess") df$delta_excess else df$delta_rapid
    df$group <- ifelse(df$dropped %in% names(x$groups), "group", "single")
    df$dropped <- factor(df$dropped, levels = rev(df$dropped[order(-df$value)]))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$dropped)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$value,
                                         yend = .data$dropped),
                            colour = "grey70") +
      ggplot2::geom_point(ggplot2::aes(shape = .data$group), size = 2.5,
                          colour = "steelblue") +
      ggplot2::scale_shape_manual(values = c(single = 16, group = 17), guide = "none") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = if (m == "excess") "Drop in RAPID excess over the chance level"
            else "Drop in RAPID",
        y = "Quasi-identifier(s) unknown to the attacker",
        title = "Risk contribution of quasi-identifiers",
        subtitle = sprintf("Target: %s", x$target_var))
    print(p)
    p
  })
  invisible(if (length(plots) == 1) plots[[1]] else plots)
}

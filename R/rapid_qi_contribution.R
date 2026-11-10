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
#' Because part of every RAPID value is a confidence component that is
#' present even when the predictions carry no information about the true
#' labels, the function also evaluates the fixed-model permutation null of
#' \code{\link{rapid_test}} for each attacker model and reports the change
#' in the excess over the null,
#' \eqn{\Delta^{exc}_j = [\mathrm{RAPID}(Q) - \bar{R}_0(Q)] -
#' [\mathrm{RAPID}(Q \setminus \{j\}) - \bar{R}_0(Q \setminus \{j\})]}.
#' This is the quantity to interpret.
#'
#' @details
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
#' quasi-identifier; differences of a few thousandths are within model
#' noise.  Repeat over several releases (or seeds) and report the spread.
#' If \code{seed} is given, every attacker model is fitted with the same
#' seed (common random numbers), which reduces this noise.
#'
#' @param original_data A data frame of original (confidential) records.
#' @param synthetic_data A data frame of synthetic records used to train
#'   the attacker model.
#' @param quasi_identifiers Character vector: the full set of
#'   quasi-identifiers known to the attacker in the reference scenario.
#' @param sensitive_attribute Character string naming the sensitive
#'   (target) variable.
#' @param groups Optional named list of character vectors; each element is
#'   a set of quasi-identifiers that is dropped jointly, in addition to the
#'   single drops.  Names are used as labels.
#' @param model_type Attacker model; see \code{\link{rapid}}.
#' @param cat_tau,num_epsilon,cat_eval_method,num_error_metric,num_epsilon_type,num_delta
#'   Risk definition; see \code{\link{rapid}}.  As in
#'   \code{\link{rapid_test}}, the default \code{cat_eval_method} is
#'   \code{"RCS_marginal"} (normalized gain over the marginal class
#'   proportion).
#' @param n_permutations Number of label permutations for the null mean
#'   of each attacker model (default 99).  Set to 0 to skip the null; the
#'   excess columns are then \code{NA}.
#' @param seed Optional random seed (see Details).
#' @param verbose Logical; print progress?
#' @param ... Additional arguments passed to the model fitting function.
#'
#' @return An object of class \code{"rapid_contribution"}, a list with
#' \describe{
#'   \item{results}{Data frame with one row for the full scenario
#'     (\code{dropped = "(none)"}) and one per dropped quasi-identifier or
#'     group: \code{dropped}, \code{n_qi}, \code{rapid}, \code{null_mean},
#'     \code{excess}, \code{delta_rapid}, \code{delta_excess}.  Rows are
#'     sorted by decreasing contribution.}
#'   \item{reference}{The row of the full scenario.}
#'   \item{quasi_identifiers, sensitive_attribute, groups, cat_tau,
#'     num_epsilon, is_categorical, n_permutations}{Settings used.}
#' }
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
#' res <- rapid_qi_contribution(X, Y,
#'                              quasi_identifiers = c("age", "sex", "region"),
#'                              sensitive_attribute = "smoker",
#'                              model_type = "cart", n_permutations = 19,
#'                              seed = 1)
#' res
#' \donttest{
#' plot(res)
#' }
#' @export
rapid_qi_contribution <- function(original_data, synthetic_data,
                                  quasi_identifiers, sensitive_attribute,
                                  groups = NULL,
                                  model_type = "rf",
                                  cat_tau = 0.3, num_epsilon = 10,
                                  cat_eval_method = "RCS_marginal",
                                  num_error_metric = "symmetric",
                                  num_epsilon_type = "percentage",
                                  num_delta = 0.01,
                                  n_permutations = 99,
                                  seed = NULL, verbose = FALSE, ...) {
  stopifnot(is.character(quasi_identifiers), length(quasi_identifiers) >= 2,
            !sensitive_attribute %in% quasi_identifiers,
            n_permutations >= 0)
  if (!is.null(groups)) {
    if (!is.list(groups) || is.null(names(groups)) || any(names(groups) == ""))
      stop("'groups' must be a named list of character vectors.")
    bad <- setdiff(unlist(groups), quasi_identifiers)
    if (length(bad))
      stop("'groups' contains variables not in 'quasi_identifiers': ",
           paste(bad, collapse = ", "))
    if (any(vapply(groups, function(g) length(setdiff(quasi_identifiers, g)) == 0,
                   logical(1))))
      stop("A group may not contain all quasi-identifiers.")
  }

  drops <- c(list("(none)" = character(0)),
             stats::setNames(as.list(quasi_identifiers), quasi_identifiers),
             groups)

  is_categorical <- NA
  one <- function(label, drop) {
    qi <- setdiff(quasi_identifiers, drop)
    if (verbose) message("  dropped: ", label)
    if (!is.null(seed)) set.seed(seed)
    prep <- .rapid_prepare_permutation(
      original_data = original_data, synthetic_data = synthetic_data,
      quasi_identifiers = qi, sensitive_attribute = sensitive_attribute,
      model_type = model_type, num_delta = num_delta, ...)
    is_categorical <<- prep$is_categorical
    stat <- function(y) .compute_rapid_permuted(
      predictions = prep$predictions, y = y,
      is_categorical = prep$is_categorical,
      cat_tau = cat_tau, num_epsilon = num_epsilon,
      cat_eval_method = cat_eval_method,
      num_error_metric = num_error_metric,
      num_epsilon_type = num_epsilon_type, num_delta = num_delta,
      original_data = prep$original_prepared,
      sensitive_attribute = sensitive_attribute)
    obs <- stat(prep$y_true)
    null_mean <- if (n_permutations > 0)
      mean(vapply(seq_len(n_permutations), function(b) stat(sample(prep$y_true)),
                  numeric(1)))
    else NA_real_
    data.frame(dropped = label, n_qi = length(qi), rapid = obs,
               null_mean = null_mean, stringsAsFactors = FALSE)
  }

  res <- do.call(rbind, Map(one, names(drops), drops))
  rownames(res) <- NULL
  res$excess <- res$rapid - res$null_mean
  ref <- res[1, ]
  res$delta_rapid <- ref$rapid - res$rapid
  res$delta_excess <- ref$excess - res$excess
  key <- if (n_permutations > 0) res$delta_excess else res$delta_rapid
  res <- res[c(1, 1 + order(key[-1], decreasing = TRUE)), ]
  rownames(res) <- NULL

  structure(list(results = res, reference = ref,
                 quasi_identifiers = quasi_identifiers,
                 sensitive_attribute = sensitive_attribute,
                 groups = groups, cat_tau = cat_tau, num_epsilon = num_epsilon,
                 is_categorical = is_categorical,
                 n_permutations = n_permutations),
            class = "rapid_contribution")
}


#' @rdname rapid_qi_contribution
#' @param x,object An object of class \code{"rapid_contribution"}.
#' @param digits Number of digits to print.
#' @export
print.rapid_contribution <- function(x, digits = 3, ...) {
  thr <- if (isTRUE(x$is_categorical)) sprintf("tau = %s", x$cat_tau)
         else sprintf("epsilon = %s", x$num_epsilon)
  cat("\n  RAPID risk contribution of quasi-identifiers (drop-one, refitted)\n\n")
  cat(sprintf("  Sensitive attribute: %s   (%s)\n", x$sensitive_attribute, thr))
  cat(sprintf("  Full scenario: %d QIs, RAPID = %.*f", length(x$quasi_identifiers),
              digits, x$reference$rapid))
  if (x$n_permutations > 0)
    cat(sprintf(", null mean = %.*f, excess = %.*f", digits, x$reference$null_mean,
                digits, x$reference$excess))
  cat("\n\n")
  tab <- x$results[-1, c("dropped", "n_qi", "rapid", "delta_rapid",
                         if (x$n_permutations > 0) c("excess", "delta_excess"))]
  num <- vapply(tab, is.double, logical(1))
  tab[num] <- lapply(tab[num], round, digits)
  print(tab, row.names = FALSE)
  if (x$n_permutations > 0)
    cat("\n  delta_excess = drop in RAPID beyond the drop of its permutation null.\n")
  cat("\n")
  invisible(x)
}


#' @rdname rapid_qi_contribution
#' @export
summary.rapid_contribution <- function(object, ...) {
  object$results
}


#' @rdname rapid_qi_contribution
#' @param measure \code{"excess"} (default; change in excess over the
#'   permutation null) or \code{"rapid"} (change in raw RAPID).
#' @export
plot.rapid_contribution <- function(x, measure = c("excess", "rapid"), ...) {
  measure <- match.arg(measure)
  if (measure == "excess" && x$n_permutations == 0) measure <- "rapid"
  df <- x$results[-1, ]
  df$value <- if (measure == "excess") df$delta_excess else df$delta_rapid
  df$group <- ifelse(df$dropped %in% names(x$groups), "group", "single")
  df$dropped <- factor(df$dropped, levels = rev(df$dropped[order(-df$value)]))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$dropped)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$value,
                                       yend = .data$dropped),
                          colour = "grey70") +
    ggplot2::geom_point(ggplot2::aes(shape = .data$group), size = 2.5,
                        colour = "steelblue") +
    ggplot2::scale_shape_manual(values = c(single = 16, group = 17), guide = "none") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = if (measure == "excess") "Drop in RAPID excess over the null"
          else "Drop in RAPID",
      y = "Quasi-identifier(s) unknown to the attacker",
      title = "Risk contribution of quasi-identifiers",
      subtitle = sprintf("Sensitive attribute: %s", x$sensitive_attribute))
}

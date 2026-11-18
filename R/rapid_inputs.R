# RAPID: conditions, argument handling and data preparation

# ---- Conditions --------------------------------------------------------------
# Every error and warning raised by RAPID code goes through these two functions,
# so that callers can catch them by class.

#' Signal a classed RAPID error
#'
#' @param message Character; pasted together without separator.
#' @param class Suffix of the condition class (`"input"`, `"model"`,
#'   `"levels"`, `"tolerance"`, `"dependency"`, `"object"`).
#' @noRd
.rapid_abort <- function(message, class = NULL) {
  stop(structure(
    class = c(if (!is.null(class)) paste0("rapid_error_", class),
              "rapid_error", "error", "condition"),
    list(message = paste(message, collapse = ""), call = NULL)
  ))
}

#' Signal a classed RAPID warning
#'
#' @param message Character; pasted together without separator.
#' @param class Suffix of the condition class (`"deprecated"`,
#'   `"ignored_argument"`, `"legacy_tolerance"`, `"skipped_model"`).
#' @noRd
.rapid_warn <- function(message, class = NULL) {
  warning(structure(
    class = c(if (!is.null(class)) paste0("rapid_warning_", class),
              "rapid_warning", "warning", "condition"),
    list(message = paste(message, collapse = ""), call = NULL)
  ))
}

# ---- Choices ---------------------------------------------------------------------

#' Match an argument against its choices, with a classed error
#'
#' Like `match.arg()`: the untouched default (all choices) gives the first one,
#' partial matching is allowed.
#' @noRd
.rapid_choice <- function(arg, choices, name) {
  if (length(arg) > 1 && identical(arg, choices)) return(choices[1])
  i <- if (is.character(arg) && length(arg) == 1) pmatch(arg, choices) else NA
  if (is.na(i)) {
    .rapid_abort(c("'", name, "' must be one of: ",
                   paste0("\"", choices, "\"", collapse = ", "), "."), "input")
  }
  choices[i]
}

# ---- Data arguments ------------------------------------------------------------

# Old argument names (RAPID paper, RAPID package) and the names used here
.rapid_aliases <- c(original_data = "X", synthetic_data = "Y",
                    quasi_identifiers = "key_vars", sensitive_attribute = "target_var")

#' Replace the old argument names in `dots` by the names used here
#'
#' @param arg Named list of the data arguments (`NULL` when not given).
#' @param dots The extra arguments of the call.
#' @return `list(arg, dots)`; `dots` without the alias names.
#' @noRd
.rapid_resolve_aliases <- function(arg, dots) {
  for (alias in intersect(names(.rapid_aliases), names(dots))) {
    name <- .rapid_aliases[[alias]]
    if (!(name %in% names(arg))) next
    if (!is.null(arg[[name]])) {
      .rapid_abort(c("'", name, "' and its alias '", alias, "' are both given. ",
                     "Use '", name, "' only."), "input")
    }
    arg[name] <- list(dots[[alias]])
    dots[[alias]] <- NULL
  }
  if (length(dots) == 0) dots <- list()
  list(arg = arg, dots = dots)
}

#' Resolve the data arguments of a RAPID function
#'
#' Accepts the four data arguments (`NULL` when not given), the old alias names
#' in `dots`, or a `synth_pair` as `X`, and validates the result.
#'
#' @return `list(X, Y, key_vars, target_var, holdout, dots)`; `dots` without the
#'   alias names.
#' @noRd
.rapid_inputs <- function(X, Y, key_vars, target_var, dots) {
  resolved <- .rapid_resolve_aliases(
    list(X = X, Y = Y, key_vars = key_vars, target_var = target_var), dots)
  arg <- resolved$arg
  dots <- resolved$dots

  holdout <- NULL
  if (inherits(arg$X, "synth_pair")) {
    pair <- arg$X
    if (!is.null(arg$Y)) {
      .rapid_abort("'Y' must not be given when 'X' is a synth_pair.", "input")
    }
    arg$X <- pair$original
    arg$Y <- pair$synthetic
    if (is.null(arg$key_vars)) arg$key_vars <- pair$key_vars
    if (is.null(arg$target_var)) arg$target_var <- pair$target_var
    holdout <- pair$holdout
    for (name in c("key_vars", "target_var")) {
      if (is.null(arg[[name]])) {
        .rapid_abort(c("The synth_pair has no '", name, "'. Set it in synth_pair() or ",
                       "pass '", name, "'."), "input")
      }
    }
  }

  if (is.null(arg$X)) {
    .rapid_abort("The original data are missing: pass 'X' (a data frame or a synth_pair).",
                 "input")
  }
  if (is.null(arg$Y)) {
    .rapid_abort("The synthetic data are missing: pass 'Y' (a data frame).", "input")
  }
  for (name in c("X", "Y")) {
    if (!is.data.frame(arg[[name]])) {
      .rapid_abort(c(name, " must be a data frame."), "input")
    }
    arg[[name]] <- as.data.frame(arg[[name]])
    if (nrow(arg[[name]]) == 0) .rapid_abort(c(name, " has no rows."), "input")
  }
  if (!is.character(arg$key_vars) || length(arg$key_vars) == 0 || anyNA(arg$key_vars)) {
    .rapid_abort("'key_vars' must be a character vector naming the quasi-identifiers.",
                 "input")
  }
  if (!is.character(arg$target_var) || length(arg$target_var) != 1 ||
      is.na(arg$target_var)) {
    .rapid_abort("'target_var' must be the name of one variable.", "input")
  }
  arg$key_vars <- unique(arg$key_vars)
  if (arg$target_var %in% arg$key_vars) {
    .rapid_abort(c("The target variable '", arg$target_var,
                   "' must not be one of the quasi-identifiers."), "input")
  }
  vars <- c(arg$key_vars, arg$target_var)
  miss <- list(X = setdiff(vars, names(arg$X)), Y = setdiff(vars, names(arg$Y)))
  miss <- miss[lengths(miss) > 0]
  if (length(miss) > 0) {
    .rapid_abort(paste(vapply(names(miss), function(d) {
      paste0("Variables missing in ", d, ": ",
             paste0("'", miss[[d]], "'", collapse = ", "), ".")
    }, character(1)), collapse = " "), "input")
  }

  list(X = arg$X, Y = arg$Y, key_vars = arg$key_vars, target_var = arg$target_var,
       holdout = holdout, dots = dots)
}

# ---- Extra arguments -------------------------------------------------------------

#' Check the arguments that are passed on to the model
#'
#' Unnamed elements are an error. So is a name that is a near miss of one of the
#' `known` arguments (edit distance of at most 2 and less than half the length
#' of the known name, so that short model arguments such as rpart's `cp` or
#' `y` are not mistaken for `X` or `Y`).
#' @noRd
.rapid_check_dots <- function(dots, known) {
  if (length(dots) == 0) return(invisible(TRUE))
  nms <- names(dots)
  if (is.null(nms) || any(!nzchar(nms))) {
    .rapid_abort(c("Unnamed extra arguments are not allowed: every argument after ",
                   "'target_var' must be named."), "input")
  }
  for (nm in setdiff(nms, known)) {
    d <- drop(utils::adist(nm, known))
    near <- d <= 2 & d < nchar(known) / 2
    if (any(near)) {
      .rapid_abort(c("Unknown argument '", nm, "'. Did you mean '",
                     known[near][which.min(d[near])], "'?"), "input")
    }
  }
  invisible(TRUE)
}

#' Refuse extra arguments in a method on a `rapid` object
#' @noRd
.rapid_check_empty_dots <- function(dots, hint) {
  if (length(dots) == 0) return(invisible(TRUE))
  nms <- names(dots)
  if (is.null(nms)) nms <- rep("", length(dots))
  nms[!nzchar(nms)] <- "<unnamed>"
  .rapid_abort(c("Unused argument(s): ", paste(nms, collapse = ", "), ". ", hint), "input")
}

# ---- Data preparation ------------------------------------------------------------

.rapid_var_kind <- function(v, name, data_name) {
  if (is.factor(v) || is.character(v) || is.logical(v)) return("categorical")
  if (is.numeric(v)) return("numeric")
  .rapid_abort(c("Variable '", name, "' in ", data_name, " has the unsupported type '",
                 class(v)[1], "'. Use a numeric, factor, character or logical variable."),
               "input")
}

# Brings one variable of X and Y to the same type. Numeric variables become
# double. Categorical variables become factors; they stay ordered only if both
# are ordered with the same levels and the variable is not the target.
.rapid_harmonise <- function(x, y, name, is_target) {
  kind_x <- .rapid_var_kind(x, name, "X")
  kind_y <- .rapid_var_kind(y, name, "Y")
  if (kind_x != kind_y) {
    .rapid_abort(c("Variable '", name, "' is ", kind_x, " in X and ", kind_y,
                   " in Y. Convert it to the same type in both data sets."), "input")
  }
  if (kind_x == "numeric") {
    return(list(x = as.numeric(x), y = as.numeric(y)))
  }
  keep_order <- !is_target && is.ordered(x) && is.ordered(y) &&
    identical(levels(x), levels(y))
  if (!keep_order) {
    x <- factor(x, levels = if (is.factor(x)) levels(x) else sort(unique(x)), ordered = FALSE)
    y <- factor(y, levels = if (is.factor(y)) levels(y) else sort(unique(y)), ordered = FALSE)
  }
  list(x = x, y = y)
}

#' Prepare original and synthetic data for the attack
#'
#' Keeps the target and the quasi-identifiers, brings both data sets to the same
#' column types, handles missing values and gives all factors the same levels in
#' both sets (levels of the original data first).
#'
#' @return `list(X, Y, is_categorical, dropped, unseen)`: `dropped` are the row
#'   indices of `X` removed by `na_strategy = "drop"`, `unseen` a named list of
#'   the quasi-identifier levels that occur in `X` but not in `Y`.
#' @noRd
.rapid_prepare_data <- function(X, Y, key_vars, target_var, na_strategy,
                                na_constant_value) {
  vars <- c(target_var, key_vars)
  X <- as.data.frame(X)[, vars, drop = FALSE]
  Y <- as.data.frame(Y)[, vars, drop = FALSE]

  for (v in vars) {
    h <- .rapid_harmonise(X[[v]], Y[[v]], v, is_target = v == target_var)
    X[[v]] <- h$x
    Y[[v]] <- h$y
  }
  is_categorical <- is.factor(X[[target_var]])

  # Missing values (a factor gets the level "missing")
  dropped <- if (!is_categorical && na_strategy == "drop") {
    which(is.na(X[[target_var]]))
  } else {
    integer(0)
  }
  na <- .handle_sensitive_na(X, Y, target_var, na_strategy, na_constant_value)
  X <- .handle_qi_na(na$original, key_vars)
  Y <- .handle_qi_na(na$synthetic, key_vars)

  # Common factor levels, those of the original data first
  unseen <- list()
  for (v in vars[vapply(X[vars], is.factor, logical(1))]) {
    if (v != target_var) {
      new <- setdiff(unique(as.character(X[[v]])), unique(as.character(Y[[v]])))
      if (length(new) > 0) unseen[[v]] <- new
    }
    lev <- union(levels(X[[v]]), levels(Y[[v]]))
    X[[v]] <- factor(X[[v]], levels = lev, ordered = is.ordered(X[[v]]))
  }
  Y <- .ensure_levels(Y, X, target_var, key_vars)

  list(X = X, Y = Y, is_categorical = is_categorical, dropped = dropped, unseen = unseen)
}

#' Refuse quasi-identifier levels that the attacker model has never seen
#' @param unseen Named list from `.rapid_prepare_data()`.
#' @param model Name of the model, for the message.
#' @noRd
.rapid_check_levels <- function(unseen, model) {
  if (length(unseen) == 0) return(invisible(TRUE))
  what <- vapply(names(unseen), function(v) {
    paste0("'", v, "' (", paste(unseen[[v]], collapse = ", "), ")")
  }, character(1))
  .rapid_abort(c("The original data contain quasi-identifier levels that do not occur in ",
                 "the synthetic data: ", paste(what, collapse = "; "), ". The model '",
                 model, "' cannot predict such records. Use a tree-based model ",
                 "(\"rf\", \"cart\", \"gbm\") or recode these levels."), "levels")
}

#' Handle NA in sensitive attribute
#' @noRd
.handle_sensitive_na <- function(original, synthetic, target_var, strategy, constant_value) {
  orig_target <- original[[target_var]]
  synth_target <- synthetic[[target_var]]

  has_na_orig <- anyNA(orig_target)
  has_na_synth <- anyNA(synth_target)

  if (!has_na_orig && !has_na_synth) {
    return(list(original = original, synthetic = synthetic))
  }

  if (is.factor(orig_target)) {
    # For factors: add "missing" level
    if (has_na_orig) {
      original[[target_var]] <- addNA(original[[target_var]])
      levels(original[[target_var]])[is.na(levels(original[[target_var]]))] <- "missing"
    }
    if (has_na_synth) {
      synthetic[[target_var]] <- addNA(synthetic[[target_var]])
      levels(synthetic[[target_var]])[is.na(levels(synthetic[[target_var]]))] <- "missing"
    }
  } else if (is.numeric(orig_target)) {
    # For numeric: apply strategy
    if (strategy == "drop") {
      keep_orig <- !is.na(original[[target_var]])
      keep_synth <- !is.na(synthetic[[target_var]])
      original <- original[keep_orig, , drop = FALSE]
      synthetic <- synthetic[keep_synth, , drop = FALSE]
    } else if (strategy == "constant") {
      original[[target_var]][is.na(original[[target_var]])] <- constant_value
      synthetic[[target_var]][is.na(synthetic[[target_var]])] <- constant_value
    } else if (strategy == "median") {
      orig_med <- median(original[[target_var]], na.rm = TRUE)
      synth_med <- median(synthetic[[target_var]], na.rm = TRUE)
      original[[target_var]][is.na(original[[target_var]])] <- orig_med
      synthetic[[target_var]][is.na(synthetic[[target_var]])] <- synth_med
    }
  }

  list(original = original, synthetic = synthetic)
}

#' Handle NA in quasi-identifiers
#' @noRd
.handle_qi_na <- function(data, qi_vars) {
  for (var in qi_vars) {
    if (!var %in% names(data)) next

    if (is.factor(data[[var]]) && anyNA(data[[var]])) {
      data[[var]] <- addNA(data[[var]])
      levels(data[[var]])[is.na(levels(data[[var]]))] <- "missing"
    } else if (is.numeric(data[[var]]) && anyNA(data[[var]])) {
      .rapid_abort(c("Variable '", var, "' is numeric and contains NA. Please handle ",
                     "missing values in the quasi-identifiers before using RAPID."),
                   "input")
    }
  }
  data
}

#' Ensure factor levels are consistent between datasets
#' @noRd
.ensure_levels <- function(synthetic, original, target_var, key_vars, verbose = FALSE) {
  # Convert character to factor if needed
  if (is.character(synthetic[[target_var]])) {
    synthetic[[target_var]] <- factor(synthetic[[target_var]])
  }
  if (is.character(original[[target_var]])) {
    original[[target_var]] <- factor(original[[target_var]])
  }

  # For categorical target, ensure all original levels exist in synthetic

  if (is.factor(original[[target_var]])) {
    orig_levels <- levels(original[[target_var]])
    synth_levels <- levels(synthetic[[target_var]])

    missing_levels <- setdiff(orig_levels, synth_levels)
    if (length(missing_levels) > 0) {
      # Add missing levels to synthetic
      levels(synthetic[[target_var]]) <- c(synth_levels, missing_levels)
    }

    # Also ensure synthetic levels match original order
    synthetic[[target_var]] <- factor(synthetic[[target_var]], levels = orig_levels)
  }

  # Do the same for categorical key variables
  for (var in key_vars) {
    if (is.factor(original[[var]])) {
      orig_levels <- levels(original[[var]])
      if (is.factor(synthetic[[var]])) {
        synth_levels <- levels(synthetic[[var]])
        all_levels <- union(orig_levels, synth_levels)
        synthetic[[var]] <- factor(synthetic[[var]], levels = all_levels,
                                   ordered = is.ordered(synthetic[[var]]))
      }
    }
  }

  synthetic
}

# ---- Tolerance and thresholds ----------------------------------------------------

#' Tolerance on the internal scale
#'
#' The tolerance of the two relative metrics is a proportion (0.1 is 10 %).
#' Legacy input is still read: `epsilon_type` (deprecated) and, without it, a
#' relative tolerance of 1 or more, which is taken as percent.
#'
#' @param metric Error metric.
#' @param epsilon User tolerance or `NULL`.
#' @param epsilon_type `NULL`, or the deprecated `"percentage"` / `"absolute"`.
#' @return `list(metric, epsilon)`.
#' @noRd
.rapid_resolve_tolerance <- function(metric, epsilon, epsilon_type = NULL) {
  if (!is.null(epsilon) &&
      (!is.numeric(epsilon) || length(epsilon) != 1 || is.na(epsilon) || epsilon <= 0)) {
    .rapid_abort("'num_epsilon' must be a single positive number.", "tolerance")
  }
  if (!is.null(epsilon_type)) {
    epsilon_type <- match.arg(epsilon_type, c("percentage", "absolute"))
    .rapid_warn(c("'num_epsilon_type' is deprecated. Choose the error metric with ",
                  "'num_error_metric' and give a relative tolerance as a proportion ",
                  "(0.1 for 10 %)."), "deprecated")
    if (epsilon_type == "absolute") {
      metric <- "absolute"
    } else {
      if (metric == "absolute") {
        .rapid_abort(c("A percentage tolerance has no meaning for the absolute error ",
                       "metric. Give 'num_epsilon' in the units of the target variable."),
                     "tolerance")
      }
      epsilon <- if (is.null(epsilon)) 0.1 else epsilon / 100
    }
  }
  if (is.null(epsilon)) {
    if (metric == "absolute") {
      .rapid_abort(c("'num_epsilon' must be given for num_error_metric = \"absolute\" ",
                     "(in the units of the target variable)."), "tolerance")
    }
    epsilon <- 0.1
  } else if (is.null(epsilon_type) && metric != "absolute" && epsilon >= 1) {
    .rapid_warn(c("num_epsilon = ", format(epsilon), " is read as ", format(epsilon),
                  " %. A relative tolerance is now a proportion: use num_epsilon = ",
                  format(epsilon / 100), "."), "legacy_tolerance")
    epsilon <- epsilon / 100
  }
  list(metric = metric, epsilon = epsilon)
}

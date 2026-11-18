# RAPID attacker models (learners) and predictive distributions
#
# A *learner* is the model an attacker trains on the synthetic data. Each learner
# is a subclass of `rapid_learner` with methods for fitting, predicting, the
# predictive distribution of a numeric target, variable importance and a
# one-line description. Adding a learner means adding one branch in
# `.rapid_learner()` and its methods here; nothing outside this file may branch
# on the name of a model.

.rapid_has_package <- function(pkg) requireNamespace(pkg, quietly = TRUE)

.rapid_learner_names <- c("rf", "lm", "cart", "gbm", "logit", "lookup")

# Arguments that travel through `...` but belong to a learner, not to its
# backend function. A learner lists those it understands in `options`.
.rapid_learner_options <- "prune"

#' Attacker model
#'
#' @param model_type Name of the model, or a `rapid_learner` (its arguments are
#'   merged with `user_args`).
#' @param user_args List of arguments for the backend function; they override
#'   the defaults of the learner.
#' @return A list of class `c("rapid_learner_<name>", "rapid_learner")` with
#'   `name`, `args`, `requires` (packages), `targets` (supported target types),
#'   `dists` (supported `num_dist` values), `options` (learner arguments it
#'   understands: `"num_scale"`, `"prune"`), `unseen_levels` (can it predict
#'   a record with a factor level that the training data never contained) and
#'   `note` (a line the print of a `rapid` object adds, or `NULL`).
#' @noRd
.rapid_learner <- function(model_type, user_args = list()) {
  if (inherits(model_type, "rapid_learner")) {
    if (length(user_args) == 0) return(model_type)
    learner <- model_type
    learner$args <- utils::modifyList(as.list(learner$args), as.list(user_args))
    return(.learner_take_options(learner))
  }
  both <- c("numeric", "categorical")
  spec <- switch(if (is.character(model_type) && length(model_type) == 1) model_type else "",
    rf = list(requires = "ranger", targets = both, dists = c("auto", "grid", "exact"),
              options = character(0), unseen_levels = TRUE),
    lm = list(requires = character(0), targets = "numeric",
              dists = c("auto", "grid", "exact"), options = "num_scale",
              unseen_levels = FALSE),
    cart = list(requires = "rpart", targets = both, dists = c("auto", "grid", "exact"),
                options = "prune", unseen_levels = TRUE),
    gbm = list(requires = "xgboost", targets = both, dists = c("auto", "grid"),
               options = character(0), unseen_levels = TRUE),
    logit = list(requires = character(0), targets = "categorical", dists = character(0),
                 options = character(0), unseen_levels = FALSE),
    lookup = list(requires = "ranger", targets = both, dists = c("auto", "exact"),
                  options = character(0), unseen_levels = TRUE,
                  note = paste("Lookup attacker: read the excess over the chance level,",
                               "not RAPID itself.")),
    .rapid_abort(c("Unknown model_type. Use one of: ",
                   paste0("\"", .rapid_learner_names, "\"", collapse = ", "), "."), "model")
  )
  learner <- c(list(name = model_type, args = as.list(user_args)), spec)
  if ("prune" %in% spec$options) learner$prune <- TRUE
  learner <- structure(learner,
                       class = c(paste0("rapid_learner_", model_type), "rapid_learner"))
  .learner_take_options(learner)
}

# Moves the learner arguments passed through `...` (`.rapid_learner_options`)
# from `args` into fields of the learner; warns about those it does not know.
.learner_take_options <- function(learner) {
  for (opt in intersect(.rapid_learner_options, names(learner$args))) {
    value <- learner$args[[opt]]
    learner$args[[opt]] <- NULL
    if (!opt %in% learner$options) {
      .rapid_warn(c("Argument '", opt, "' is ignored: it has no meaning for model_type = \"",
                    learner$name, "\"."), "ignored_argument")
    } else {
      if (!is.logical(value) || length(value) != 1 || is.na(value)) {
        .rapid_abort(c("'", opt, "' must be TRUE or FALSE."), "input")
      }
      learner[[opt]] <- value
    }
  }
  if (length(learner$args) == 0) learner$args <- list()
  learner
}

# Names of the built-in learners that support a target type
.rapid_learners_for <- function(type, exclude = character(0)) {
  ok <- vapply(.rapid_learner_names, function(m) type %in% .rapid_learner(m)$targets,
               logical(1))
  paste(setdiff(.rapid_learner_names[ok], exclude), collapse = ", ")
}

#' Check that a learner can be used
#'
#' Classed errors for a missing backend package (`"dependency"`) and for a
#' target type or distribution the learner does not support (`"model"`); a
#' classed warning when `num_scale` is given for a learner that has no scale.
#' @noRd
.learner_check <- function(learner, is_categorical, num_dist, num_scale = "auto") {
  for (pkg in learner$requires) {
    if (!.rapid_has_package(pkg)) {
      .rapid_abort(c("Package '", pkg, "' is required for model_type = \"", learner$name,
                     "\". Install it with install.packages(\"", pkg, "\")."), "dependency")
    }
  }
  type <- if (is_categorical) "categorical" else "numeric"
  if (!type %in% learner$targets) {
    .rapid_abort(c("model_type = \"", learner$name, "\" does not support a ", type,
                   " target. Models that do: ", .rapid_learners_for(type), "."), "model")
  }
  if (!is_categorical && !num_dist %in% c(learner$dists, "none")) {
    .rapid_abort(c("num_dist = \"", num_dist, "\" is not available for model_type = \"",
                   learner$name, "\". Available: ",
                   paste0("\"", learner$dists, "\"", collapse = ", "), "."), "model")
  }
  if (!identical(num_scale, "auto") && !"num_scale" %in% learner$options) {
    .rapid_warn(c("Argument 'num_scale' is ignored: it has no meaning for model_type = \"",
                  learner$name, "\"."), "ignored_argument")
  }
  invisible(TRUE)
}

# ---- Generics --------------------------------------------------------------------

#' Fit the attacker model on the synthetic data
#' @param num_dist `"auto"`, `"grid"`, `"exact"`, or `"none"` when no predictive
#'   distribution will be asked for.
#' @param num_scale,metric Working scale and error metric (used by lm).
#' @noRd
.learner_fit <- function(learner, formula, data, target_var, num_dist = "auto",
                         num_scale = "auto", metric = "symmetric") {
  UseMethod(".learner_fit")
}

#' Predict: probability matrix (categorical) or point predictions (numeric)
#' @noRd
.learner_predict <- function(learner, fit, newdata, target_var) UseMethod(".learner_predict")

#' Predictive distribution of a numeric target, one entry per record
#' @noRd
.learner_pdist <- function(learner, fit, newdata, train, target_var, num_dist = "auto",
                           n_quantiles = 100) {
  UseMethod(".learner_pdist")
}

#' Variable importance: named numeric vector or NULL
#' @noRd
.learner_importance <- function(learner, fit) UseMethod(".learner_importance")

#' One line with the settings of the fitted attacker that drive the result
#' @noRd
.learner_describe <- function(learner, fit) UseMethod(".learner_describe")

#' @exportS3Method .learner_pdist rapid_learner
#' @noRd
.learner_pdist.rapid_learner <- function(learner, fit, newdata, train, target_var,
                                         num_dist = "auto", n_quantiles = 100) {
  .rapid_abort(c("model_type = \"", learner$name, "\" has no predictive distribution ",
                 "for a numeric target."), "model")
}

#' @exportS3Method .learner_importance rapid_learner
#' @noRd
.learner_importance.rapid_learner <- function(learner, fit) NULL

#' @exportS3Method .learner_describe rapid_learner
#' @noRd
.learner_describe.rapid_learner <- function(learner, fit) learner$name

# ---- Shared helpers ----------------------------------------------------------------

# Midpoint probabilities of the grid: 0.5 %, 1.5 %, ..., 99.5 % for 100 quantiles
.pdist_probs <- function(n_quantiles) {
  if (!is.numeric(n_quantiles) || length(n_quantiles) != 1 || is.na(n_quantiles) ||
      n_quantiles < 2) {
    .rapid_abort("'n_quantiles' must be a single number of at least 2.", "input")
  }
  (seq_len(n_quantiles) - 0.5) / n_quantiles
}

# Stores the class levels of the target with the fit
.learner_tag <- function(fit, y) {
  attr(fit, "y_levels") <- if (is.factor(y)) levels(y) else NULL
  fit
}

# Classes without a training record are dropped for the fit (ranger warns about
# them, rpart cannot predict with them); `.learner_probs()` gives them
# probability 0 again.
.learner_drop_unused <- function(data, target_var) {
  if (is.factor(data[[target_var]])) data[[target_var]] <- droplevels(data[[target_var]])
  data
}

# Probability matrix with one column per class level, in the order of the levels
.learner_probs <- function(P, levels) {
  P <- as.matrix(P)
  if (is.null(colnames(P)) && ncol(P) == length(levels)) colnames(P) <- levels
  out <- matrix(0, nrow = nrow(P), ncol = length(levels), dimnames = list(NULL, levels))
  known <- intersect(levels, colnames(P))
  out[, known] <- P[, known]
  out
}

.lm_importance <- function(fit) {
  coefs <- stats::coef(fit)[-1]
  abs(coefs) / sum(abs(coefs))
}

# ---- Random forest (ranger) ----------------------------------------------------------

.ranger_predict <- function(fit, newdata) {
  p <- stats::predict(fit, data = newdata)$predictions
  levels <- attr(fit, "y_levels")
  if (is.null(levels)) as.numeric(p) else .learner_probs(p, levels)
}

# Terminal nodes of the records, one column per tree
.ranger_leaves <- function(fit, data) {
  stats::predict(fit, data = data, type = "terminalNodes")$predictions
}

# Exact predictive distribution of a forest: the training values that share a
# terminal node with the record, averaged over the trees
.forest_pdist <- function(fit, newdata, train, target_var) {
  .pdist_forest(.ranger_leaves(fit, newdata), .ranger_leaves(fit, train),
                train[[target_var]])
}

#' @exportS3Method .learner_fit rapid_learner_rf
#' @noRd
.learner_fit.rapid_learner_rf <- function(learner, formula, data, target_var,
                                          num_dist = "auto", num_scale = "auto",
                                          metric = "symmetric") {
  y <- data[[target_var]]
  data <- .learner_drop_unused(data, target_var)
  defaults <- list(formula = formula, data = data, probability = is.factor(y),
                   importance = "impurity")
  # The grid of a numeric target comes from a quantile regression forest
  if (!is.factor(y) && num_dist %in% c("auto", "grid")) defaults$quantreg <- TRUE
  fit <- do.call(ranger::ranger, utils::modifyList(defaults, learner$args))
  .learner_tag(fit, y)
}

#' @exportS3Method .learner_predict rapid_learner_rf
#' @noRd
.learner_predict.rapid_learner_rf <- function(learner, fit, newdata, target_var) {
  .ranger_predict(fit, newdata)
}

#' @exportS3Method .learner_pdist rapid_learner_rf
#' @noRd
.learner_pdist.rapid_learner_rf <- function(learner, fit, newdata, train, target_var,
                                            num_dist = "auto", n_quantiles = 100) {
  if (num_dist == "exact") return(.forest_pdist(fit, newdata, train, target_var))
  probs <- .pdist_probs(n_quantiles)
  if (is.null(fit$random.node.values)) {
    .rapid_abort(c("The forest was fitted without quantreg = TRUE, so it has no ",
                   "predictive distribution."), "model")
  }
  q <- stats::predict(fit, data = newdata, type = "quantiles",
                      quantiles = probs)$predictions
  .pdist_grid(unname(as.matrix(q)))
}

#' @exportS3Method .learner_importance rapid_learner_rf
#' @noRd
.learner_importance.rapid_learner_rf <- function(learner, fit) {
  tryCatch(ranger::importance(fit), error = function(e) NULL)
}

#' @exportS3Method .learner_describe rapid_learner_rf
#' @noRd
.learner_describe.rapid_learner_rf <- function(learner, fit) {
  paste0("rf: ", fit$num.trees, " trees, min.node.size ", fit$min.node.size)
}

# ---- Linear model ----------------------------------------------------------------------

#' @exportS3Method .learner_fit rapid_learner_lm
#' @noRd
.learner_fit.rapid_learner_lm <- function(learner, formula, data, target_var,
                                          num_dist = "auto", num_scale = "auto",
                                          metric = "symmetric") {
  y <- data[[target_var]]
  fit_on <- function(log_scale) {
    d <- data
    if (log_scale) d[[target_var]] <- log(y)
    do.call(stats::lm, utils::modifyList(list(formula = formula, data = d), learner$args))
  }
  n_nonpositive <- sum(y <= 0, na.rm = TRUE)

  if (num_scale == "log") {
    if (n_nonpositive > 0) {
      .rapid_abort(c("num_scale = \"log\" needs a positive target, but ", n_nonpositive,
                     " synthetic value(s) are zero or negative. Use num_scale = ",
                     "\"identity\" or \"auto\"."), "model")
    }
    fit <- fit_on(TRUE)
    scale <- "log"
  } else {
    fit <- fit_on(FALSE)
    scale <- "identity"
    if (num_scale == "auto" && metric %in% c("symmetric", "stabilised_relative") &&
        n_nonpositive == 0) {
      # Compare the Gaussian likelihoods on the original scale: the log model
      # pays the Jacobian sum(log(y)).
      fit_log <- fit_on(TRUE)
      jacobian <- sum(stats::model.response(stats::model.frame(fit_log)))
      if (as.numeric(stats::logLik(fit_log)) - jacobian > as.numeric(stats::logLik(fit))) {
        fit <- fit_log
        scale <- "log"
      }
    }
  }
  attr(fit, "rapid_scale") <- scale
  attr(fit, "rapid_scale_auto") <- num_scale == "auto"
  .learner_tag(fit, y)
}

#' @exportS3Method .learner_predict rapid_learner_lm
#' @noRd
.learner_predict.rapid_learner_lm <- function(learner, fit, newdata, target_var) {
  p <- as.numeric(stats::predict(fit, newdata = newdata))
  # On the log scale the point prediction is the predictive median
  if (identical(attr(fit, "rapid_scale"), "log")) exp(p) else p
}

#' @exportS3Method .learner_pdist rapid_learner_lm
#' @noRd
.learner_pdist.rapid_learner_lm <- function(learner, fit, newdata, train, target_var,
                                            num_dist = "auto", n_quantiles = 100) {
  p <- stats::predict(fit, newdata = newdata, se.fit = TRUE)
  mean <- as.numeric(p$fit)
  sd <- as.numeric(sqrt(p$se.fit^2 + summary(fit)$sigma^2))
  log_scale <- identical(attr(fit, "rapid_scale"), "log")
  if (num_dist == "grid") {
    probs <- .pdist_probs(n_quantiles)
    D <- outer(mean, rep(1, n_quantiles)) + outer(sd, stats::qnorm(probs))
    return(.pdist_grid(unname(if (log_scale) exp(D) else D)))
  }
  if (log_scale) .pdist_lognormal(mean, sd) else .pdist_normal(mean, sd)
}

#' @exportS3Method .learner_importance rapid_learner_lm
#' @noRd
.learner_importance.rapid_learner_lm <- function(learner, fit) .lm_importance(fit)

#' @exportS3Method .learner_describe rapid_learner_lm
#' @noRd
.learner_describe.rapid_learner_lm <- function(learner, fit) {
  paste0("lm: ", attr(fit, "rapid_scale"), " scale",
         if (isTRUE(attr(fit, "rapid_scale_auto"))) " (chosen automatically)")
}

# ---- CART (rpart) ------------------------------------------------------------------------

#' @exportS3Method .learner_fit rapid_learner_cart
#' @noRd
.learner_fit.rapid_learner_cart <- function(learner, formula, data, target_var,
                                            num_dist = "auto", num_scale = "auto",
                                            metric = "symmetric") {
  y <- data[[target_var]]
  data <- .learner_drop_unused(data, target_var)
  defaults <- list(formula = formula, data = data,
                   method = if (is.factor(y)) "class" else "anova")
  # Grown large and pruned by cross-validation; a control object of the user
  # replaces these three values altogether.
  if (is.null(learner$args$control)) {
    defaults <- c(defaults, list(cp = 0, minbucket = 5, xval = 10))
  }
  fit <- do.call(rpart::rpart, utils::modifyList(defaults, learner$args))

  pruned <- isTRUE(learner$prune)
  if (pruned) {
    if (fit$control$xval < 1) {
      .rapid_abort(c("prune = TRUE needs cross-validation, but xval = 0. Set xval to a ",
                     "positive number or use prune = FALSE."), "model")
    }
    tab <- fit$cptable
    # A tree without a split has nothing to prune
    if (nrow(tab) > 1 && "xerror" %in% colnames(tab) && any(is.finite(tab[, "xerror"]))) {
      fit <- rpart::prune(fit, cp = tab[which.min(tab[, "xerror"]), "CP"])
    }
  }
  attr(fit, "rapid_pruned") <- pruned
  .learner_tag(fit, y)
}

#' @exportS3Method .learner_predict rapid_learner_cart
#' @noRd
.learner_predict.rapid_learner_cart <- function(learner, fit, newdata, target_var) {
  levels <- attr(fit, "y_levels")
  if (is.null(levels)) {
    as.numeric(stats::predict(fit, newdata = newdata))
  } else {
    .learner_probs(stats::predict(fit, newdata = newdata, type = "prob"), levels)
  }
}

#' @exportS3Method .learner_pdist rapid_learner_cart
#' @noRd
.learner_pdist.rapid_learner_cart <- function(learner, fit, newdata, train, target_var,
                                              num_dist = "auto", n_quantiles = 100) {
  y_train <- if (!is.null(fit$y)) fit$y else train[[target_var]]
  if (length(y_train) != length(fit$where)) {
    .rapid_abort("The training values of the tree are not available.", "model")
  }
  # Predicting with the row numbers of the frame as fitted values returns the
  # leaf of each record.
  leaf_fit <- fit
  leaf_fit$frame$yval <- seq_len(nrow(fit$frame))
  leaf_new <- as.integer(stats::predict(leaf_fit, newdata = newdata))
  values <- lapply(split(y_train, fit$where), sort)
  leaf <- match(leaf_new, as.integer(names(values)))

  if (num_dist == "grid") {
    probs <- .pdist_probs(n_quantiles)
    # Type-1 quantiles keep point masses such as zeros
    leaf_q <- vapply(values, function(v) stats::quantile(v, probs, type = 1, names = FALSE),
                     numeric(n_quantiles))
    return(.pdist_grid(unname(t(leaf_q[, leaf, drop = FALSE]))))
  }
  .pdist_leaf(leaf, unname(values))
}

#' @exportS3Method .learner_importance rapid_learner_cart
#' @noRd
.learner_importance.rapid_learner_cart <- function(learner, fit) fit$variable.importance

#' @exportS3Method .learner_describe rapid_learner_cart
#' @noRd
.learner_describe.rapid_learner_cart <- function(learner, fit) {
  n_leaves <- sum(fit$frame$var == "<leaf>")
  paste0("cart: ", n_leaves, if (n_leaves == 1) " leaf" else " leaves",
         if (isTRUE(attr(fit, "rapid_pruned"))) ", pruned by cross-validation" else
           ", not pruned")
}

# ---- Gradient boosting (xgboost) -----------------------------------------------------------

# Design matrix of the quasi-identifiers (no response needed in `data`)
.gbm_matrix <- function(formula, data) {
  tt <- stats::delete.response(stats::terms(formula, data = data))
  stats::model.matrix(tt, data = data)[, -1, drop = FALSE]
}

# User arguments are booster parameters, except those of xgb.train() itself
.gbm_split_args <- function(args) {
  train_names <- c("nrounds", "early_stopping_rounds", "evals", "print_every_n",
                   "maximize", "callbacks", "verbose")
  is_train <- names(args) %in% train_names
  list(params = args[!is_train],
       train = utils::modifyList(list(nrounds = 100, verbose = 0), args[is_train]))
}

#' @exportS3Method .learner_fit rapid_learner_gbm
#' @noRd
.learner_fit.rapid_learner_gbm <- function(learner, formula, data, target_var,
                                           num_dist = "auto", num_scale = "auto",
                                           metric = "symmetric") {
  y <- data[[target_var]]
  args <- .gbm_split_args(learner$args)
  if (is.factor(y)) {
    K <- nlevels(y)
    if (K < 2) {
      .rapid_abort("model_type = \"gbm\" needs a target with at least two classes.", "model")
    }
    label <- as.integer(y) - 1
    base <- if (K == 2) {
      list(objective = "binary:logistic", eval_metric = "logloss", verbosity = 0)
    } else {
      list(objective = "multi:softprob", num_class = K, eval_metric = "mlogloss",
           verbosity = 0)
    }
  } else {
    label <- y
    base <- list(objective = "reg:squarederror", eval_metric = "rmse", verbosity = 0)
  }
  dtrain <- xgboost::xgb.DMatrix(data = .gbm_matrix(formula, data), label = label)
  fit <- do.call(xgboost::xgb.train,
                 c(list(params = utils::modifyList(base, args$params), data = dtrain),
                   args$train))
  # Attributes, not `$<-`: assigning into the booster fails with xgboost 3
  attr(fit, "x_formula") <- formula
  attr(fit, "rapid_nrounds") <- args$train$nrounds
  .learner_tag(fit, y)
}

#' @exportS3Method .learner_predict rapid_learner_gbm
#' @noRd
.learner_predict.rapid_learner_gbm <- function(learner, fit, newdata, target_var) {
  X <- .gbm_matrix(attr(fit, "x_formula"), newdata)
  p <- stats::predict(fit, xgboost::xgb.DMatrix(X))
  levels <- attr(fit, "y_levels")
  if (is.null(levels)) return(as.numeric(p))
  K <- length(levels)
  if (K == 2) {
    p <- as.numeric(p)
    P <- cbind(1 - p, p)
  } else {
    # xgboost 3 returns an n x K matrix; older versions a vector by record
    P <- if (is.matrix(p)) p else matrix(p, ncol = K, byrow = TRUE)
  }
  dimnames(P) <- list(NULL, levels)
  P
}

#' @exportS3Method .learner_pdist rapid_learner_gbm
#' @noRd
.learner_pdist.rapid_learner_gbm <- function(learner, fit, newdata, train, target_var,
                                             num_dist = "auto", n_quantiles = 100) {
  probs <- .pdist_probs(n_quantiles)
  formula <- attr(fit, "x_formula")
  args <- .gbm_split_args(learner$args)
  # A second model with the quantile loss, same settings and rounds as the
  # point model
  params <- utils::modifyList(args$params, list(objective = "reg:quantileerror",
                                                quantile_alpha = probs, verbosity = 0))
  dtrain <- xgboost::xgb.DMatrix(data = .gbm_matrix(formula, train),
                                 label = train[[target_var]])
  q_fit <- xgboost::xgb.train(params = params, data = dtrain,
                              nrounds = args$train$nrounds, verbose = 0)
  q <- stats::predict(q_fit, xgboost::xgb.DMatrix(.gbm_matrix(formula, newdata)))
  q <- if (is.matrix(q)) q else matrix(q, ncol = n_quantiles, byrow = TRUE)
  # Crossing quantiles are sorted per record
  .pdist_grid(unname(t(apply(q, 1, sort))))
}

#' @exportS3Method .learner_importance rapid_learner_gbm
#' @noRd
.learner_importance.rapid_learner_gbm <- function(learner, fit) {
  imp <- tryCatch(xgboost::xgb.importance(model = fit), error = function(e) NULL)
  if (is.null(imp) || nrow(imp) == 0) return(NULL)
  stats::setNames(imp$Gain, imp$Feature)
}

#' @exportS3Method .learner_describe rapid_learner_gbm
#' @noRd
.learner_describe.rapid_learner_gbm <- function(learner, fit) {
  paste0("gbm: ", attr(fit, "rapid_nrounds"), " rounds")
}

# ---- Logistic regression ---------------------------------------------------------------------

#' @exportS3Method .learner_fit rapid_learner_logit
#' @noRd
.learner_fit.rapid_learner_logit <- function(learner, formula, data, target_var,
                                             num_dist = "auto", num_scale = "auto",
                                             metric = "symmetric") {
  y <- data[[target_var]]
  if (!is.factor(y) || nlevels(y) != 2) {
    .rapid_abort(c("model_type = \"logit\" needs a target with exactly two classes; '",
                   target_var, "' has ", if (is.factor(y)) nlevels(y) else "no",
                   " classes. Models for any number of classes: ",
                   .rapid_learners_for("categorical", exclude = "logit"), "."), "model")
  }
  fit <- do.call(stats::glm, utils::modifyList(
    list(formula = formula, data = data, family = stats::binomial()), learner$args))
  .learner_tag(fit, y)
}

#' @exportS3Method .learner_predict rapid_learner_logit
#' @noRd
.learner_predict.rapid_learner_logit <- function(learner, fit, newdata, target_var) {
  p <- as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
  P <- cbind(1 - p, p)
  dimnames(P) <- list(NULL, attr(fit, "y_levels"))
  P
}

#' @exportS3Method .learner_importance rapid_learner_logit
#' @noRd
.learner_importance.rapid_learner_logit <- function(learner, fit) .lm_importance(fit)

# ---- Exact predictive distributions ------------------------------------------------------------

#' Normal predictive distribution (linear model)
#' @noRd
.pdist_normal <- function(mean, sd) {
  structure(list(mean = mean, sd = sd, n = length(mean)),
            class = c("rapid_pdist_normal", "rapid_pdist"))
}

#' @exportS3Method .pdist_mass rapid_pdist_normal
#' @noRd
.pdist_mass.rapid_pdist_normal <- function(pd, lo, hi) {
  stats::pnorm(hi, pd$mean, pd$sd) - stats::pnorm(lo, pd$mean, pd$sd)
}

#' Lognormal predictive distribution (linear model on the log scale)
#' @noRd
.pdist_lognormal <- function(meanlog, sdlog) {
  structure(list(meanlog = meanlog, sdlog = sdlog, n = length(meanlog)),
            class = c("rapid_pdist_lognormal", "rapid_pdist"))
}

#' @exportS3Method .pdist_mass rapid_pdist_lognormal
#' @noRd
.pdist_mass.rapid_pdist_lognormal <- function(pd, lo, hi) {
  # No mass at or below zero
  cdf <- function(v) {
    v <- rep_len(v, pd$n)
    out <- numeric(pd$n)
    pos <- v > 0
    out[pos] <- stats::pnorm((log(v[pos]) - pd$meanlog[pos]) / pd$sdlog[pos])
    out
  }
  cdf(hi) - cdf(lo)
}

#' Predictive distribution of a tree: the training values in the leaf
#' @param leaf Index into `values`, one per record.
#' @param values List of sorted training values, one element per leaf.
#' @noRd
.pdist_leaf <- function(leaf, values) {
  structure(list(leaf = leaf, values = values, n = length(leaf)),
            class = c("rapid_pdist_leaf", "rapid_pdist"))
}

#' @exportS3Method .pdist_mass rapid_pdist_leaf
#' @noRd
.pdist_mass.rapid_pdist_leaf <- function(pd, lo, hi) {
  lo <- rep_len(lo, pd$n)
  hi <- rep_len(hi, pd$n)
  out <- numeric(pd$n)
  for (k in unique(pd$leaf)) {
    i <- which(pd$leaf == k)
    v <- pd$values[[k]]
    out[i] <- (findInterval(hi[i], v, left.open = TRUE) - findInterval(lo[i], v)) / length(v)
  }
  out
}

#' Predictive distribution of a forest: Meinshausen weights
#'
#' Per tree the record has the training values of its terminal node, each with
#' weight 1 over the node size; the mass of a band is the mean over the trees
#' of the weight inside it. All training values in the node count, in-bag or not.
#' @param leaf_new,leaf_train Terminal node ids, `n x T` and `m x T`.
#' @param y_train The `m` training values.
#' @noRd
.pdist_forest <- function(leaf_new, leaf_train, y_train) {
  leaf_new <- as.matrix(leaf_new)
  leaf_train <- as.matrix(leaf_train)
  ys <- sort(unique(y_train))
  width <- length(ys) + 1
  rank <- match(y_train, ys)
  # One sorted composite key (node, rank of the value) per tree
  keys <- lapply(seq_len(ncol(leaf_train)), function(t) sort(leaf_train[, t] * width + rank))
  structure(list(leaf = leaf_new, keys = keys, ys = ys, width = width, n = nrow(leaf_new)),
            class = c("rapid_pdist_forest", "rapid_pdist"))
}

#' @exportS3Method .pdist_mass rapid_pdist_forest
#' @noRd
.pdist_mass.rapid_pdist_forest <- function(pd, lo, hi) {
  lo <- rep_len(lo, pd$n)
  hi <- rep_len(hi, pd$n)
  # Values below hi and values up to lo: the band is open
  k_hi <- findInterval(hi, pd$ys, left.open = TRUE)
  k_lo <- findInterval(lo, pd$ys)
  total <- numeric(pd$n)
  for (t in seq_along(pd$keys)) {
    keys <- pd$keys[[t]]
    base <- pd$leaf[, t] * pd$width
    below <- findInterval(base, keys)
    size <- findInterval(base + pd$width, keys) - below
    inside <- findInterval(base + k_hi, keys) - findInterval(base + k_lo, keys)
    total <- total + ifelse(size > 0, pmax(inside, 0) / size, 0)
  }
  total / length(pd$keys)
}

# ---- Lookup attacker (one fully grown tree) -------------------------------------------------

#' @exportS3Method .learner_fit rapid_learner_lookup
#' @noRd
.learner_fit.rapid_learner_lookup <- function(learner, formula, data, target_var,
                                              num_dist = "auto", num_scale = "auto",
                                              metric = "symmetric") {
  y <- data[[target_var]]
  data <- .learner_drop_unused(data, target_var)
  n_keys <- length(attr(stats::terms(formula, data = data), "term.labels"))
  defaults <- list(formula = formula, data = data, probability = is.factor(y),
                   num.trees = 1, replace = FALSE, sample.fraction = 1, min.node.size = 1,
                   mtry = n_keys)
  fit <- do.call(ranger::ranger, utils::modifyList(defaults, learner$args))
  .learner_tag(fit, y)
}

#' @exportS3Method .learner_predict rapid_learner_lookup
#' @noRd
.learner_predict.rapid_learner_lookup <- function(learner, fit, newdata, target_var) {
  .ranger_predict(fit, newdata)
}

#' @exportS3Method .learner_pdist rapid_learner_lookup
#' @noRd
.learner_pdist.rapid_learner_lookup <- function(learner, fit, newdata, train, target_var,
                                                num_dist = "auto", n_quantiles = 100) {
  .forest_pdist(fit, newdata, train, target_var)
}

#' @exportS3Method .learner_describe rapid_learner_lookup
#' @noRd
.learner_describe.rapid_learner_lookup <- function(learner, fit) {
  n_leaves <- sum(ranger::treeInfo(fit, tree = 1)$terminal)
  paste0("lookup: 1 fully grown tree, ", n_leaves, " leaves")
}

# ---- User-defined learners -------------------------------------------------------------------

#' User-Defined Attacker Model for RAPID
#'
#' Builds an attacker model from functions of the user, for the argument
#' \code{model_type} of \code{\link{rapid}}, \code{\link{rapid_suite}} and the
#' functions that fit through \code{rapid()} (\code{\link{rapid_test}},
#' \code{\link{rapid_threshold_select}}, \code{\link{rapid_qi_contribution}},
#' \code{\link{rapid_synthesizer_cv}}).
#'
#' The functions are called with their arguments by position, in the order
#' given in the argument descriptions below. An argument name of that list
#' that a function has at another position is an error of class
#' \code{"rapid_error_input"}.
#'
#' \code{fit(formula, data, ...)} trains the attacker on the synthetic data:
#' \code{formula} is \code{target ~ key_var1 + key_var2 + ...}, \code{data}
#' holds the prepared synthetic records (classes of a categorical target that no
#' synthetic record has are dropped from its levels), and \code{...} are the
#' arguments given in \code{...} of \code{rapid()}.
#'
#' \code{predict_prob(fit, newdata)} returns, for a categorical target, a
#' matrix with one row per record of \code{newdata} and the class names as
#' column names. A class missing from the columns gets probability 0.
#' \code{predict_point(fit, newdata)} returns, for a numeric target, one point
#' prediction per record. \code{predict_quantiles(fit, newdata, probs)} returns
#' the predictive distribution as an \code{nrow(newdata) x length(probs)}
#' matrix of quantiles. The band rule (\code{num_eval_method = "band"}, the
#' default) needs it, the point rule does not. \code{num_dist = "auto"} and
#' \code{"grid"} both use these quantiles, and \code{"exact"} is not available.
#'
#' The output of the functions is checked (dimensions, names, no \code{NA}),
#' and a violation is an error of class \code{"rapid_error_model"}. Errors raised
#' inside the functions are passed on unchanged. Records with a
#' quasi-identifier level that the synthetic data do not contain are refused
#' with an error of class \code{"rapid_error_levels"}, before the functions are
#' called. Variable importance (\code{plot(x, which = 4)}) is not available.
#'
#' The functions are stored with the resulting \code{rapid} object as given,
#' including the environment of a closure, so define them at top level or keep
#' them free of large captured objects.
#'
#' @param name character, the name of the model; it is the \code{model_type}
#'   of the \code{rapid} object and the name of the model in a suite
#' @param fit function \code{(formula, data)} returning the fitted model
#' @param predict_prob function \code{(fit, newdata)} for a categorical target,
#'   or \code{NULL}
#' @param predict_point function \code{(fit, newdata)} for a numeric target, or
#'   \code{NULL}
#' @param predict_quantiles function \code{(fit, newdata, probs)} giving the
#'   predictive distribution of a numeric target, or \code{NULL}
#' @param requires character vector of packages the functions need; a missing
#'   one is an error of class \code{"rapid_error_dependency"} (in
#'   \code{rapid_suite()} the model is skipped)
#'
#' @return \code{rapid_learner()}: an object of class
#'   \code{c("rapid_learner_custom", "rapid_learner")}. Its field
#'   \code{targets} is \code{"categorical"} when \code{predict_prob} is given
#'   and \code{"numeric"} when \code{predict_point} is given. \code{format()}
#'   returns the lines that \code{print()} shows (name, target types, required
#'   packages and the functions given), \code{print()} the object, invisibly.
#'
#' @seealso \code{\link{rapid}}, \code{\link{rapid_suite}}
#' @family rapid
#' @examples
#' set.seed(1)
#' X <- data.frame(age = sample(20:60, 100, replace = TRUE),
#'                 hours = sample(10:50, 100, replace = TRUE))
#' X$income <- 20000 + 600 * X$age + 300 * X$hours + rnorm(100, 0, 5000)
#' Y <- X
#' Y$income <- Y$income + rnorm(100, 0, 3000)
#'
#' # A linear model with its normal predictive distribution
#' my_lm <- rapid_learner(
#'   "my_lm",
#'   fit = stats::lm,
#'   predict_point = function(fit, newdata) predict(fit, newdata),
#'   predict_quantiles = function(fit, newdata, probs) {
#'     p <- predict(fit, newdata, se.fit = TRUE)
#'     sd <- sqrt(p$se.fit^2 + p$residual.scale^2)
#'     outer(p$fit, rep(1, length(probs))) + outer(sd, qnorm(probs))
#'   }
#' )
#' my_lm
#' format(my_lm)
#' rapid(X, Y, key_vars = c("age", "hours"), target_var = "income",
#'       model_type = my_lm, n_perm = 3, seed = 1)
#'
#' # Without predict_quantiles only the point rule is available
#' my_point <- rapid_learner("my_point", fit = stats::lm,
#'                           predict_point = function(fit, newdata) predict(fit, newdata))
#' rapid(X, Y, key_vars = c("age", "hours"), target_var = "income",
#'       model_type = my_point, num_eval_method = "point", n_perm = 3, seed = 1)
#'
#' # A categorical target: logistic regression with class probabilities
#' X$high <- factor(ifelse(X$income > median(X$income), "yes", "no"))
#' Y$high <- X$high
#' my_logit <- rapid_learner(
#'   "my_logit",
#'   fit = function(formula, data) stats::glm(formula, data = data, family = binomial()),
#'   predict_prob = function(fit, newdata) {
#'     p <- predict(fit, newdata, type = "response")
#'     cbind(no = 1 - p, yes = p)
#'   }
#' )
#' rapid(X, Y, key_vars = c("age", "hours"), target_var = "high",
#'       model_type = my_logit, n_perm = 3, seed = 1)
#' @export
rapid_learner <- function(name, fit, predict_prob = NULL, predict_point = NULL,
                          predict_quantiles = NULL, requires = character()) {
  if (!is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name)) {
    .rapid_abort("'name' must be a single non-empty character string.", "input")
  }
  .learner_check_fun(fit, "fit", c("formula", "data"))
  .learner_check_fun(predict_prob, "predict_prob", c("fit", "newdata"), optional = TRUE)
  .learner_check_fun(predict_point, "predict_point", c("fit", "newdata"), optional = TRUE)
  .learner_check_fun(predict_quantiles, "predict_quantiles", c("fit", "newdata", "probs"),
                     optional = TRUE)
  if (is.null(predict_prob) && is.null(predict_point)) {
    .rapid_abort(c("Give 'predict_prob' (categorical target), 'predict_point' (numeric ",
                   "target), or both."), "input")
  }
  if (!is.null(predict_quantiles) && is.null(predict_point)) {
    .rapid_abort(c("'predict_quantiles' belongs to a numeric target and needs ",
                   "'predict_point' as well."), "input")
  }
  if (!is.character(requires) || anyNA(requires)) {
    .rapid_abort("'requires' must be a character vector of package names.", "input")
  }
  structure(
    list(name = name, args = list(), requires = requires,
         targets = c(if (!is.null(predict_point)) "numeric",
                     if (!is.null(predict_prob)) "categorical"),
         dists = c("auto", "grid"), options = character(0), unseen_levels = FALSE,
         note = NULL, fit = fit, predict_prob = predict_prob, predict_point = predict_point,
         predict_quantiles = predict_quantiles),
    class = c("rapid_learner_custom", "rapid_learner")
  )
}

# A function of the user must have the named arguments at their positions, or `...`
.learner_check_fun <- function(f, what, args, optional = FALSE) {
  if (optional && is.null(f)) return(invisible(TRUE))
  if (!is.function(f)) {
    .rapid_abort(c("'", what, "' must be a function", if (optional) " or NULL", "."), "input")
  }
  formals_f <- names(formals(args(f)))
  missing_args <- setdiff(args, formals_f)
  if (length(missing_args) > 0 && !"..." %in% formals_f) {
    .rapid_abort(c("'", what, "' must be a function (", paste(args, collapse = ", "),
                   "); it has no argument ", paste0("'", missing_args, "'", collapse = ", "),
                   "."), "input")
  }
  # The call is by position: a required name it has must sit at its position
  position <- match(args, formals_f)
  misplaced <- args[!is.na(position) & position != seq_along(args)]
  if (length(misplaced) > 0) {
    .rapid_abort(c("'", what, "' is called by position as (", paste(args, collapse = ", "),
                   "); its argument ", paste0("'", misplaced, "'", collapse = ", "),
                   " is not at that position."), "input")
  }
  invisible(TRUE)
}

#' @rdname rapid_learner
#' @param x a \code{rapid_learner} object
#' @param ... ignored
#' @export
format.rapid_learner <- function(x, ...) {
  c(paste0("RAPID attacker model: ", x$name),
    paste0("  Targets:  ", paste(x$targets, collapse = ", ")),
    paste0("  Requires: ", if (length(x$requires) > 0) paste(x$requires, collapse = ", ")
                           else "no package"))
}

#' @rdname rapid_learner
#' @export
format.rapid_learner_custom <- function(x, ...) {
  given <- c("fit", "predict_prob", "predict_point", "predict_quantiles")
  given <- given[!vapply(x[given], is.null, logical(1))]
  c(NextMethod(), paste0("  User-defined: ", paste(given, collapse = ", ")))
}

#' @rdname rapid_learner
#' @export
print.rapid_learner <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

# Output of a user's function that breaks its contract
.learner_bad_output <- function(learner, what, problem) {
  .rapid_abort(c("'", what, "' of the learner \"", learner$name, "\" ", problem, "."),
               "model")
}

#' @exportS3Method .learner_fit rapid_learner_custom
#' @noRd
.learner_fit.rapid_learner_custom <- function(learner, formula, data, target_var,
                                              num_dist = "auto", num_scale = "auto",
                                              metric = "symmetric") {
  y <- data[[target_var]]
  data <- .learner_drop_unused(data, target_var)
  fit <- do.call(learner$fit, c(list(formula, data), learner$args))
  .learner_tag(fit, y)
}

#' @exportS3Method .learner_predict rapid_learner_custom
#' @noRd
.learner_predict.rapid_learner_custom <- function(learner, fit, newdata, target_var) {
  n <- nrow(newdata)
  levels <- attr(fit, "y_levels")
  if (is.null(levels)) {
    p <- learner$predict_point(fit, newdata)
    if (!is.numeric(p) || length(p) != n) {
      .learner_bad_output(learner, "predict_point",
                          paste0("must return one number per record of newdata (", n, ")"))
    }
    if (anyNA(p)) .learner_bad_output(learner, "predict_point", "returned NA")
    return(as.numeric(p))
  }
  P <- learner$predict_prob(fit, newdata)
  if (is.data.frame(P)) P <- as.matrix(P)
  if (!is.matrix(P) || !is.numeric(P) || nrow(P) != n) {
    .learner_bad_output(learner, "predict_prob",
                        paste0("must return a numeric matrix with one row per record ",
                               "of newdata (", n, ")"))
  }
  if (is.null(colnames(P))) {
    .learner_bad_output(learner, "predict_prob",
                        "must return a matrix with the class names as column names")
  }
  unknown <- setdiff(colnames(P), levels)
  if (length(unknown) > 0 || anyDuplicated(colnames(P)) > 0) {
    .learner_bad_output(learner, "predict_prob",
                        paste0("must have distinct class names of the target as column ",
                               "names; unknown: ", paste(unknown, collapse = ", ")))
  }
  if (anyNA(P)) .learner_bad_output(learner, "predict_prob", "returned NA")
  .learner_probs(P, levels)
}

#' @exportS3Method .learner_pdist rapid_learner_custom
#' @noRd
.learner_pdist.rapid_learner_custom <- function(learner, fit, newdata, train, target_var,
                                                num_dist = "auto", n_quantiles = 100) {
  if (is.null(learner$predict_quantiles)) {
    # Only a rule that asks for the distribution fails (the point rule does not)
    return(.pdist_unavailable(nrow(newdata), c(
      "The band rule needs a predictive distribution, but the learner \"", learner$name,
      "\" has no 'predict_quantiles'. Give it to rapid_learner() or use ",
      "num_eval_method = \"point\".")))
  }
  probs <- .pdist_probs(n_quantiles)
  Q <- learner$predict_quantiles(fit, newdata, probs)
  if (is.data.frame(Q)) Q <- as.matrix(Q)
  if (!is.matrix(Q) || !is.numeric(Q) || nrow(Q) != nrow(newdata) ||
      ncol(Q) != length(probs)) {
    .learner_bad_output(learner, "predict_quantiles",
                        paste0("must return a numeric matrix with one row per record (",
                               nrow(newdata), ") and one column per probability (",
                               length(probs), ")"))
  }
  if (anyNA(Q)) .learner_bad_output(learner, "predict_quantiles", "returned NA")
  .pdist_grid(unname(Q))
}

#' @exportS3Method .learner_describe rapid_learner_custom
#' @noRd
.learner_describe.rapid_learner_custom <- function(learner, fit) {
  paste0(learner$name, ": user-defined learner")
}

#' A predictive distribution that is not available
#'
#' Stands in for the distribution of a learner that has none, so that rules
#' without a distribution work; the mass raises `msg` as a `"model"` error.
#' @noRd
.pdist_unavailable <- function(n, msg) {
  structure(list(n = n, msg = msg), class = c("rapid_pdist_unavailable", "rapid_pdist"))
}

#' @exportS3Method .pdist_mass rapid_pdist_unavailable
#' @noRd
.pdist_mass.rapid_pdist_unavailable <- function(pd, lo, hi) .rapid_abort(pd$msg, "model")

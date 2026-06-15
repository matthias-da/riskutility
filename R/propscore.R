#' Propensity score utility measure
#'
#' How well can you predict if an observation is from which data set?
#'
#' @param X data frame
#' @param Y data frame with the same structure as X
#' @param form formula. If NULL all variables are used as predictors
#' @param method method for propensity score estimation: \code{"rf"} (default,
#' uses \code{randomForest}), \code{"ranger"} (uses \code{ranger} with
#' proximity-based structural metrics), or \code{"logreg"} (logistic
#' regression).
#' @param adjust_size for very unbalanced sizes of original and synthetic data.
#' Instead of using a constant c, observations are drawn from the smaller
#' data set so that the both data sets has the same size. If a cluster
#' structure is present (e.g. person in housholds), draws are taken from the
#' clusters (e.g. households).
#' @param cluster vector specifiying the cluster structure. Should be NULL
#' if no cluster structure is present in the data.
#' @param na missing value treatment. Either stop, remove or impute
#' (using a kNN from R package VIM).
#' @param proximity character; proximity computation mode for
#' \code{method = "ranger"}: \code{"summary"} (default) returns aggregate
#' within-class and cross-class proximity statistics, \code{"full"} also
#' stores the full proximity matrix, \code{"none"} skips proximity
#' computation. Ignored for other methods.
#' @param importance logical; whether to compute variable importance for
#' \code{method = "ranger"}. Default \code{TRUE}. Ignored for other methods.
#' @param ... additional arguments passed to methods. For
#' \code{method = "ranger"}, extra arguments are forwarded to
#' \code{\link[ranger]{ranger}} via \code{modifyList}.
#' @importFrom randomForest randomForest
#' @importFrom stats as.formula formula
#' @importFrom stats terms.formula
#' @importFrom stats complete.cases
#' @author Matthias Templ
#' @return An S3 object of class \code{"propscore"} containing:
#' \describe{
#'   \item{predictions}{Data frame of predicted propensity scores for all records.}
#'   \item{ps_ratio}{Propensity score ratio (mean propensity of synthetic / original).}
#'   \item{ps_score}{Propensity score statistic.}
#'   \item{cr}{Classification rate.}
#'   \item{mean_ps_x}{Mean propensity score for original records.}
#'   \item{mean_ps_y}{Mean propensity score for synthetic records.}
#'   \item{density_ratio}{Density ratio of propensity score distributions.}
#'   \item{density_ratio_bayes}{Bayesian density ratio.}
#'   \item{kl}{Kullback-Leibler divergence between propensity distributions.}
#'   \item{kl_bayes}{Bayesian KL divergence.}
#'   \item{mean_ratio}{Mean ratio of propensity densities.}
#'   \item{sd_ratio}{Standard deviation of density ratio.}
#'   \item{mean_ratio_bayes}{Mean Bayesian density ratio.}
#'   \item{sd_ratio_bayes}{Standard deviation of Bayesian density ratio.}
#'   \item{points}{Grid points used for density estimation.}
#'   \item{denX}{Density estimates for original data propensity scores.}
#'   \item{denY}{Density estimates for synthetic data propensity scores.}
#'   \item{bayesspace}{Logical, whether Bayesian space was used.}
#'   \item{n_x}{Number of original records.}
#'   \item{n_y}{Number of synthetic records.}
#'   \item{method}{Method used for propensity score estimation.}
#'   \item{oob_error}{OOB prediction error (ranger method only).}
#'   \item{var_importance}{Named numeric vector of variable importance (ranger method only).}
#'   \item{within_orig_prox}{Mean within-original proximity (ranger with proximity != "none").}
#'   \item{within_synth_prox}{Mean within-synthetic proximity (ranger with proximity != "none").}
#'   \item{cross_prox}{Mean cross-class proximity (ranger with proximity != "none").}
#'   \item{structure_ratio}{Ratio of cross-class to within-class proximity (ranger with proximity != "none").}
#'   \item{proximity_matrix}{Full proximity matrix (ranger with proximity = "full" only).}
#' }
#' @family utility
#' @export
#' @references
#' Templ, M. Statistical Disclosure Control for Microdata: Methods and Applications in R.
#' \emph{Springer International Publishing}, 287 pages, 2017. ISBN 978-3-319-50272-4. \doi{10.1007/978-3-319-50272-4}
#' @examples
#' # Simple example with synthetic data
#' set.seed(123)
#' X <- data.frame(
#'   age = sample(20:60, 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   income = rnorm(100, 50000, 10000)
#' )
#' Y <- data.frame(
#'   age = sample(20:60, 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   income = rnorm(100, 48000, 12000)
#' )
#' propscore(X, Y, form = ~ age + gender + income)
#'
#' \donttest{
#' # Extended example using simPop and sdcMicro (requires these packages)
#' if (requireNamespace("simPop", quietly = TRUE) &&
#'     requireNamespace("sdcMicro", quietly = TRUE)) {
#'   data(eusilc13puf, package = "simPop")
#'   eusilc13puf$age <- as.numeric(as.character(eusilc13puf$age))
#'   keyvars <- c("age", "rb090", "db040", "pl031", "pb220a")
#'   sdc <- sdcMicro::createSdcObj(eusilc13puf,
#'                                 keyVars = keyvars,
#'                                 numVars = "pgrossIncome",
#'                                 w = "rb050",
#'                                 hhId = "db030")
#'   sdc <- sdcMicro::localSuppression(sdc)
#'   sdc <- sdcMicro::microaggregation(sdc)
#'   eusilc13puf_anon <- sdcMicro::extractManipData(sdc)
#'   propscore(eusilc13puf, eusilc13puf_anon, na = "remove",
#'      form = ~ age + rb090 + pl031 + db040 + pb220a + pgrossIncome)
#' }
#' }
#'
#' # Account for a survey cluster structure (e.g. households) and keep the
#' # original/synthetic sizes as-is via adjust_size = FALSE
#' set.seed(123)
#' Xc <- data.frame(
#'   age = sample(20:60, 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   income = rnorm(100, 50000, 10000),
#'   hid = factor(sample(1:25, 100, replace = TRUE))
#' )
#' Yc <- data.frame(
#'   age = sample(20:60, 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   income = rnorm(100, 48000, 12000),
#'   hid = factor(sample(1:25, 100, replace = TRUE))
#' )
#' propscore(Xc, Yc, form = ~ age + gender + income,
#'           cluster = "hid", adjust_size = FALSE)
propscore <- function(X, ...) {
  UseMethod("propscore")
}

#' @rdname propscore
#' @export
propscore.synth_pair <- function(X, form = NULL, ...) {
  # Build formula from vars if not provided
  if (is.null(form) && length(X$vars) > 0) {
    form_str <- paste("~", paste(X$vars, collapse = " + "))
    form <- as.formula(form_str)
  }

  propscore.default(
    X = X$original,
    Y = X$synthetic,
    form = form,
    ...
  )
}

#' @rdname propscore
#' @export
propscore.default <- function(X, Y,
                              form = NULL,
                              method = c("rf", "ranger", "logreg"),
                              adjust_size = TRUE,
                              cluster = NULL,
                              na = "impute",
                              proximity = c("summary", "full", "none"),
                              importance = TRUE,
                              ...) {
  # Capture missing() before match.arg resolves defaults
  proximity_supplied <- !missing(proximity)
  importance_supplied <- !missing(importance)
  method <- match.arg(method)
  proximity <- match.arg(proximity)
  # Check if X and Y are data frames
  if (!is.data.frame(X)) {
    stop("X must be a data frame.")
  }
  if (!is.data.frame(Y)) {
    stop("Y must be a data frame.")
  }

  # Helper function to check and modify the formula
  check_and_modify_formula <- function(form) {
    form_str <- deparse(form)
    if (!grepl("^group ~", form_str)) {
      form_str <- paste("group ~", sub("^~", "", form_str))
      form <- as.formula(form_str)
    }
    return(form)
  }

  if (is.null(form)) {
    # If no formula is provided, check if X and Y have the same structure
    if (!all(names(X) == names(Y))) {
      stop("X and Y must have the same structure.")
    }
    if (!all(sapply(X, class) == sapply(Y, class))) {
      stop("X and Y must have the same column classes.")
    }
    #form <- paste("group", paste(colnames(X), collapse = " + "), collapse = " ~ ")
    form <- formula("group ~ .")
  } else {
    # Check if form is of class 'formula'
    if (!inherits(form, "formula")) {
      stop("form must be a valid formula.")
    }

    form <- check_and_modify_formula(form)

    # Extract variable names from the formula
    formula_vars <- all.vars(form)[-1]

    # Check if variables in the formula appear in X (and Y)
    missing_vars <- setdiff(formula_vars, names(X))
    if (length(missing_vars) > 0) {
      stop(paste("The following variables from the formula are missing in X (and Y):", paste(missing_vars, collapse = ", ")))
    }

    # Check if the columns in X and Y for the formula variables have the same structure
    for (var in formula_vars) {
      if (!identical(class(X[[var]]), class(Y[[var]]))) {
        stop(paste("The variable", var, "must have the same class in both X and Y."))
      }
    }
  }

  # Check the cluster argument
  if (!is.null(cluster)) {
    if (!is.character(cluster) || length(cluster) != 1) {
      stop("cluster must be a single string representing a variable name in X.")
    }
    if (!(cluster %in% names(X))) {
      stop(paste("The cluster variable", cluster, "is not found in X."))
    }
  }

  # Check the na argument
  valid_na_values <- c("impute", "remove", "stop")
  if (!(na %in% valid_na_values)) {
    stop(paste("na must be one of", paste(valid_na_values, collapse = ", "), "."))
  }

  # Extract variable names; handle "group ~ ." formula specially
  if (identical(deparse(form), "group ~ .")) {
    vars <- names(X)
  } else {
    vars <- attributes(terms.formula(form))$term.labels
  }
  if(any(is.na(X[, vars]))){
    if(na == "stop") stop("data contains missing values")
    if(na == "impute"){
      X[, vars] <- VIM::kNN(X[, vars], imp_var = FALSE)
      Y[, vars] <- VIM::kNN(Y[, vars], imp_var = FALSE)
    }
    if(na == "remove"){
      X <- X[complete.cases(X[, vars]), ]
      Y <- Y[complete.cases(Y[, vars]), ]
    }
  }

  # Function to sample data based on the cluster variable
  sample_by_cluster <- function(X, Y, cluster) {
    X_clusters <- unique(X[[cluster]])
    Y_clusters <- unique(Y[[cluster]])

    if (nrow(X) > nrow(Y)) {
      # Sample households from X to match the number of persons in Y
      sampled_clusters <- sample(X_clusters, length(Y_clusters))
      X_sampled <- X[X[[cluster]] %in% sampled_clusters, ]
      return(list(X = X_sampled, Y = Y))
    } else {
      # Sample households from Y to match the number of persons in X
      sampled_clusters <- sample(Y_clusters, length(X_clusters))
      Y_sampled <- Y[Y[[cluster]] %in% sampled_clusters, ]
      return(list(X = X, Y = Y_sampled))
    }
  }

  # Sample the data based on the cluster variable if cluster is provided
  if (!is.null(cluster)) {
    sampled_data <- sample_by_cluster(X, Y, cluster)
    X <- sampled_data$X
    Y <- sampled_data$Y
  }

  mi_x <- n_x <- nrow(X)
  mi_y <- n_y <- nrow(Y)

  cr <- n_y / (n_x + n_y)

  if(adjust_size & !(identical(n_x, n_y))){
#    ma <- max(c(n_x, n_y))
    if(n_x > n_y){
      index0 <- sample(1:n_x, size = n_y, replace = FALSE)
      X <- X[index0, ]
    } else {
      index0 <- sample(1:n_y, size = n_x, replace = FALSE)
      Y <- Y[index0, ]
    }
    cr <- nrow(Y) / (nrow(X) + nrow(Y))
    mi_x <- mi_y <- min(c(n_x, n_y))
  }

  # Warn if ranger-only params used with non-ranger method
  if (method != "ranger" && (proximity_supplied || importance_supplied)) {
    message("'proximity' and 'importance' are only used with method = 'ranger'.")
  }

  # --- ranger branch: uses .rf_proximity() engine, returns early ---
  if (method == "ranger") {
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("Package 'ranger' required for propscore(method = 'ranger'). ",
           "Install with install.packages('ranger')", call. = FALSE)
    }

    # Use .rf_proximity() engine
    rf_res <- .rf_proximity(X[, vars, drop = FALSE], Y[, vars, drop = FALSE],
                            vars = vars, n_trees = 500L,
                            importance = importance, ...)

    # OOB propensity scores
    p_hat <- rf_res$forest$predictions[, 2]  # P(synthetic)
    # Handle OOB prediction NAs
    if (any(is.na(p_hat))) {
      p_inbag <- predict(rf_res$forest,
                         rbind(X[, vars, drop = FALSE],
                               Y[, vars, drop = FALSE]))$predictions[, 2]
      p_hat[is.na(p_hat)] <- p_inbag[is.na(p_hat)]
    }

    n_total <- nrow(X) + nrow(Y)
    pmse <- (1 / n_total) * sum((p_hat - cr)^2)

    # Proximity-based structural metrics
    prox_fields <- list(
      within_orig_prox = NULL,
      within_synth_prox = NULL,
      cross_prox = NULL,
      structure_ratio = NULL,
      proximity_matrix = NULL
    )
    if (proximity != "none") {
      n1 <- rf_res$n1
      n2 <- rf_res$n2
      idx_orig <- seq_len(n1)
      idx_synth <- n1 + seq_len(n2)

      # Within-original proximity (mean of upper triangle)
      if (n1 > 1) {
        wo <- .proximity_from_nodes(rf_res$terminal_nodes, idx_orig, idx_orig)
        prox_fields$within_orig_prox <- (sum(wo) - n1) / (n1 * (n1 - 1))
      } else {
        prox_fields$within_orig_prox <- NA_real_
      }

      # Within-synthetic proximity
      if (n2 > 1) {
        ws <- .proximity_from_nodes(rf_res$terminal_nodes, idx_synth, idx_synth)
        prox_fields$within_synth_prox <- (sum(ws) - n2) / (n2 * (n2 - 1))
      } else {
        prox_fields$within_synth_prox <- NA_real_
      }

      # Cross-class proximity
      cross <- .proximity_from_nodes(rf_res$terminal_nodes, idx_orig, idx_synth)
      prox_fields$cross_prox <- mean(cross)

      # Structure ratio
      denom <- (prox_fields$within_orig_prox +
                prox_fields$within_synth_prox) / 2
      prox_fields$structure_ratio <- if (!is.na(denom) && denom > 0) {
        prox_fields$cross_prox / denom
      } else {
        NA_real_
      }

      if (proximity == "full") {
        n_combined <- n1 + n2
        if (n_combined > 10000) {
          mem_mb <- round(n_combined^2 * 8 / 1e6)
          warning("Full proximity matrix is ", n_combined, " x ",
                  n_combined, " (~", mem_mb, " MB). ",
                  "Consider proximity = 'summary'.", call. = FALSE)
        }
        all_idx <- seq_len(n_combined)
        prox_fields$proximity_matrix <- .proximity_from_nodes(
          rf_res$terminal_nodes, all_idx, all_idx
        )
      }
    }

    oob_error <- rf_res$oob_error
    if (oob_error > 0.45) {
      message("OOB error > 0.45: forest has little discriminative power.")
    }

    # Build result — ranger bypasses KDE and returns directly
    p <- data.frame(
      prediction = p_hat,
      group = rep(c("real", "synth"), times = c(nrow(X), nrow(Y)))
    )

    results <- list(
      predictions = p,
      ps_ratio = NA_real_,
      ps_score = pmse,
      cr = cr,
      mean_ps_x = mean(p_hat[seq_len(nrow(X))], na.rm = TRUE),
      mean_ps_y = mean(p_hat[nrow(X) + seq_len(nrow(Y))], na.rm = TRUE),
      density_ratio = NULL,
      density_ratio_bayes = NULL,
      kl = NA_real_,
      kl_bayes = NA_real_,
      mean_ratio = NA_real_,
      sd_ratio = NA_real_,
      mean_ratio_bayes = NA_real_,
      sd_ratio_bayes = NA_real_,
      points = NULL,
      denX = NULL,
      denY = NULL,
      bayesspace = FALSE,
      n_x = n_x,
      n_y = n_y,
      method = "ranger",
      oob_error = oob_error,
      var_importance = rf_res$importance
    )

    # Add proximity fields
    for (nm in names(prox_fields)) results[[nm]] <- prox_fields[[nm]]

    class(results) <- "propscore"
    return(results)
  }

    Z <- rbind(X[, vars], Y[, vars])
    Z$group <- factor(rep(c(0,1), times = c(mi_x, mi_y)))
    if(method == "logreg"){
      res <- glm(form, data = Z, family = binomial())
      p <- data.frame("prediction" = predict(res, type = "response"),
                      "group" = rep(c("real", "synth"), times = c(mi_x, mi_y)))
    } else{
      res <- randomForest::randomForest(form, data = Z)
      # p <- data.frame("prediction" = predict(res, data = Z, predict.all = TRUE)$pred,
      #                 "group" = rep(c("real", "synth"), times = c(mi_x, mi_y)))
      p <- data.frame("prediction" = res$votes[,2],
                      "group" = rep(c("real", "synth"), times = c(mi_x, mi_y)))
    }

    ps <- 1 / (2 * mi_x) * sum((p$prediction - cr)^2)
    t1 <- as.numeric(table(p$prediction < cr))
    ps_ratio <- t1[1] / t1[2]

    # Kernel Density Estimation, difference in distributions
    kde_x <- density(p$prediction[1:mi_x])
    kde_y <- density(p$prediction[(mi_x+1):length(ps)])
    # Define a sequence of points where we want to compute the density ratio
    points <- seq(min(kde_x$x, kde_y$x), max(kde_x$x, kde_y$x), length.out = 1000)
    # Interpolate densities at the sequence of points
    density_X <- density_X_orig <- approx(kde_x$x, kde_x$y, xout = points)$y
    density_Y <- density_Y_orig <- approx(kde_y$x, kde_y$y, xout = points)$y
    distance <- sqrt(sum((density_X - density_Y)^2, na.rm = TRUE))
    # Compute the density ratio
    density_ratio <- density_X / density_Y
    kl <- sum(density_X * log(density_ratio), na.rm = TRUE)
    # Geometric mean (inline to avoid robCompositions dependency)
    gm <- function(x) exp(mean(log(x[x > 0]), na.rm = TRUE))
    density_X <- log(density_X / gm(density_X))
    density_Y <- log(density_Y / gm(density_Y))
    distance <- sqrt(sum((density_X - density_Y)^2, na.rm=TRUE))
    # Compute the density ratio
    density_ratio_bayes <- density_X - density_Y
    kl_bayes <- length(density_X_orig) / 2 * log(mean(density_X_orig / density_Y_orig, na.rm = TRUE) * mean(density_Y_orig / density_X_orig, na.rm = TRUE))
    mean_ratio <- mean(density_ratio, na.rm = TRUE)
    sd_ratio <- sd(density_ratio, na.rm = TRUE)
    mean_ratio_bayes <- mean(density_ratio_bayes, na.rm = TRUE)
    sd_ratio_bayes <- sd(density_ratio_bayes, na.rm = TRUE)


    results <- list("predictions" = p,
                "ps_ratio" = ps_ratio,
                "ps_score" = ps,
                "cr" = cr,
                "mean_ps_x" = mean(p$prediction[1:mi_x], na.rm = TRUE),
                "mean_ps_y" = mean(p$prediction[(mi_x + 1):length(p$prediction)], na.rm = TRUE),
                "density_ratio" = density_ratio,
                "density_ratio_bayes" = density_ratio_bayes,
                "kl" = kl,
                "kl_bayes" = kl_bayes,
                "mean_ratio" = mean_ratio,
                "sd_ratio" = sd_ratio,
                "mean_ratio_bayes" = mean_ratio_bayes,
                "sd_ratio_bayes" = sd_ratio_bayes,
                "points" = points,
                "denX" = density_X_orig,
                "denY" = density_Y_orig,
                "bayesspace" = FALSE,
                "n_x" = n_x,
                "n_y" = n_y,
                "method" = method)
    class(results) <- "propscore"
    return(results)

}

#' Print method for propscore objects
#'
#' @param x an object of class "propscore"
#' @param ... additional arguments passed to the print method
#' @return The input object, invisibly.
#' @export
print.propscore <- function(x, ...){
  cat("mean propensity scores for x: ", x$mean_ps_x)
  cat("\nmean propensity scores for y: ", x$mean_ps_y)
  cat("\npropensity score statistics: ", x$ps_score, "\n")
}

#' Summary method for propscore objects
#'
#' @param object an object of class "propscore"
#' @param ... additional arguments passed to the summary method
#' @return An object of class \code{summary.propscore} with the following components:
#' \describe{
#'   \item{ps_score}{The propensity score statistic (pMSE).}
#'   \item{ps_ratio}{Ratio of predictions below vs above the class ratio.}
#'   \item{cr}{Class ratio (proportion of synthetic records in the combined set).}
#'   \item{mean_ps_x}{Mean predicted propensity for original records.}
#'   \item{mean_ps_y}{Mean predicted propensity for synthetic records.}
#'   \item{kl}{KL divergence between propensity density estimates.}
#'   \item{kl_bayes}{KL divergence in Bayes (log-ratio) space.}
#'   \item{mean_ratio}{Mean density ratio.}
#'   \item{sd_ratio}{Standard deviation of density ratio.}
#'   \item{n_x}{Number of original records.}
#'   \item{n_y}{Number of synthetic records.}
#'   \item{method}{Method used for propensity estimation.}
#' }
#' @export
summary.propscore <- function(object, ...) {
  res <- list(
    ps_score = object$ps_score,
    ps_ratio = object$ps_ratio,
    cr = object$cr,
    mean_ps_x = object$mean_ps_x,
    mean_ps_y = object$mean_ps_y,
    kl = object$kl,
    kl_bayes = object$kl_bayes,
    mean_ratio = object$mean_ratio,
    sd_ratio = object$sd_ratio,
    mean_ratio_bayes = object$mean_ratio_bayes,
    sd_ratio_bayes = object$sd_ratio_bayes,
    n_x = object$n_x,
    n_y = object$n_y,
    method = object$method
  )
  if (!is.null(object$oob_error)) {
    res$oob_error <- object$oob_error
    res$var_importance <- object$var_importance
  }
  if (!is.null(object$structure_ratio)) {
    res$within_orig_prox <- object$within_orig_prox
    res$within_synth_prox <- object$within_synth_prox
    res$cross_prox <- object$cross_prox
    res$structure_ratio <- object$structure_ratio
  }
  class(res) <- "summary.propscore"
  res
}

#' Print method for summary.propscore objects
#'
#' @param x an object of class "summary.propscore"
#' @param ... additional arguments passed to the print method
#' @return The input object, invisibly.
#' @export
print.summary.propscore <- function(x, ...) {
  cat("Propensity Score Utility Summary\n")
  cat("================================\n")
  cat("Method:          ", x$method, "\n")
  cat("Sample sizes:     n_original =", x$n_x, ", n_synthetic =", x$n_y, "\n")
  cat("Class ratio (cr):", round(x$cr, 4), "\n\n")
  cat("Propensity Score Statistic (pMSE):", format(x$ps_score, digits = 4), "\n")
  if (!is.na(x$ps_ratio)) {
    cat("PS ratio (below/above cr):       ", round(x$ps_ratio, 4), "\n\n")
  } else {
    cat("\n")
  }
  cat("Mean propensity (original): ", round(x$mean_ps_x, 4), "\n")
  cat("Mean propensity (synthetic):", round(x$mean_ps_y, 4), "\n\n")
  if (!is.null(x$oob_error)) {
    cat("OOB error:                  ", round(x$oob_error, 4), "\n")
  }
  if (!is.null(x$structure_ratio)) {
    cat("\nProximity structure:\n")
    cat("  Within-original:  ", round(x$within_orig_prox, 4), "\n")
    cat("  Within-synthetic: ", round(x$within_synth_prox, 4), "\n")
    cat("  Cross-class:      ", round(x$cross_prox, 4), "\n")
    cat("  Structure ratio:  ",
        if (is.na(x$structure_ratio)) "NA"
        else round(x$structure_ratio, 4),
        " (1 = indistinguishable)\n")
  }
  if (!is.na(x$kl)) {
    cat("\nDensity diagnostics:\n")
    cat("  KL divergence:              ", format(x$kl, digits = 4), "\n")
    cat("  KL divergence (Bayes space):", format(x$kl_bayes, digits = 4), "\n")
    cat("  Mean density ratio:         ", round(x$mean_ratio, 4),
        " (sd:", round(x$sd_ratio, 4), ")\n")
    cat("  Mean density ratio (Bayes): ", round(x$mean_ratio_bayes, 4),
        " (sd:", round(x$sd_ratio_bayes, 4), ")\n")
  }
  invisible(x)
}

#' Plot method for propscore objects
#'
#' @param x an object of class "propscore"
#' @param y not used
#' @param ... additional arguments passed to the plot method
#' @param which which plot to show: 1 for density, 2 for density ratio,
#' 3 for proximity structure (ranger only), 4 for variable importance
#' (ranger only)
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.propscore <- function(x, y = NULL, ..., which = 1){

  is_list_of_lists <- function(obj) {
    # Check if obj is a list
    if (!is.list(obj)) {
      return(FALSE)
    }

    # Check if every element of obj is a list
    for (element in obj) {
      if (!is.list(element)) {
        return(FALSE)
      }
    }

    return(TRUE)
  }

  show <- rep(FALSE, 4)
  show[which] <- TRUE

  if(is_list_of_lists(x)){
    if(x[[1]]$bayesspace){
      ylab1 <- "densities (Bayes space)"
      ylab2 <- "density ratio (Bayes space)"
    } else {
      lab1 <- "densities"
      lab2 <- "density ratio"
    }
    extractFromList <- function(x, what){
      df <- data.frame(sapply(x, function(x) x[[what]]))
      colnames(df) <- levels(x[[1]]$strata_x)
      df <- reshape2::melt(df)
      df$what <- what
      colnames(df)[1] <- "strata"
      return(df)
    }
    appendMe <- function(dfNames) {
      do.call(rbind, lapply(dfNames, function(x) {
        cbind(get(x), source = x)
      }))
    }

    denX <- extractFromList(x, what = "denX")
    denY <- extractFromList(x, what = "denY")
    points <- extractFromList(x, what = "points")
    denX$x <- points$value
    denY$x <- points$value
    #    den <- appendMe(c("denX", "denY"))
    den <- rbind(denX, denY)
    denRatio <- extractFromList(x, what = "density_ratio")
    denRatio$x <- points$value

    if(show[1L]){
      print(ggplot(den, aes(x = x, y = value, colour = what)) +
              geom_line() +
              facet_wrap(~ strata) +
              ylab(ylab2) +
              theme_minimal() +
              theme(legend.position="none") +
              geom_hline(yintercept = 0, color = "grey", linetype = "dashed"))
    }
    if(show[2L]){
      print(ggplot(denRatio, aes(x = x, y = value, colour = what)) +
              geom_line() +
              facet_wrap(~ strata) +
              ylab(ylab2) +
              theme_minimal() +
              theme(legend.position="none") +
              geom_hline(yintercept = 0, color = "grey", linetype = "dashed"))
    }

  } else {
    # Guard for ranger method: density/ratio plots need KDE fields
    if (is.null(x$points) && (show[1L] || show[2L])) {
      message("Density plots not available for method = 'ranger'. ",
              "Use which = 3 (proximity) or which = 4 (importance).")
      show[1L] <- FALSE
      show[2L] <- FALSE
    }

    if(x$bayesspace){
      ylab1 <- "densities (Bayes space)"
      ylab2 <- "density ratio (Bayes space)"
    } else {
      ylab1 <- "density"
      ylab2 <- "density ratio"
    }
    if(sum(show[1:2]) > 1){
      op <- par(mfrow = c(1,2))
      on.exit(par(op))
    }
    if(show[1L]){
      plot(x = x$points, y = x$denX, type = "l", ylab = ylab1, xlab = "")
      lines(x = x$points, y = x$denY, col = "red")
      abline(h = 1, lty = 2, col = "gray")
    }
    if(show[2L]){
      plot(x = x$points,
           y = x$density_ratio,
           type = "l",
           ylab = ylab2, xlab = "")
      abline(h = 1, lty = 2, col = "gray")
    }

    if (show[3L]) {
      # Proximity structure bar chart
      if (is.null(x$within_orig_prox)) {
        message("No proximity data (method != 'ranger' or proximity = 'none')")
      } else {
        grp_levels <- c("Within original", "Within synthetic", "Cross-class")
        df <- data.frame(
          group = factor(grp_levels, levels = grp_levels),
          proximity = c(x$within_orig_prox, x$within_synth_prox, x$cross_prox)
        )
        p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group,
                                               y = .data$proximity)) +
          ggplot2::geom_col(fill = "steelblue") +
          ggplot2::labs(title = "Proximity Structure",
                        x = NULL, y = "Mean proximity",
                        caption = if (is.na(x$structure_ratio)) "Structure ratio: NA"
                                  else sprintf("Structure ratio: %.3f",
                                               x$structure_ratio)) +
          ggplot2::theme_minimal()
        print(p)
      }
    }

    if (show[4L]) {
      # Variable importance bar chart
      if (is.null(x$var_importance)) {
        message("No importance data (importance = FALSE or method != 'ranger')")
      } else {
        imp <- sort(x$var_importance, decreasing = TRUE)
        df <- data.frame(var = factor(names(imp), levels = names(imp)),
                         importance = imp)
        p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$var,
                                               y = .data$importance)) +
          ggplot2::geom_col(fill = "darkgreen") +
          ggplot2::coord_flip() +
          ggplot2::labs(title = "Variable Importance",
                        x = NULL, y = "Importance") +
          ggplot2::theme_minimal()
        print(p)
      }
    }

  }
}







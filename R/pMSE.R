#' Propensity Score Mean Squared Error (pMSE) Utility Measure
#'
#' Computes the propensity score mean squared error (pMSE) and related utility
#' statistics for synthetic data. This is a comprehensive implementation of the
#' propensity score approach to measuring synthetic data utility, following
#' Snoke et al. (2018) and Woo et al. (2009).
#'
#' @param X data frame of original data, or a \code{\link{synth_pair}} object
#' @param Y data frame of synthetic data (not needed if X is a synth_pair)
#' @param vars character vector of variable names to use in the propensity model.
#'   If NULL, all common variables are used.
#' @param method character, propensity score estimation method:
#'   \itemize{
#'     \item "logit" - logistic regression (default, allows analytical null)
#'     \item "cart" - CART classification tree (requires rpart)
#'     \item "rf" - random forest
#'   }
#' @param maxorder integer, maximum order of interactions for logit method (default: 1).
#'   Use 0 for main effects only, 1 for two-way interactions, etc.
#' @param tree.method character, tree method when method="cart": "rpart" (default) or "ctree"
#' @param nperms integer, number of permutations for null distribution when using
#'   cart or rf methods (default: 20). Set to 0 to skip permutation test.
#' @param break.vars character vector of numeric variables to discretize for modeling.
#'   If NULL (default), continuous variables are used as-is for logit/rf, or
#'   automatically binned for cart.
#' @param ngroups integer, number of groups for discretizing continuous variables
#'   (default: 5)
#' @param na.rm logical, remove records with NA values (default: TRUE).
#'   Alternatively, use na="indicator" to add indicator variables for NAs.
#' @param na character, NA handling: "remove" (default), "indicator" (add indicators),
#'   or "impute" (impute with median/mode)
#' @param seed integer, random seed for reproducibility (default: NULL)
#' @param cp numeric, complexity parameter for rpart (default: 1e-8, smaller = more complex tree)
#' @param minbucket integer, minimum observations in terminal nodes for rpart (default: 5)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "pMSE" containing:
#' \itemize{
#'   \item pMSE: the propensity score mean squared error
#'   \item pMSE_null: expected pMSE under correct synthesis (NULL hypothesis)
#'   \item pMSE_ratio: ratio of observed to expected pMSE (values ~1 are good)
#'   \item S_pMSE: standardized pMSE (z-score, values ~0 are good)
#'   \item SPECKS: Kolmogorov-Smirnov statistic on propensity scores
#'   \item PO50: percent correctly classified at 0.5 threshold (values ~50% are good)
#'   \item ks_pvalue: p-value from KS test
#'   \item propensity_original: propensity scores for original records
#'   \item propensity_synthetic: propensity scores for synthetic records
#'   \item c: proportion of synthetic records in combined data
#'   \item k: number of parameters in propensity model (for logit)
#'   \item method: method used
#'   \item utility_interpretation: character, interpretation of results
#' }
#'
#' @details
#' The propensity score approach to utility assessment works by combining original
#' and synthetic data, fitting a model to predict which dataset each record belongs
#' to, and measuring how distinguishable the datasets are.
#'
#' \strong{The pMSE statistic:}
#' \deqn{pMSE = \frac{1}{N} \sum_{i=1}^{N} (\hat{p}_i - c)^2}
#' where \eqn{\hat{p}_i} is the estimated propensity score (probability of being
#' synthetic) and \eqn{c = n_Y/(n_X + n_Y)} is the true proportion of synthetic records.
#'
#' \strong{Null expectation (for logistic regression):}
#' \deqn{E[pMSE] = \frac{(k-1)(1-c)^2 c}{N}}
#' where \eqn{k} is the number of parameters and \eqn{N = n_X + n_Y}.
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item \strong{pMSE ratio ~ 1}: Correct synthesis - synthetic indistinguishable from original
#'   \item \strong{pMSE ratio < 1}: Acceptable (random variation)
#'   \item \strong{pMSE ratio > 2}: Poor utility - notable differences
#'   \item \strong{pMSE ratio > 5}: Very poor utility - substantial differences
#'   \item \strong{S_pMSE ~ 0}: Good utility
#'   \item \strong{|S_pMSE| > 2}: Significant difference from null
#'   \item \strong{SPECKS < 0.1}: Good utility
#'   \item \strong{PO50 ~ 50%}: Good utility (random classification)
#' }
#'
#' \strong{Methods:}
#' \itemize{
#'   \item \strong{logit}: Analytical null distribution available. Most interpretable.
#'     Use maxorder to control model complexity (0=main effects, 1=two-way interactions).
#'   \item \strong{cart}: Uses CART trees. Null estimated via permutation. Good for
#'     complex relationships but can be unstable.
#'   \item \strong{rf}: Random forest. Most flexible but computationally intensive.
#'     Null estimated via permutation.
#' }
#'
#' @seealso \code{\link{specks}} for a simpler SPECKS-focused function,
#'   \code{\link{propscore}} for an alternative implementation
#'
#' @references
#' Snoke, J., Raab, G.M., Nowok, B., Dibben, C., & Slavkovic, A. (2018).
#' General and specific utility measures for synthetic data.
#' \emph{Journal of the Royal Statistical Society: Series A}, 181(3), 663-688.
#' \doi{10.1111/rssa.12358}
#'
#' Woo, M.J., Reiter, J.P., Oganian, A., & Karr, A.F. (2009).
#' Global measures of data utility for microdata masked for disclosure limitation.
#' \emph{Journal of Privacy and Confidentiality}, 1(1), 111-124.
#'
#' @author Matthias Templ
#' @export
#' @importFrom stats glm binomial predict complete.cases as.formula ks.test median coef terms.formula quantile sd
#' @examples
#' # Create example data
#' set.seed(123)
#' n <- 500
#' original <- data.frame(
#'   age = sample(18:80, n, replace = TRUE),
#'   gender = factor(sample(c("M", "F"), n, replace = TRUE)),
#'   region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE)),
#'   income = exp(rnorm(n, 10, 1))
#' )
#'
#' # Good synthetic data (similar distributions)
#' synthetic_good <- data.frame(
#'   age = sample(18:80, n, replace = TRUE),
#'   gender = factor(sample(c("M", "F"), n, replace = TRUE)),
#'   region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE)),
#'   income = exp(rnorm(n, 10, 1))
#' )
#'
#' # Compute pMSE with logistic regression
#' result <- pMSE(original, synthetic_good, method = "logit", maxorder = 1)
#' print(result)
#' summary(result)
#' plot(result)
#'
#' # Poor synthetic data
#' synthetic_poor <- data.frame(
#'   age = sample(40:60, n, replace = TRUE),
#'   gender = factor(sample(c("M", "F"), n, replace = TRUE, prob = c(0.9, 0.1))),
#'   region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE, prob = c(0.7, 0.1, 0.1, 0.1))),
#'   income = exp(rnorm(n, 11, 0.5))
#' )
#'
#' result_poor <- pMSE(original, synthetic_poor, method = "logit")
#' print(result_poor)
#'
#' \donttest{
#' # Using CART method with permutation test
#' result_cart <- pMSE(original, synthetic_good, method = "cart", nperms = 10)
#' print(result_cart)
#' }
pMSE <- function(X, ...) {
  UseMethod("pMSE")
}

#' @rdname pMSE
#' @export
pMSE.synth_pair <- function(X, ...) {
  pMSE.default(
    X = X$original,
    Y = X$synthetic,
    vars = X$vars,
    ...
  )
}

#' @rdname pMSE
#' @export
pMSE.default <- function(X, Y,
                         vars = NULL,
                         method = c("logit", "cart", "rf"),
                         maxorder = 1,
                         tree.method = c("rpart", "ctree"),
                         nperms = 20,
                         break.vars = NULL,
                         ngroups = 5,
                         na = c("remove", "indicator", "impute"),
                         na.rm = TRUE,
                         seed = NULL,
                         cp = 1e-8,
                         minbucket = 5,
                         ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

  method <- match.arg(method)
  tree.method <- match.arg(tree.method)
  na <- match.arg(na)

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Determine variables to use
  if (is.null(vars)) {
    vars <- intersect(names(X), names(Y))
  }

  if (length(vars) == 0) {
    stop("No common variables found between X and Y.")
  }

  # Check variables exist
  missing_X <- setdiff(vars, names(X))
  missing_Y <- setdiff(vars, names(Y))
  if (length(missing_X) > 0) {
    stop(paste("Variables missing in X:", paste(missing_X, collapse = ", ")))
  }
  if (length(missing_Y) > 0) {
    stop(paste("Variables missing in Y:", paste(missing_Y, collapse = ", ")))
  }

  # Subset to selected variables
  X_sub <- X[, vars, drop = FALSE]
  Y_sub <- Y[, vars, drop = FALSE]

  # Combine datasets
  n_X <- nrow(X_sub)
  n_Y <- nrow(Y_sub)
  N <- n_X + n_Y

  combined <- rbind(X_sub, Y_sub)
  combined$..synth.. <- factor(c(rep(0, n_X), rep(1, n_Y)), levels = c("0", "1"))

  # Handle missing values
  if (na == "remove" || na.rm) {
    complete_idx <- complete.cases(combined[, vars, drop = FALSE])
    if (sum(!complete_idx) > 0) {
      message(sprintf("Removing %d records with missing values.", sum(!complete_idx)))
      combined <- combined[complete_idx, , drop = FALSE]
      # Recalculate sizes
      n_X <- sum(combined$..synth.. == "0")
      n_Y <- sum(combined$..synth.. == "1")
      N <- n_X + n_Y
    }
  } else if (na == "indicator") {
    # Add indicator variables for NAs (Rosenbaum & Rubin approach)
    for (v in vars) {
      if (any(is.na(combined[[v]]))) {
        ind_name <- paste0(v, "_NA")
        combined[[ind_name]] <- as.integer(is.na(combined[[v]]))
        # Impute the actual variable
        if (is.numeric(combined[[v]])) {
          combined[[v]][is.na(combined[[v]])] <- median(combined[[v]], na.rm = TRUE)
        } else {
          combined[[v]] <- addNA(as.factor(combined[[v]]))
        }
      }
    }
    vars <- setdiff(names(combined), "..synth..")
  } else if (na == "impute") {
    for (v in vars) {
      if (any(is.na(combined[[v]]))) {
        if (is.numeric(combined[[v]])) {
          combined[[v]][is.na(combined[[v]])] <- median(combined[[v]], na.rm = TRUE)
        } else {
          # Mode imputation for factors
          tbl <- table(combined[[v]], useNA = "no")
          mode_val <- names(tbl)[which.max(tbl)]
          combined[[v]][is.na(combined[[v]])] <- mode_val
        }
      }
    }
  }

  if (n_X == 0) stop("No original records remaining after NA handling.")
  if (n_Y == 0) stop("No synthetic records remaining after NA handling.")

  # Proportion of synthetic records
  c_prop <- n_Y / N

  # Discretize continuous variables if requested
  if (!is.null(break.vars)) {
    for (v in break.vars) {
      if (v %in% names(combined) && is.numeric(combined[[v]])) {
        combined[[v]] <- cut(combined[[v]],
                              breaks = quantile(combined[[v]], probs = seq(0, 1, length.out = ngroups + 1), na.rm = TRUE),
                              include.lowest = TRUE)
      }
    }
  }

  # Build formula
  model_vars <- setdiff(names(combined), "..synth..")

  if (method == "logit") {
    if (maxorder == 0) {
      form_str <- paste("..synth.. ~", paste(model_vars, collapse = " + "))
    } else {
      form_str <- paste("..synth.. ~ (", paste(model_vars, collapse = " + "), ")^", min(maxorder + 1, length(model_vars)))
    }
    form <- as.formula(form_str)

    # Fit logistic regression
    model <- tryCatch({
      glm(form, data = combined, family = binomial())
    }, error = function(e) {
      warning("Full model failed, falling back to main effects only: ", e$message)
      glm(as.formula(paste("..synth.. ~", paste(model_vars, collapse = " + "))),
          data = combined, family = binomial())
    })

    propensity <- predict(model, type = "response")
    k <- length(coef(model))

    # Analytical null distribution for logit
    pMSE_null <- (k - 1) * (1 - c_prop)^2 * c_prop / N
    pMSE_var_null <- 2 * (k - 1) * (1 - c_prop)^4 * c_prop^2 / N^2

  } else if (method == "cart") {
    if (!requireNamespace("rpart", quietly = TRUE)) {
      stop("Package 'rpart' required for method='cart'. Please install it.")
    }

    form <- as.formula(paste("..synth.. ~", paste(model_vars, collapse = " + ")))

    if (tree.method == "rpart") {
      model <- rpart::rpart(form, data = combined, method = "class",
                            control = rpart::rpart.control(cp = cp, minbucket = minbucket))
      propensity <- predict(model, type = "prob")[, "1"]
    } else {
      if (!requireNamespace("partykit", quietly = TRUE)) {
        stop("Package 'partykit' required for tree.method='ctree'. Please install it.")
      }
      model <- partykit::ctree(form, data = combined)
      propensity <- predict(model, type = "prob")[, "1"]
    }

    k <- NA

    # Permutation null for CART
    if (nperms > 0) {
      pMSE_perms <- numeric(nperms)
      for (p in seq_len(nperms)) {
        perm_data <- combined
        perm_data$..synth.. <- sample(perm_data$..synth..)
        if (tree.method == "rpart") {
          perm_model <- rpart::rpart(form, data = perm_data, method = "class",
                                     control = rpart::rpart.control(cp = cp, minbucket = minbucket))
          perm_prop <- predict(perm_model, type = "prob")[, "1"]
        } else {
          perm_model <- partykit::ctree(form, data = perm_data)
          perm_prop <- predict(perm_model, type = "prob")[, "1"]
        }
        pMSE_perms[p] <- mean((perm_prop - c_prop)^2)
      }
      pMSE_null <- mean(pMSE_perms)
      pMSE_var_null <- var(pMSE_perms)
    } else {
      pMSE_null <- NA
      pMSE_var_null <- NA
    }

  } else if (method == "rf") {
    form <- as.formula(paste("..synth.. ~", paste(model_vars, collapse = " + ")))
    model <- randomForest::randomForest(form, data = combined, ntree = 500)
    propensity <- predict(model, type = "prob")[, "1"]
    k <- NA

    # Permutation null for RF
    if (nperms > 0) {
      pMSE_perms <- numeric(nperms)
      for (p in seq_len(nperms)) {
        perm_data <- combined
        perm_data$..synth.. <- sample(perm_data$..synth..)
        perm_model <- randomForest::randomForest(form, data = perm_data, ntree = 200)
        perm_prop <- predict(perm_model, type = "prob")[, "1"]
        pMSE_perms[p] <- mean((perm_prop - c_prop)^2)
      }
      pMSE_null <- mean(pMSE_perms)
      pMSE_var_null <- var(pMSE_perms)
    } else {
      pMSE_null <- NA
      pMSE_var_null <- NA
    }
  }

  # Compute pMSE
  pMSE_val <- mean((propensity - c_prop)^2)

  # Compute pMSE ratio and standardized pMSE
  if (!is.na(pMSE_null) && pMSE_null > 0) {
    pMSE_ratio <- pMSE_val / pMSE_null
  } else {
    pMSE_ratio <- NA
  }

  if (!is.na(pMSE_var_null) && pMSE_var_null > 0) {
    S_pMSE <- (pMSE_val - pMSE_null) / sqrt(pMSE_var_null)
  } else {
    S_pMSE <- NA
  }

  # Split propensity scores
  synth_idx <- combined$..synth.. == "1"
  prop_original <- propensity[!synth_idx]
  prop_synthetic <- propensity[synth_idx]

  # Compute SPECKS (KS statistic)
  ks_result <- ks.test(prop_original, prop_synthetic)
  SPECKS <- as.numeric(ks_result$statistic)
  ks_pvalue <- ks_result$p.value

  # Compute PO50 (percent correctly classified at 0.5)
  predicted_class <- ifelse(propensity > 0.5, 1, 0)
  actual_class <- as.numeric(as.character(combined$..synth..))
  correct <- sum(predicted_class == actual_class)
  PO50 <- 100 * correct / N

  # Utility interpretation
  interpretation <- interpret_pMSE(pMSE_ratio, S_pMSE, SPECKS, PO50)

  # Build results
  results <- list(
    pMSE = pMSE_val,
    pMSE_null = pMSE_null,
    pMSE_var_null = pMSE_var_null,
    pMSE_ratio = pMSE_ratio,
    S_pMSE = S_pMSE,
    SPECKS = SPECKS,
    ks_pvalue = ks_pvalue,
    PO50 = PO50,
    propensity_original = prop_original,
    propensity_synthetic = prop_synthetic,
    c = c_prop,
    k = k,
    N = N,
    n_original = n_X,
    n_synthetic = n_Y,
    method = method,
    maxorder = if (method == "logit") maxorder else NA,
    nperms = if (method != "logit") nperms else NA,
    vars = model_vars,
    utility_interpretation = interpretation
  )

  class(results) <- "pMSE"
  return(results)
}


#' Interpret pMSE results
#' @keywords internal
interpret_pMSE <- function(pMSE_ratio, S_pMSE, SPECKS, PO50) {
  issues <- character(0)
  quality <- "GOOD"

  # Check pMSE ratio
  if (!is.na(pMSE_ratio)) {
    if (pMSE_ratio > 10) {
      issues <- c(issues, "pMSE ratio very high (>10)")
      quality <- "POOR"
    } else if (pMSE_ratio > 5) {
      issues <- c(issues, "pMSE ratio high (>5)")
      if (quality != "POOR") quality <- "FAIR"
    } else if (pMSE_ratio > 2) {
      issues <- c(issues, "pMSE ratio elevated (>2)")
      if (quality == "GOOD") quality <- "MODERATE"
    }
  }

  # Check S_pMSE
  if (!is.na(S_pMSE)) {
    if (abs(S_pMSE) > 3) {
      issues <- c(issues, "Standardized pMSE indicates significant difference")
      if (quality == "GOOD") quality <- "MODERATE"
    }
  }

  # Check SPECKS
  if (!is.na(SPECKS)) {
    if (SPECKS > 0.3) {
      issues <- c(issues, "SPECKS high (>0.3)")
      if (quality != "POOR") quality <- "FAIR"
    } else if (SPECKS > 0.2) {
      issues <- c(issues, "SPECKS elevated (>0.2)")
      if (quality == "GOOD") quality <- "MODERATE"
    }
  }

  # Check PO50
  if (!is.na(PO50)) {
    if (abs(PO50 - 50) > 20) {
      issues <- c(issues, "Classification accuracy far from 50%")
      if (quality != "POOR") quality <- "FAIR"
    } else if (abs(PO50 - 50) > 10) {
      issues <- c(issues, "Classification accuracy somewhat elevated")
      if (quality == "GOOD") quality <- "MODERATE"
    }
  }

  if (length(issues) == 0) {
    return(paste(quality, "- Synthetic data has good utility; distributions similar to original."))
  } else {
    return(paste(quality, "-", paste(issues, collapse = "; ")))
  }
}


#' Print method for pMSE objects
#'
#' @param x an object of class "pMSE"
#' @param ... additional arguments (ignored)
#' @export
print.pMSE <- function(x, ...) {
  cat("Propensity Score Mean Squared Error (pMSE) Utility\n")
  cat("===================================================\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original:", x$n_original, "| Synthetic:", x$n_synthetic, "\n")
  cat("  Combined (N):", x$N, "| Proportion synthetic (c):", sprintf("%.3f", x$c), "\n")
  cat("  Variables:", length(x$vars), "\n")
  cat("  Method:", x$method)
  if (!is.na(x$maxorder)) cat(" (maxorder =", x$maxorder, ")")
  if (!is.na(x$nperms) && x$nperms > 0) cat(" (", x$nperms, " permutations)", sep = "")
  cat("\n\n")

  cat("Utility Statistics:\n")
  cat("  pMSE:          ", sprintf("%.6f", x$pMSE), "\n")
  if (!is.na(x$pMSE_null)) {
    cat("  pMSE (null):   ", sprintf("%.6f", x$pMSE_null), "\n")
    cat("  pMSE ratio:    ", sprintf("%.2f", x$pMSE_ratio),
        ifelse(x$pMSE_ratio <= 2, " (good)", ifelse(x$pMSE_ratio <= 5, " (moderate)", " (poor)")), "\n")
  }
  if (!is.na(x$S_pMSE)) {
    cat("  S_pMSE:        ", sprintf("%.2f", x$S_pMSE),
        ifelse(abs(x$S_pMSE) <= 2, " (good)", " (elevated)"), "\n")
  }
  cat("  SPECKS:        ", sprintf("%.4f", x$SPECKS),
      ifelse(x$SPECKS <= 0.1, " (good)", ifelse(x$SPECKS <= 0.2, " (moderate)", " (poor)")), "\n")
  cat("  KS p-value:    ", sprintf("%.4f", x$ks_pvalue), "\n")
  cat("  PO50:          ", sprintf("%.1f%%", x$PO50),
      ifelse(abs(x$PO50 - 50) <= 10, " (good)", " (elevated)"), "\n\n")

  cat("Interpretation:\n")
  cat(" ", x$utility_interpretation, "\n")

  invisible(x)
}


#' Summary method for pMSE objects
#'
#' @param object an object of class "pMSE"
#' @param ... additional arguments (ignored)
#' @export
summary.pMSE <- function(object, ...) {
  prop_orig <- object$propensity_original
  prop_synth <- object$propensity_synthetic

  summ <- list(
    pMSE = object$pMSE,
    pMSE_null = object$pMSE_null,
    pMSE_ratio = object$pMSE_ratio,
    S_pMSE = object$S_pMSE,
    SPECKS = object$SPECKS,
    ks_pvalue = object$ks_pvalue,
    PO50 = object$PO50,
    prop_original_summary = c(
      mean = mean(prop_orig),
      sd = sd(prop_orig),
      min = min(prop_orig),
      q25 = unname(quantile(prop_orig, 0.25)),
      median = median(prop_orig),
      q75 = unname(quantile(prop_orig, 0.75)),
      max = max(prop_orig)
    ),
    prop_synthetic_summary = c(
      mean = mean(prop_synth),
      sd = sd(prop_synth),
      min = min(prop_synth),
      q25 = unname(quantile(prop_synth, 0.25)),
      median = median(prop_synth),
      q75 = unname(quantile(prop_synth, 0.75)),
      max = max(prop_synth)
    ),
    c = object$c,
    k = object$k,
    method = object$method,
    vars = object$vars,
    n_original = object$n_original,
    n_synthetic = object$n_synthetic,
    utility_interpretation = object$utility_interpretation
  )

  class(summ) <- "summary.pMSE"
  return(summ)
}


#' Print method for summary.pMSE objects
#'
#' @param x an object of class "summary.pMSE"
#' @param ... additional arguments (ignored)
#' @export
print.summary.pMSE <- function(x, ...) {
  cat("Summary: Propensity Score Mean Squared Error (pMSE)\n")
  cat("====================================================\n\n")

  cat("Method:", x$method, "\n")
  cat("Variables (", length(x$vars), "):", paste(x$vars, collapse = ", "), "\n")
  cat("Sample sizes: Original =", x$n_original, ", Synthetic =", x$n_synthetic, "\n\n")

  cat("Main Utility Statistics:\n")
  cat("  pMSE:        ", sprintf("%.6f", x$pMSE), "\n")
  if (!is.na(x$pMSE_null)) cat("  pMSE (null): ", sprintf("%.6f", x$pMSE_null), "\n")
  if (!is.na(x$pMSE_ratio)) cat("  pMSE ratio:  ", sprintf("%.2f", x$pMSE_ratio), "\n")
  if (!is.na(x$S_pMSE)) cat("  S_pMSE:      ", sprintf("%.2f", x$S_pMSE), "\n")
  cat("  SPECKS:      ", sprintf("%.4f", x$SPECKS), "\n")
  cat("  PO50:        ", sprintf("%.1f%%", x$PO50), "\n\n")

  cat("Propensity Score Distribution (Original):\n")
  cat("  Mean:", sprintf("%.4f", x$prop_original_summary["mean"]),
      " SD:", sprintf("%.4f", x$prop_original_summary["sd"]), "\n")
  cat("  Range: [", sprintf("%.4f", x$prop_original_summary["min"]),
      ",", sprintf("%.4f", x$prop_original_summary["max"]), "]\n")
  cat("  Quartiles:", sprintf("%.4f", x$prop_original_summary["q25"]),
      sprintf("%.4f", x$prop_original_summary["median"]),
      sprintf("%.4f", x$prop_original_summary["q75"]), "\n\n")

  cat("Propensity Score Distribution (Synthetic):\n")
  cat("  Mean:", sprintf("%.4f", x$prop_synthetic_summary["mean"]),
      " SD:", sprintf("%.4f", x$prop_synthetic_summary["sd"]), "\n")
  cat("  Range: [", sprintf("%.4f", x$prop_synthetic_summary["min"]),
      ",", sprintf("%.4f", x$prop_synthetic_summary["max"]), "]\n")
  cat("  Quartiles:", sprintf("%.4f", x$prop_synthetic_summary["q25"]),
      sprintf("%.4f", x$prop_synthetic_summary["median"]),
      sprintf("%.4f", x$prop_synthetic_summary["q75"]), "\n\n")

  cat("Interpretation:\n")
  cat(" ", x$utility_interpretation, "\n")

  invisible(x)
}


#' Plot method for pMSE objects
#'
#' @param x an object of class "pMSE"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer vector, which plots to show:
#'   \itemize{
#'     \item 1: Propensity score density comparison (default)
#'     \item 2: Propensity score ECDF comparison
#'     \item 3: Boxplot of propensity scores by group
#'     \item 4: Histogram comparison
#'   }
#' @importFrom graphics hist lines legend abline par boxplot
#' @importFrom stats density ecdf
#' @importFrom grDevices adjustcolor
#' @export
plot.pMSE <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 4)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    if (n_plots == 2) {
      op <- par(mfrow = c(1, 2))
    } else if (n_plots == 3) {
      op <- par(mfrow = c(1, 3))
    } else {
      op <- par(mfrow = c(2, 2))
    }
    on.exit(par(op))
  }

  prop_orig <- x$propensity_original
  prop_synth <- x$propensity_synthetic

  if (show[1]) {
    # Density comparison
    dens_orig <- density(prop_orig, from = 0, to = 1)
    dens_synth <- density(prop_synth, from = 0, to = 1)

    plot(dens_orig, main = "Propensity Score Densities",
         xlab = "Propensity Score", ylab = "Density",
         col = "steelblue", lwd = 2,
         xlim = c(0, 1),
         ylim = c(0, max(c(dens_orig$y, dens_synth$y)) * 1.1), ...)
    lines(dens_synth, col = "coral", lwd = 2)
    abline(v = x$c, col = "gray40", lty = 2, lwd = 1.5)
    legend("topright",
           legend = c("Original", "Synthetic", paste0("c = ", round(x$c, 2))),
           col = c("steelblue", "coral", "gray40"),
           lty = c(1, 1, 2), lwd = c(2, 2, 1.5), cex = 0.8)
  }

  if (show[2]) {
    # ECDF comparison
    ecdf_orig <- ecdf(prop_orig)
    ecdf_synth <- ecdf(prop_synth)

    plot(ecdf_orig, main = paste0("ECDF (SPECKS = ", round(x$SPECKS, 3), ")"),
         xlab = "Propensity Score", ylab = "Cumulative Probability",
         col = "steelblue", lwd = 2,
         xlim = c(0, 1), ...)
    lines(ecdf_synth, col = "coral", lwd = 2)
    legend("bottomright",
           legend = c("Original", "Synthetic"),
           col = c("steelblue", "coral"),
           lty = 1, lwd = 2, cex = 0.8)
  }

  if (show[3]) {
    # Boxplot
    all_props <- c(prop_orig, prop_synth)
    groups <- factor(c(rep("Original", length(prop_orig)),
                       rep("Synthetic", length(prop_synth))))

    boxplot(all_props ~ groups,
            main = paste0("Propensity Scores (PO50 = ", round(x$PO50, 1), "%)"),
            ylab = "Propensity Score",
            col = c("steelblue", "coral"),
            ylim = c(0, 1), ...)
    abline(h = x$c, col = "gray40", lty = 2, lwd = 1.5)
    abline(h = 0.5, col = "darkgreen", lty = 3, lwd = 1.5)
  }

  if (show[4]) {
    # Histogram comparison
    breaks <- seq(0, 1, by = 0.05)

    hist(prop_orig, breaks = breaks, col = adjustcolor("steelblue", 0.5),
         main = "Propensity Score Distributions",
         xlab = "Propensity Score", ylab = "Frequency",
         xlim = c(0, 1), ...)
    hist(prop_synth, breaks = breaks, col = adjustcolor("coral", 0.5), add = TRUE)
    abline(v = x$c, col = "gray40", lty = 2, lwd = 2)
    legend("topright",
           legend = c("Original", "Synthetic"),
           fill = adjustcolor(c("steelblue", "coral"), 0.5), cex = 0.8)
  }
}

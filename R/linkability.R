#' Linkability Risk
#'
#' Computes the Linkability Risk for synthetic data, the second of three explicit
#' GDPR anonymization failure criteria (Article 29 Working Party). An attacker
#' who holds two disjoint subsets of columns from the original data tries to link
#' records across these subsets using the synthetic data as a bridge.
#'
#' @param X data frame of original/training data
#' @param Y data frame of synthetic data
#' @param holdout data frame of holdout data (optional). If NULL, a holdout set
#'   is automatically created by splitting X.
#' @param holdout_fraction numeric, fraction of X to use as holdout if holdout
#'   is NULL (default: 0.5)
#' @param n_attacks integer, number of synthetic records to use for the attack
#'   (default: 2000). If larger than nrow(Y), all synthetic records are used.
#' @param n_neighbors integer, number of nearest neighbors to consider for a
#'   successful link (default: 1). With \code{n_neighbors = 1}, only exact NN
#'   matches count. With \code{n_neighbors > 1}, a link succeeds if the
#'   auxiliary NN is among the top-k secret NNs.
#' @param aux_cols character vector of column names for the auxiliary (attacker-known)
#'   subset. If NULL (default), columns are randomly split into two equal halves.
#' @param vars character vector of variable names to use. If NULL (default),
#'   all common variables between X and Y are used.
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param seed integer, random seed for reproducibility (default: NULL)
#' @param confidence_level numeric, confidence level for Wilson score intervals
#'   (default: 0.95)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "linkability" containing:
#' \itemize{
#'   \item risk: residual risk score, bounded between 0 and 1
#'   \item risk_ci: confidence interval for risk
#'   \item risk_attack: attack success rate (fraction of links in original)
#'   \item risk_attack_ci: Wilson CI for attack success rate
#'   \item risk_control: control success rate (fraction of links in holdout, NA if none)
#'   \item risk_control_ci: Wilson CI for control success rate (NA if none)
#'   \item n_attacks: number of attacks performed
#'   \item n_success: number of successful links in original
#'   \item n_control_success: number of successful links in holdout
#'   \item links: logical vector of link successes per attack (original)
#'   \item links_control: logical vector of link successes per attack (holdout)
#'   \item privacy_pass: logical, risk <= 0.1
#'   \item n_original: number of original records (training portion)
#'   \item n_synthetic: number of synthetic records
#'   \item n_holdout: number of holdout records
#'   \item aux_cols: auxiliary columns used
#'   \item secret_cols: secret columns used
#'   \item n_neighbors: number of neighbors considered
#'   \item vars: all variables used
#'   \item confidence_level: confidence level used
#' }
#'
#' @details
#' The linkability attack simulates an adversary who holds two disjoint sets of
#' columns from the original data (e.g., demographics and medical records) and
#' attempts to re-link them using the synthetic data.
#'
#' For each sampled synthetic record:
#' \enumerate{
#'   \item Find its nearest neighbor in the original data using only auxiliary columns
#'   \item Find its nearest neighbor in the original data using only secret columns
#'   \item If both nearest neighbors point to the same original record, the link
#'     is successful
#' }
#'
#' With \code{n_neighbors > 1}, a link succeeds if the auxiliary nearest neighbor
#' is among the top-k secret nearest neighbors, making the attack more powerful.
#'
#' The risk score uses Wilson score intervals for robust estimation:
#' \deqn{r_{attack} = n_{success} / n_{attacks}}
#' \deqn{r_{control} = n_{control\_success} / n_{attacks}}
#' \deqn{risk = (r_{attack} - r_{control}) / (1 - r_{control})}
#'
#' Without holdout data, the absolute attack rate is reported as the risk.
#'
#' Gower distance is used for nearest-neighbor computation, which handles mixed
#' (numeric and categorical) variable types.
#'
#' @section Column splitting:
#' When \code{aux_cols} is not specified, the available columns are randomly split
#' into two approximately equal groups. The split is controlled by \code{seed}.
#' For reproducible results, always set the seed.
#'
#' When \code{aux_cols} is specified, the remaining columns form the secret set.
#' At least one column must be in each group.
#'
#' @section Holdout splitting:
#' When no external holdout is provided, the original data is split internally.
#' The holdout serves as a control: links found in the holdout represent baseline
#' linkability due to data structure rather than information leakage. The residual
#' risk subtracts this baseline.
#'
#' @seealso \code{\link{singling_out}} for singling out risk,
#'   \code{\link{dcap}} for differential correct attribution probability,
#'   \code{\link{tcap}} for targeted CAP
#'
#' @references
#' Giomi, M., Boenisch, F., Wehmeyer, C. & Tasnadi, B. (2023).
#' A Unified Framework for Quantifying Privacy Risk in Synthetic Data.
#' \emph{Proceedings on Privacy Enhancing Technologies (PoPETs)}, 2023(2), 312--328.
#' \doi{10.56553/popets-2023-0055}
#'
#' Article 29 Data Protection Working Party (2014).
#' Opinion 05/2014 on Anonymisation Techniques. WP216.
#'
#' @author Matthias Templ
#' @family privacy-models
#' @export
#' @importFrom VIM gowerD
#' @importFrom stats complete.cases median qnorm
#' @importFrom graphics hist abline legend barplot par text
#' @examples
#' # Create example data
#' set.seed(123)
#' n <- 200
#' X <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   income = rnorm(n, 50000, 15000),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
#' )
#'
#' # Good synthetic data (independent, no memorization)
#' Y <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   income = rnorm(n, 50000, 15000),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
#' )
#'
#' result <- linkability(X, Y, n_attacks = 500, seed = 42)
#' print(result)
#' summary(result)
#'
#' \donttest{
#' # Memorized data (Y is copy of X) - should show high risk
#' Y_copy <- X[sample(nrow(X), n, replace = TRUE), ]
#' result_bad <- linkability(X, Y_copy, n_attacks = 500, seed = 42)
#' print(result_bad)
#' }
linkability <- function(X, ...) {
  UseMethod("linkability")
}

#' @rdname linkability
#' @export
linkability.synth_pair <- function(X, ...) {
  linkability.default(
    X = X$original,
    Y = X$synthetic,
    holdout = X$holdout,
    vars = X$vars,
    ...
  )
}

#' @rdname linkability
#' @export
linkability.default <- function(X, Y,
                                 holdout = NULL,
                                 holdout_fraction = 0.5,
                                 n_attacks = 2000,
                                 n_neighbors = 1,
                                 aux_cols = NULL,
                                 vars = NULL,
                                 na.rm = TRUE,
                                 seed = NULL,
                                 confidence_level = 0.95,
                                 ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")
  if (!is.null(holdout) && !is.data.frame(holdout)) {
    stop("holdout must be a data frame or NULL.")
  }

  if (!is.numeric(n_attacks) || length(n_attacks) != 1 || n_attacks < 1) {
    stop("n_attacks must be a positive integer.")
  }
  n_attacks <- as.integer(n_attacks)

  if (!is.numeric(n_neighbors) || length(n_neighbors) != 1 || n_neighbors < 1) {
    stop("n_neighbors must be a positive integer.")
  }
  n_neighbors <- as.integer(n_neighbors)

  if (!is.numeric(confidence_level) || length(confidence_level) != 1 ||
      confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be a number between 0 and 1.")
  }

  # Determine variables to use
  if (is.null(vars)) {
    vars <- intersect(names(X), names(Y))
    if (!is.null(holdout)) {
      vars <- intersect(vars, names(holdout))
    }
  }

  if (length(vars) < 2) {
    stop("At least 2 variables are required for linkability analysis.")
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

  # Check variable types match
  for (var in vars) {
    if (!identical(class(X[[var]]), class(Y[[var]]))) {
      stop(paste("Variable", var, "has different class in X and Y."))
    }
  }

  # Subset to selected variables
  X <- X[, vars, drop = FALSE]
  Y <- Y[, vars, drop = FALSE]

  # Handle missing values
  if (na.rm) {
    complete_X <- complete.cases(X)
    complete_Y <- complete.cases(Y)
    X <- X[complete_X, , drop = FALSE]
    Y <- Y[complete_Y, , drop = FALSE]
  }

  if (nrow(X) == 0) stop("No complete cases in X after removing NAs.")
  if (nrow(Y) == 0) stop("No complete cases in Y after removing NAs.")

  # Split columns into aux and secret
  if (is.null(aux_cols)) {
    if (!is.null(seed)) set.seed(seed)
    n_aux <- max(1, floor(length(vars) / 2))
    aux_cols <- sample(vars, n_aux)
  } else {
    # Validate user-specified aux_cols
    missing_aux <- setdiff(aux_cols, vars)
    if (length(missing_aux) > 0) {
      stop(paste("aux_cols not in vars:", paste(missing_aux, collapse = ", ")))
    }
  }

  secret_cols <- setdiff(vars, aux_cols)
  if (length(secret_cols) == 0) {
    stop("secret_cols is empty: aux_cols must not cover all variables.")
  }
  if (length(aux_cols) == 0) {
    stop("aux_cols must contain at least one variable.")
  }

  # Create or validate holdout
  if (is.null(holdout)) {
    if (!is.null(seed)) set.seed(seed + 1L)
    n_holdout <- max(1, floor(nrow(X) * holdout_fraction))
    holdout_idx <- sample(nrow(X), n_holdout)
    holdout <- X[holdout_idx, , drop = FALSE]
    train <- X[-holdout_idx, , drop = FALSE]
  } else {
    holdout <- holdout[, vars, drop = FALSE]
    if (na.rm) {
      complete_H <- complete.cases(holdout)
      holdout <- holdout[complete_H, , drop = FALSE]
    }
    if (nrow(holdout) == 0) stop("No complete cases in holdout after removing NAs.")
    train <- X
  }

  n_train <- nrow(train)
  n_synthetic <- nrow(Y)
  n_holdout_final <- nrow(holdout)

  # Cap n_neighbors at dataset sizes
  n_neighbors <- min(n_neighbors, n_train, n_holdout_final)

  # Sample synthetic records for attack
  if (!is.null(seed)) set.seed(seed + 2L)
  if (n_attacks >= n_synthetic) {
    attack_idx <- seq_len(n_synthetic)
    n_attacks <- n_synthetic
  } else {
    attack_idx <- sample(n_synthetic, n_attacks)
  }
  Y_attack <- Y[attack_idx, , drop = FALSE]

  # --- Attack on training data ---
  links <- .linkability_attack(Y_attack, train, aux_cols, secret_cols,
                               n_neighbors)

  # --- Attack on holdout (control) ---
  links_control <- .linkability_attack(Y_attack, holdout, aux_cols, secret_cols,
                                       n_neighbors)

  n_success <- sum(links)
  n_control_success <- sum(links_control)

  # Compute residual risk with Wilson CIs
  rr <- .residual_risk(n_success, n_control_success, n_attacks, confidence_level)

  results <- list(
    risk = rr$risk,
    risk_ci = rr$risk_ci,
    risk_attack = rr$risk_attack,
    risk_attack_ci = rr$risk_attack_ci,
    risk_control = rr$risk_control,
    risk_control_ci = rr$risk_control_ci,
    n_attacks = n_attacks,
    n_success = n_success,
    n_control_success = n_control_success,
    links = links,
    links_control = links_control,
    privacy_pass = rr$privacy_pass,
    n_original = n_train,
    n_synthetic = n_synthetic,
    n_holdout = n_holdout_final,
    aux_cols = aux_cols,
    secret_cols = secret_cols,
    n_neighbors = n_neighbors,
    vars = vars,
    confidence_level = confidence_level
  )

  class(results) <- "linkability"
  return(results)
}

# --- Internal helpers ---

#' Run linkability attack: find NNs via aux and secret columns, check for links
#'
#' @param Y_attack data frame of synthetic attack records
#' @param ref_data data frame of reference data (training or holdout)
#' @param aux_cols character vector of auxiliary column names
#' @param secret_cols character vector of secret column names
#' @param n_neighbors integer, number of neighbors to consider
#' @return logical vector of length nrow(Y_attack), TRUE if link succeeded
#' @keywords internal
.linkability_attack <- function(Y_attack, ref_data, aux_cols, secret_cols,
                                n_neighbors) {
  n <- nrow(Y_attack)
  n_ref <- nrow(ref_data)

  if (n_ref == 0) return(rep(FALSE, n))

  # Compute Gower distances for aux columns
  dist_aux <- VIM::gowerD(Y_attack[, aux_cols, drop = FALSE],
                           ref_data[, aux_cols, drop = FALSE])

  # Compute Gower distances for secret columns
  dist_secret <- VIM::gowerD(Y_attack[, secret_cols, drop = FALSE],
                              ref_data[, secret_cols, drop = FALSE])

  if (n_neighbors == 1) {
    # Simple case: check if NN from aux == NN from secret
    nn_aux <- apply(dist_aux, 1, which.min)
    nn_secret <- apply(dist_secret, 1, which.min)
    return(nn_aux == nn_secret)
  }

  # Top-k case: check if aux NN is among top-k secret NNs
  nn_aux <- apply(dist_aux, 1, which.min)
  k <- min(n_neighbors, n_ref)
  vapply(seq_len(n), function(i) {
    top_k <- order(dist_secret[i, ])[seq_len(k)]
    nn_aux[i] %in% top_k
  }, logical(1))
}

# --- S3 methods ---

#' Print method for linkability objects
#'
#' @param x an object of class "linkability"
#' @param ... additional arguments (ignored)
#' @export
print.linkability <- function(x, ...) {
  cat("Linkability Risk Assessment\n")
  cat("===========================\n")
  cat("Auxiliary columns:", paste(x$aux_cols, collapse = ", "), "\n")
  cat("Secret columns:", paste(x$secret_cols, collapse = ", "), "\n")
  cat("Neighbors considered:", x$n_neighbors, "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original (training):", x$n_original, "\n")
  cat("  Holdout:", x$n_holdout, "\n")
  cat("  Synthetic:", x$n_synthetic, "\n\n")

  cat("Attack Results (", x$n_attacks, " attacks):\n", sep = "")
  cat("  Successful links in original:", x$n_success,
      sprintf("(%.1f%%)", 100 * x$risk_attack), "\n")
  cat("  Successful links in holdout:", x$n_control_success,
      sprintf("(%.1f%%)", 100 * x$risk_control), "\n\n")

  cat("Risk Score:\n")
  cat("  Residual risk:", round(x$risk, 4),
      sprintf("[%.4f, %.4f]", x$risk_ci["lower"], x$risk_ci["upper"]),
      sprintf("(%d%% CI)\n", round(100 * x$confidence_level)))

  cat("\nPrivacy Assessment:")
  if (x$privacy_pass) {
    cat(" PASS\n")
    cat("  Linkability risk is within acceptable bounds (<= 0.1).\n")
  } else {
    cat(" WARNING\n")
    cat("  Elevated linkability risk detected (> 0.1).\n")
    cat("  Synthetic data may enable cross-linking of individual records.\n")
  }

  invisible(x)
}

#' Summary method for linkability objects
#'
#' @param object an object of class "linkability"
#' @param ... additional arguments (ignored)
#' @export
summary.linkability <- function(object, ...) {
  summ <- list(
    risk = object$risk,
    risk_ci = object$risk_ci,
    risk_attack = object$risk_attack,
    risk_attack_ci = object$risk_attack_ci,
    risk_control = object$risk_control,
    risk_control_ci = object$risk_control_ci,
    n_attacks = object$n_attacks,
    n_success = object$n_success,
    n_control_success = object$n_control_success,
    success_rate_attack = object$n_success / object$n_attacks,
    success_rate_control = object$n_control_success / object$n_attacks,
    privacy_pass = object$privacy_pass,
    n_original = object$n_original,
    n_synthetic = object$n_synthetic,
    n_holdout = object$n_holdout,
    aux_cols = object$aux_cols,
    secret_cols = object$secret_cols,
    n_neighbors = object$n_neighbors,
    vars = object$vars,
    confidence_level = object$confidence_level
  )

  class(summ) <- "summary.linkability"
  return(summ)
}

#' Print method for summary.linkability objects
#'
#' @param x an object of class "summary.linkability"
#' @param ... additional arguments (ignored)
#' @export
print.summary.linkability <- function(x, ...) {
  cat("Summary: Linkability Risk Assessment\n")
  cat("====================================\n")
  cat("Auxiliary columns (", length(x$aux_cols), "):",
      paste(x$aux_cols, collapse = ", "), "\n")
  cat("Secret columns (", length(x$secret_cols), "):",
      paste(x$secret_cols, collapse = ", "), "\n")
  cat("Neighbors:", x$n_neighbors, "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original:", x$n_original, "| Holdout:", x$n_holdout,
      "| Synthetic:", x$n_synthetic, "\n\n")

  cat("Key Metrics:\n")
  cat("  Residual risk:", round(x$risk, 4),
      sprintf("[%.4f, %.4f]", x$risk_ci["lower"], x$risk_ci["upper"]), "\n")
  cat("  Attack rate:", sprintf("%.4f", x$risk_attack),
      sprintf("[%.4f, %.4f]", x$risk_attack_ci["lower"],
              x$risk_attack_ci["upper"]), "\n")
  cat("  Control rate:", sprintf("%.4f", x$risk_control),
      sprintf("[%.4f, %.4f]", x$risk_control_ci["lower"],
              x$risk_control_ci["upper"]), "\n")
  cat("  Privacy:", ifelse(x$privacy_pass, "PASS", "WARNING"), "\n\n")

  cat("Link Results:\n")
  cat("  Links in original:", x$n_success, "/", x$n_attacks,
      sprintf("(%.1f%%)\n", 100 * x$success_rate_attack))
  cat("  Links in holdout:", x$n_control_success, "/", x$n_attacks,
      sprintf("(%.1f%%)\n", 100 * x$success_rate_control))

  invisible(x)
}

#' Plot method for linkability objects
#'
#' @param x an object of class "linkability"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot(s) to show: 1 = Risk comparison barplot (attack vs control),
#'   2 = Success rate by column group (aux vs secret column counts)
#' @export
plot.linkability <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Risk comparison barplot
    vals <- c(x$risk_attack, x$risk_control, x$risk)
    names_vals <- c("Attack", "Control", "Residual")
    cols <- c("coral", "steelblue", "darkgreen")

    bp <- barplot(vals, names.arg = names_vals,
                  main = paste("Linkability Risk\n(",
                               x$n_attacks, "attacks)"),
                  ylab = "Rate",
                  col = cols,
                  ylim = c(0, max(0.2, max(vals) * 1.3)), ...)
    abline(h = 0.1, col = "red", lwd = 2, lty = 2)
    text(bp, vals + max(vals) * 0.05,
         labels = round(vals, 3), cex = 0.9)
    legend("topright", "threshold = 0.1", col = "red", lty = 2, lwd = 2,
           cex = 0.8)
  }

  if (show[2]) {
    # Column group comparison
    vals <- c(length(x$aux_cols), length(x$secret_cols))
    names_vals <- c("Auxiliary", "Secret")
    cols <- c("coral", "steelblue")

    bp <- barplot(vals, names.arg = names_vals,
                  main = paste("Column Split\n",
                               length(x$aux_cols), "aux +",
                               length(x$secret_cols), "secret"),
                  ylab = "Number of columns",
                  col = cols, ...)
    text(bp, vals + 0.2, labels = vals, cex = 0.9)
  }
}

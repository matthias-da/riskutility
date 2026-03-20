#' Nearest-Neighbor Adversarial Accuracy (NNAA)
#'
#' Computes the Nearest-Neighbor Adversarial Accuracy and Privacy Loss metrics
#' for synthetic data evaluation. NNAA uses a 1-nearest-neighbor two-sample test
#' to assess whether synthetic data is distinguishable from real data.
#' Privacy Loss compares adversarial accuracy on training versus holdout data
#' to detect memorization.
#'
#' @param X data frame of original/training data
#' @param Y data frame of synthetic data
#' @param holdout data frame of holdout data (optional). If NULL, a holdout set
#'   is automatically created by splitting X.
#' @param holdout_fraction numeric, fraction of X to use as holdout if holdout
#'   is NULL (default: 0.5)
#' @param vars character vector of variable names to use for distance calculation.
#'   If NULL (default), all common variables between X, Y, and holdout are used.
#' @param method character, distance method: "gower" (default, handles mixed types)
#'   or "euclidean" (numerical variables only)
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param seed integer, random seed for holdout sampling (default: NULL)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "nnaa" containing:
#' \itemize{
#'   \item aa_train: adversarial accuracy on training data vs synthetic
#'   \item aa_holdout: adversarial accuracy on holdout data vs synthetic (NA if no holdout)
#'   \item privacy_loss: aa_holdout - aa_train (NA if no holdout)
#'   \item privacy_pass: logical, TRUE if privacy_loss <= 0.1
#'   \item aa_train_left: proportion where real points are closer to other real than to synthetic
#'   \item aa_train_right: proportion where synthetic points are closer to real than to other synthetic
#'   \item aa_holdout_left: holdout left component (NA if no holdout)
#'   \item aa_holdout_right: holdout right component (NA if no holdout)
#'   \item d_TS: nearest-neighbor distances from training to synthetic (per record)
#'   \item d_TT: nearest-neighbor distances within training, self-excluded (per record)
#'   \item n_train, n_synthetic, n_holdout: dataset sizes
#'   \item method, vars, holdout_fraction: parameters used
#' }
#'
#' @details
#' Adversarial Accuracy (AA) measures how well a 1-nearest-neighbor classifier
#' can distinguish real from synthetic data. It is defined as:
#'
#' \deqn{AA(T, S) = 0.5 \cdot (\text{share}(d_{TS} > d_{TT}) + \text{share}(d_{ST} > d_{SS}))}
#'
#' where:
#' \itemize{
#'   \item \eqn{d_{TS}(i)}: distance from real point i to its nearest synthetic neighbor
#'   \item \eqn{d_{TT}(i)}: distance from real point i to its nearest real neighbor (self-excluded)
#'   \item \eqn{d_{ST}(i)}: distance from synthetic point i to its nearest real neighbor
#'   \item \eqn{d_{SS}(i)}: distance from synthetic point i to its nearest synthetic neighbor (self-excluded)
#' }
#'
#' Interpretation of AA:
#' \itemize{
#'   \item \strong{AA ~ 0.5}: Ideal - real and synthetic are indistinguishable (good utility)
#'   \item \strong{AA < 0.5}: Synthetic is too close to real - possible memorization (privacy concern)
#'   \item \strong{AA > 0.5}: Synthetic differs from real - poor utility
#' }
#'
#' Privacy Loss compares AA on training versus holdout data:
#' \deqn{Privacy Loss = AA(holdout, S) - AA(train, S)}
#'
#' Interpretation of Privacy Loss:
#' \itemize{
#'   \item \strong{~0}: Good privacy - synthetic is equally close to training and holdout
#'   \item \strong{> 0}: Privacy concern - synthetic is closer to training than holdout
#'   \item \strong{~ 0.5}: Severe memorization detected
#' }
#'
#' @section Holdout splitting:
#' When no external holdout is provided, NNAA internally splits the original data.
#' With default \code{holdout_fraction = 0.5}, half of X becomes the holdout and
#' half is used as the training set. This is standard practice for NNAA.
#' For best results, provide an actual holdout set from the synthesis process.
#'
#' @seealso \code{\link{dcr}} for distance to closest record,
#'   \code{\link{nndr}} for nearest neighbor distance ratio,
#'   \code{\link{ims}} for identical match share
#'
#' @references
#' Yale, A., Dash, S., Dutta, R., Guyon, I., Pavao, A. & Bennett, K.P. (2020).
#' Generation and Evaluation of Privacy Preserving Synthetic Health Data.
#' \emph{Neurocomputing}, 416, 244--255.
#' \doi{10.1016/j.neucom.2019.12.136}
#'
#' @author Matthias Templ
#' @family distance-risk
#' @export
#' @importFrom VIM gowerD
#' @importFrom stats complete.cases quantile sd median
#' @importFrom graphics hist abline legend barplot par
#' @examples
#' # Create example data
#' set.seed(123)
#' n <- 200
#' X <- data.frame(
#'   age = rnorm(n, 40, 10),
#'   income = rnorm(n, 50000, 15000),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
#' )
#'
#' # Good synthetic data (random, no memorization)
#' Y_good <- data.frame(
#'   age = rnorm(n, 40, 10),
#'   income = rnorm(n, 50000, 15000),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
#' )
#'
#' # Compute NNAA
#' result <- nnaa(X, Y_good, seed = 42)
#' print(result)
#' summary(result)
nnaa <- function(X, ...) {
  UseMethod("nnaa")
}

#' @rdname nnaa
#' @export
nnaa.synth_pair <- function(X, ...) {
  if (!is.null(X$source) && X$source == "sdcMicro") {
    stop("Package 'sdcMicro' is required for nnaa() but is only loaded for ecosystem integration.\n",
         "Please install it: install.packages('sdcMicro')\n",
         "Alternatives: dcr(), nndr(), ims() for distance-based privacy evaluation",
         call. = FALSE)
  }
  nnaa.default(
    X = X$original,
    Y = X$synthetic,
    holdout = X$holdout,
    vars = X$vars,
    ...
  )
}

#' @rdname nnaa
#' @export
nnaa.default <- function(X, Y,
                         holdout = NULL,
                         holdout_fraction = 0.5,
                         vars = NULL,
                         method = c("gower", "euclidean"),
                         na.rm = TRUE,
                         seed = NULL,
                         ...) {

  method <- match.arg(method)

  # Shared input validation, variable intersection, NA removal, holdout split
  v <- .validate_pair_inputs(X, Y, holdout = holdout,
                             holdout_fraction = holdout_fraction,
                             vars = vars, na.rm = na.rm, seed = seed)
  Y       <- v$Y
  vars    <- v$vars
  train   <- v$train
  holdout <- v$holdout

  n_synthetic <- nrow(Y)
  n_train <- nrow(train)
  n_holdout <- nrow(holdout)

  # Helper: get min distance per row from a distance matrix
  row_min <- function(mat) {
    apply(mat, 1, min, na.rm = TRUE)
  }

  # Helper: get min distance per row with self-exclusion (set diag to Inf)
  row_min_self_excluded <- function(mat) {
    diag(mat) <- Inf
    apply(mat, 1, min, na.rm = TRUE)
  }

  # Compute all needed distance matrices
  if (method == "gower") {
    # --- Training AA ---
    # d_TS: training -> synthetic (n_train x n_synthetic)
    dist_TS <- VIM::gowerD(train, Y)
    # d_TT: training -> training (n_train x n_train), self-excluded
    dist_TT <- VIM::gowerD(train, train)
    # d_ST: synthetic -> training (n_synthetic x n_train)
    # This is the transpose of dist_TS
    dist_ST <- t(dist_TS)
    # d_SS: synthetic -> synthetic (n_synthetic x n_synthetic), self-excluded
    dist_SS <- VIM::gowerD(Y, Y)

    d_TS <- row_min(dist_TS)
    d_TT <- row_min_self_excluded(dist_TT)
    d_ST <- row_min(dist_ST)
    d_SS <- row_min_self_excluded(dist_SS)

    # --- Holdout AA ---
    # d_HS: holdout -> synthetic (n_holdout x n_synthetic)
    dist_HS <- VIM::gowerD(holdout, Y)
    # d_HH: holdout -> holdout (n_holdout x n_holdout), self-excluded
    dist_HH <- VIM::gowerD(holdout, holdout)
    # d_SH: synthetic -> holdout (n_synthetic x n_holdout)
    dist_SH <- t(dist_HS)

    d_HS <- row_min(dist_HS)
    d_HH <- row_min_self_excluded(dist_HH)
    d_SH <- row_min(dist_SH)
    # Reuse d_SS for holdout AA right component

  } else if (method == "euclidean") {
    # .normalize_and_split / .euclidean_dist are defined in R/utils_internal.R
    norms <- .normalize_and_split(train, holdout, Y)
    train_norm <- norms[[1]]
    holdout_norm <- norms[[2]]
    Y_norm <- norms[[3]]

    # Training AA distances
    dist_TS <- .euclidean_dist(train_norm, Y_norm)
    dist_TT <- .euclidean_dist(train_norm, train_norm)
    dist_ST <- t(dist_TS)
    dist_SS <- .euclidean_dist(Y_norm, Y_norm)

    d_TS <- row_min(dist_TS)
    d_TT <- row_min_self_excluded(dist_TT)
    d_ST <- row_min(dist_ST)
    d_SS <- row_min_self_excluded(dist_SS)

    # Holdout AA distances
    dist_HS <- .euclidean_dist(holdout_norm, Y_norm)
    dist_HH <- .euclidean_dist(holdout_norm, holdout_norm)
    dist_SH <- t(dist_HS)

    d_HS <- row_min(dist_HS)
    d_HH <- row_min_self_excluded(dist_HH)
    d_SH <- row_min(dist_SH)
  }

  # Compute Adversarial Accuracy components
  # AA(T,S) = 0.5 * (share(d_TS > d_TT) + share(d_ST > d_SS))
  aa_train_left <- mean(d_TS > d_TT)
  aa_train_right <- mean(d_ST > d_SS)
  aa_train <- 0.5 * (aa_train_left + aa_train_right)

  # AA(H,S) = 0.5 * (share(d_HS > d_HH) + share(d_SH > d_SS))
  aa_holdout_left <- mean(d_HS > d_HH)
  aa_holdout_right <- mean(d_SH > d_SS)
  aa_holdout <- 0.5 * (aa_holdout_left + aa_holdout_right)

  # Privacy Loss = AA(holdout, S) - AA(train, S)
  privacy_loss <- aa_holdout - aa_train

  # Privacy pass: common threshold is 0.1
  privacy_pass <- privacy_loss <= 0.1

  results <- list(
    aa_train = aa_train,
    aa_holdout = aa_holdout,
    privacy_loss = privacy_loss,
    privacy_pass = privacy_pass,
    aa_train_left = aa_train_left,
    aa_train_right = aa_train_right,
    aa_holdout_left = aa_holdout_left,
    aa_holdout_right = aa_holdout_right,
    d_TS = d_TS,
    d_TT = d_TT,
    d_ST = d_ST,
    d_SS = d_SS,
    d_HS = d_HS,
    d_HH = d_HH,
    d_SH = d_SH,
    n_train = n_train,
    n_synthetic = n_synthetic,
    n_holdout = n_holdout,
    method = method,
    vars = vars,
    holdout_fraction = if (is.null(holdout)) holdout_fraction else NA
  )

  class(results) <- "nnaa"
  return(results)
}

#' Print method for nnaa objects
#'
#' @param x an object of class "nnaa"
#' @param ... additional arguments (ignored)
#' @export
print.nnaa <- function(x, ...) {
  cat("Nearest-Neighbor Adversarial Accuracy (NNAA)\n")
  cat("=============================================\n")
  cat("Method:", x$method, "\n")
  cat("Variables used:", length(x$vars), "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Training records:", x$n_train, "\n")
  cat("  Holdout records:", x$n_holdout, "\n")
  cat("  Synthetic records:", x$n_synthetic, "\n\n")

  cat("Adversarial Accuracy (ideal ~ 0.5):\n")
  cat("  AA (train vs synth):", round(x$aa_train, 4), "\n")
  cat("    Left  (real NN > synth NN):", round(x$aa_train_left, 4), "\n")
  cat("    Right (synth NN > real NN):", round(x$aa_train_right, 4), "\n")

  if (!is.na(x$aa_holdout)) {
    cat("  AA (holdout vs synth):", round(x$aa_holdout, 4), "\n")
    cat("    Left  (hold NN > synth NN):", round(x$aa_holdout_left, 4), "\n")
    cat("    Right (synth NN > hold NN):", round(x$aa_holdout_right, 4), "\n\n")

    cat("Privacy Loss (ideal ~ 0):\n")
    cat("  Privacy Loss:", round(x$privacy_loss, 4),
        if (x$privacy_loss > 0.1) " (privacy concern)" else "", "\n\n")
  } else {
    cat("\n")
  }

  cat("Privacy Assessment:")
  if (x$privacy_pass) {
    cat(" PASS\n")
    cat("  No evidence of training data memorization.\n")
  } else {
    cat(" WARNING\n")
    cat("  Synthetic data may have memorized training records.\n")
    cat("  Privacy Loss > 0.1 indicates training data is closer to\n")
    cat("  synthetic than unseen holdout data.\n")
  }

  invisible(x)
}

#' Summary method for nnaa objects
#'
#' @param object an object of class "nnaa"
#' @param ... additional arguments (ignored)
#' @export
summary.nnaa <- function(object, ...) {
  summ <- list(
    aa_train = object$aa_train,
    aa_holdout = object$aa_holdout,
    privacy_loss = object$privacy_loss,
    privacy_pass = object$privacy_pass,
    aa_train_left = object$aa_train_left,
    aa_train_right = object$aa_train_right,
    aa_holdout_left = object$aa_holdout_left,
    aa_holdout_right = object$aa_holdout_right,
    quantiles_d_TS = quantile(object$d_TS,
                              probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                              na.rm = TRUE),
    quantiles_d_TT = quantile(object$d_TT,
                              probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                              na.rm = TRUE),
    quantiles_d_SS = quantile(object$d_SS,
                              probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                              na.rm = TRUE),
    mean_d_TS = mean(object$d_TS, na.rm = TRUE),
    mean_d_TT = mean(object$d_TT, na.rm = TRUE),
    mean_d_ST = mean(object$d_ST, na.rm = TRUE),
    mean_d_SS = mean(object$d_SS, na.rm = TRUE),
    n_train = object$n_train,
    n_synthetic = object$n_synthetic,
    n_holdout = object$n_holdout,
    method = object$method,
    vars = object$vars
  )

  class(summ) <- "summary.nnaa"
  return(summ)
}

#' Print method for summary.nnaa objects
#'
#' @param x an object of class "summary.nnaa"
#' @param ... additional arguments (ignored)
#' @export
print.summary.nnaa <- function(x, ...) {
  cat("Summary: Nearest-Neighbor Adversarial Accuracy (NNAA)\n")
  cat("=====================================================\n")
  cat("Method:", x$method, "\n")
  cat("Variables:", length(x$vars), "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Training:", x$n_train, "| Holdout:", x$n_holdout,
      "| Synthetic:", x$n_synthetic, "\n\n")

  cat("Key Metrics:\n")
  cat("  AA (train):", round(x$aa_train, 4), "(ideal: ~0.5)\n")
  cat("  AA (holdout):", round(x$aa_holdout, 4), "(ideal: ~0.5)\n")
  cat("  Privacy Loss:", round(x$privacy_loss, 4), "(ideal: ~0.0)\n")
  cat("  Privacy:", ifelse(x$privacy_pass, "PASS", "WARNING"), "\n\n")

  cat("AA Components:\n")
  cat("  Train:   left =", round(x$aa_train_left, 4),
      " | right =", round(x$aa_train_right, 4), "\n")
  cat("  Holdout: left =", round(x$aa_holdout_left, 4),
      " | right =", round(x$aa_holdout_right, 4), "\n\n")

  cat("NN Distances (Training -> Synthetic):\n")
  cat("  Mean:", round(x$mean_d_TS, 4), "\n")
  cat("  Quantiles:\n")
  print(round(x$quantiles_d_TS, 4))
  cat("\n")

  cat("NN Distances (Training -> Training, self-excluded):\n")
  cat("  Mean:", round(x$mean_d_TT, 4), "\n")
  cat("  Quantiles:\n")
  print(round(x$quantiles_d_TT, 4))
  cat("\n")

  cat("NN Distances (Synthetic -> Synthetic, self-excluded):\n")
  cat("  Mean:", round(x$mean_d_SS, 4), "\n")
  cat("  Quantiles:\n")
  print(round(x$quantiles_d_SS, 4))

  invisible(x)
}

#' Plot method for nnaa objects
#'
#' @param x an object of class "nnaa"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot(s) to show: 1 = AA comparison barplot (train vs holdout),
#'   2 = NN distance distributions (overlay histograms of d_TS vs d_TT)
#' @export
plot.nnaa <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # AA comparison barplot
    vals <- c(x$aa_train, x$aa_holdout)
    names_vals <- c("AA (train)", "AA (holdout)")
    cols <- c("coral", "steelblue")

    bp <- barplot(vals, names.arg = names_vals,
                  main = paste("Adversarial Accuracy\nPrivacy Loss:",
                               round(x$privacy_loss, 4)),
                  ylab = "Adversarial Accuracy",
                  col = cols, ylim = c(0, 1), ...)
    abline(h = 0.5, col = "grey40", lwd = 2, lty = 2)
    text(bp, vals + 0.03, labels = round(vals, 3), cex = 0.9)
  }

  if (show[2]) {
    # NN distance distributions: d_TS vs d_TT
    all_dists <- c(x$d_TS, x$d_TT)
    breaks <- seq(0, max(all_dists, na.rm = TRUE) * 1.1, length.out = 31)

    hist(x$d_TT, breaks = breaks, col = rgb(0, 0, 1, 0.5),
         main = "NN Distance Distributions",
         xlab = "Nearest-Neighbor Distance", ...)
    hist(x$d_TS, breaks = breaks, col = rgb(1, 0, 0, 0.5), add = TRUE)
    legend("topright",
           legend = c(
             paste("Train->Train (mean:", round(mean(x$d_TT), 3), ")"),
             paste("Train->Synth (mean:", round(mean(x$d_TS), 3), ")")
           ),
           fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)),
           cex = 0.8)
  }
}

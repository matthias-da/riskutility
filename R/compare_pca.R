#' Compare Principal Component Analysis (PCA) between Two Datasets with Separate Loadings
#'
#' This function performs Principal Component Analysis (PCA) on two datasets (X and Y) for a specified set of numeric variables.
#' When side_by_side is FALSE, a combined PCA is performed on the union of the data for visualization of the point scores,
#' while separate PCA analyses are computed for X and Y to obtain loadings. The function then overlays loadings arrows from X (blue)
#' and Y (red) onto the combined PCA biplot. When side_by_side is TRUE, separate biplots for X and Y are produced and arranged side by side.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the anonymized/synthetic dataset.
#' @param vars A character vector specifying the numeric variables to include in the PCA.
#' @param center Logical. Should the variables be centered (zero mean)? Default is TRUE.
#' @param scale Logical. Should the variables be scaled to unit variance? Default is TRUE.
#' @param biplot Logical. If TRUE, loadings are added as arrows. Default is TRUE.
#' @param side_by_side Logical. If TRUE, separate PCA biplots are produced for X and Y and arranged side by side.
#'        If FALSE (default), a combined PCA biplot is produced with separate loadings for X and Y.
#'
#' @param ... additional arguments passed to methods
#'
#' @return A list with two elements:
#' \describe{
#'   \item{pca}{If side_by_side = FALSE, a list with three PCA objects: combined, X, and Y.
#'              If side_by_side = TRUE, a list with two PCA objects: X and Y.}
#'   \item{plot}{A ggplot2 object: either a combined biplot with separate loadings or separate plots arranged side by side.}
#' }
#'
#' @details When side_by_side is FALSE, the function first adds a dataset indicator to X and Y, combines them,
#' and performs PCA on the combined dataset for the point scores. Then, PCA is computed separately for X and Y.
#' Loadings from these separate analyses are scaled using the range of the combined PCA scores and overlaid on the
#' combined biplot with blue arrows for X and red arrows for Y.
#'
#' @importFrom data.table as.data.table rbindlist
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_text labs theme_minimal scale_color_manual
#' @family comparison
#' @author Matthias Templ
#' @export
#'
#' @examples
#' set.seed(123)
#' X <- data.frame(
#'   income = rnorm(500, mean = 50000, sd = 10000),
#'   age = rnorm(500, mean = 40, sd = 10)
#' )
#' # Add a positively correlated variable
#' X$expenses <- X$income * 0.5 + rnorm(n = 500, mean = 1000, sd = 500)
#' Y <- data.frame(
#'   income = rnorm(1000, mean = 48000, sd = 12000),
#'   age = rnorm(1000, mean = 42, sd = 11)
#' )
#' Y$expenses <- Y$income * 0.5 + rnorm(n = 1000, mean = 1000, sd = 500)
#' # Combined PCA biplot with separate loadings for X (blue) and Y (red)
#' res_combined <- compare_pca(X, Y,
#'   vars = c("income", "age", "expenses"),
#'   biplot = TRUE, side_by_side = FALSE)
#' print(res_combined$pca)
#' print(res_combined$plot)
#'
#' # Separate PCA biplots for X and Y arranged side by side
#' res_separate <- compare_pca(X, Y, vars = c("income", "age"), biplot = TRUE, side_by_side = TRUE)
#' print(res_separate$pca)
#' print(res_separate$plot)
compare_pca <- function(X, ...) {
  UseMethod("compare_pca")
}

#' @rdname compare_pca
#' @export
compare_pca.synth_pair <- function(X, ...) {
  compare_pca.default(X = X$original, Y = X$synthetic, ...)
}

#' @rdname compare_pca
#' @export
compare_pca.default <- function(X, Y, vars, center = TRUE, scale = TRUE, biplot = TRUE, side_by_side = FALSE, ...) {
  create_biplot_arrows <- function(pca_result, scores_df) {
    arrow_multiplier <- min(
      (max(scores_df$PC1) - min(scores_df$PC1)) / (max(pca_result$rotation[, 1]) - min(pca_result$rotation[, 1])),
      (max(scores_df$PC2) - min(scores_df$PC2)) / (max(pca_result$rotation[, 2]) - min(pca_result$rotation[, 2]))
    )
    loadings <- as.data.frame(pca_result$rotation[, 1:2])
    loadings$var <- rownames(loadings)
    loadings$PC1 <- loadings[, 1] * arrow_multiplier
    loadings$PC2 <- loadings[, 2] * arrow_multiplier
    return(loadings)
  }

  if (!side_by_side) {
    X_dt <- as.data.table(X)[, ..vars]
    X_dt[, .dataset := "Original (X)"]
    Y_dt <- as.data.table(Y)[, ..vars]
    Y_dt[, .dataset := "Anonymized/Synthetic (Y)"]

    combined <- rbindlist(list(X_dt, Y_dt), use.names = TRUE, fill = TRUE)
    combined$.dataset <- factor(combined$.dataset, levels = c("Original (X)", "Anonymized/Synthetic (Y)"))

    pca_combined <- prcomp(combined[, ..vars], center = center, scale. = scale)
    scores_combined <- as.data.frame(pca_combined$x)
    scores_combined$.dataset <- combined$.dataset

    pca_X <- prcomp(X_dt[, ..vars], center = center, scale. = scale)
    pca_Y <- prcomp(Y_dt[, ..vars], center = center, scale. = scale)

    loadings_X <- create_biplot_arrows(pca_X, scores_combined)
    loadings_Y <- create_biplot_arrows(pca_Y, scores_combined)

    p <- ggplot(scores_combined, aes(x = PC1, y = PC2, color = .dataset)) +
      geom_point(alpha = 0.7) +
      labs(x = "PC1", y = "PC2", title = "Combined PCA Biplot with Separate Loadings") +
      theme_minimal() +
      scale_color_manual(values = c("blue", "red"))

    if (biplot) {
      p <- p +
        geom_segment(data = loadings_X, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                     arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
        geom_text(data = loadings_X, aes(x = PC1, y = PC2, label = var), color = "blue", vjust = 1.5, size = 3) +
        geom_segment(data = loadings_Y, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                     arrow = arrow(length = unit(0.2, "cm")), color = "red") +
        geom_text(data = loadings_Y, aes(x = PC1, y = PC2, label = var), color = "red", vjust = -0.5, size = 3)
    }

    return(list(pca = list(combined = pca_combined, X = pca_X, Y = pca_Y), plot = p))

  } else {
    X_dt <- as.data.table(X)[, ..vars]
    Y_dt <- as.data.table(Y)[, ..vars]

    pca_X <- prcomp(X_dt, center = center, scale. = scale)
    pca_Y <- prcomp(Y_dt, center = center, scale. = scale)

    scores_X <- as.data.frame(pca_X$x)
    scores_Y <- as.data.frame(pca_Y$x)

    pX <- ggplot(scores_X, aes(x = PC1, y = PC2)) +
      geom_point(alpha = 0.7, color = "blue") +
      labs(x = "PC1", y = "PC2", title = "PCA Biplot: Original (X)") +
      theme_minimal()

    pY <- ggplot(scores_Y, aes(x = PC1, y = PC2)) +
      geom_point(alpha = 0.7, color = "red") +
      labs(x = "PC1", y = "PC2", title = "PCA Biplot: Anonymized/Synthetic (Y)") +
      theme_minimal()

    if (biplot) {
      loadings_X <- create_biplot_arrows(pca_X, scores_X)
      loadings_Y <- create_biplot_arrows(pca_Y, scores_Y)

      pX <- pX +
        geom_segment(data = loadings_X, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                     arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
        geom_text(data = loadings_X, aes(x = PC1, y = PC2, label = var), color = "blue", vjust = 1.5, size = 3)

      pY <- pY +
        geom_segment(data = loadings_Y, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                     arrow = arrow(length = unit(0.2, "cm")), color = "red") +
        geom_text(data = loadings_Y, aes(x = PC1, y = PC2, label = var), color = "red", vjust = -0.5, size = 3)
    }

    if (!requireNamespace("gridExtra", quietly = TRUE)) {
      stop("Package 'gridExtra' is required for side_by_side plots. Please install it.")
    }
    combined_plot <- gridExtra::grid.arrange(pX, pY, ncol = 2)

    return(list(pca = list(X = pca_X, Y = pca_Y), plot = combined_plot))
  }
}


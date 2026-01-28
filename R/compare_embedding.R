#' Compare Dimensionality Reduction Methods (t-SNE, UMAP, MDS-Sammon) for Two Datasets
#'
#' This function compares two datasets using dimensionality reduction methods: t-distributed
#' Stochastic Neighbor Embedding (t-SNE), Uniform Manifold Approximation and Projection (UMAP),
#' and Multidimensional Scaling (MDS) with Sammon's mapping. It allows visualizing complex
#' data structures, such as clusters or outliers, in two-dimensional space.
#'
#' @param X,Y data.frame or data.table: datasets to compare.
#' @param vars Character vector: numeric variables for dimensionality reduction.
#' @param method Character: one of 'tsne', 'umap', or 'mds'.
#' @param perplexity Numeric: t-SNE perplexity (default: 30).
#' @param n_neighbors Numeric: UMAP neighbors (default: 15).
#' @param side_by_side Logical: if TRUE, plots embeddings for X and Y side by side.
#' @param seed Integer: random seed for reproducibility of stochastic methods (t-SNE, UMAP).
#'   Default NULL uses current random state.
#'
#' @return A list containing the embeddings and a ggplot2 visualization.
#' @importFrom MASS sammon
#' @import ggplot2
#' @import data.table
#' @export
#'
#' @examples
#' # MDS (Sammon) embedding - doesn't require optional packages
#' set.seed(123)
#' X <- data.frame(
#'   var1 = rnorm(100, 50, 10),
#'   var2 = rnorm(100, 40, 8),
#'   var3 = rnorm(100, 30, 5)
#' )
#' Y <- data.frame(
#'   var1 = rnorm(100, 48, 12),
#'   var2 = rnorm(100, 42, 9),
#'   var3 = rnorm(100, 28, 6)
#' )
#' res_mds <- compare_embedding(X, Y, vars = c("var1", "var2", "var3"),
#'                              method = 'mds')
#'
#' \donttest{
#' # t-SNE (requires Rtsne package)
#' if (requireNamespace("Rtsne", quietly = TRUE)) {
#'   res_tsne <- compare_embedding(X, Y, vars = c("var1", "var2", "var3"),
#'                                 method = 'tsne')
#' }
#'
#' # UMAP (requires uwot package)
#' if (requireNamespace("uwot", quietly = TRUE)) {
#'   res_umap <- compare_embedding(X, Y, vars = c("var1", "var2", "var3"),
#'                                 method = 'umap')
#' }
#' }
compare_embedding <- function(X, Y, vars, method = 'tsne', perplexity = 30, n_neighbors = 15, side_by_side = FALSE, seed = NULL) {
  # Set seed for reproducibility of stochastic methods
  if (!is.null(seed)) set.seed(seed)

  X <- as.data.table(X)[, ..vars]
  Y <- as.data.table(Y)[, ..vars]

  get_embedding <- function(data, method) {
    if (method == 'tsne') {
      if (!requireNamespace("Rtsne", quietly = TRUE)) {
        stop("Package 'Rtsne' is required for this function. Please install it.")
      }
      emb <- Rtsne::Rtsne(data, perplexity = perplexity)$Y
    } else if (method == 'umap') {
      if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("Package 'uwot' is required for this function. Please install it.")
      }
      emb <- uwot::umap(data, n_neighbors = n_neighbors)
    } else if (method == 'mds') {
      dist_matrix <- dist(scale(data))
      emb <- MASS::sammon(dist_matrix)$points
    }
    return(as.data.table(emb)[, .(Dim1 = V1, Dim2 = V2)])
  }

  emb_X <- get_embedding(X, method)
  emb_Y <- get_embedding(Y, method)

  emb_X[, dataset := 'Original (X)']
  emb_Y[, dataset := 'Anonymized/Synthetic (Y)']

  combined_emb <- rbind(emb_X, emb_Y)

  plot_embedding <- function(data, title) {
    ggplot(data, aes(Dim1, Dim2, color = dataset)) +
      geom_point(alpha = 0.7) +
      theme_minimal() +
      labs(title = title, color = "Dataset") +
      scale_color_manual(values = c("blue", "red"))
  }

  if (!side_by_side) {
    plot <- plot_embedding(combined_emb, paste(toupper(method), "Embedding (Combined)"))
  } else {
    plot_X <- plot_embedding(emb_X, paste(toupper(method), "Embedding: Original (X)"))
    plot_Y <- plot_embedding(emb_Y, paste(toupper(method), "Embedding: Synthetic (Y)"))
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
      stop("Package 'gridExtra' is required for this function. Please install it.")
    }
    plot <- gridExtra::grid.arrange(plot_X, plot_Y, ncol = 2)
  }

  list(embedding = list(X = emb_X, Y = emb_Y), plot = plot)
}

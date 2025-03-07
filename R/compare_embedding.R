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
#'
#' @return A list containing the embeddings and a ggplot2 visualization.
#' @importFrom Rtsne Rtsne
#' @importFrom uwot umap
#' @importFrom MASS sammon
#' @import ggplot2
#' @import data.table
#' @import gridExtra
#' @export
#'
#' @examples
#' library(data.table)
#' set.seed(123)
#'
#' vars <- c("income", "age", "expenses", "savings", "debt")
#'
#' mu_X <- c(50000, 40, 30000, 20000, 10000)
#' mu_Y <- c(48000, 42, 28000, 22000, 12000)
#'
#' Sigma <- matrix(c(
#'   1.0, 0.6, 0.4, -0.3, 0.2,
#'   0.6, 1.0, 0.5, -0.2, 0.3,
#'   0.4, 0.5, 1.0, -0.4, 0.3,
#'  -0.3, -0.2, -0.4, 1.0, -0.6,
#'   0.2, 0.3, 0.3, -0.6, 1.0
#' ), byrow = TRUE, ncol = 5) * 10000
#'
#' cluster_X1 <- MASS::mvrnorm(250, mu_X, Sigma)
#' cluster_X2 <- MASS::mvrnorm(250, mu_X + 10000, Sigma)
#' cluster_Y1 <- MASS::mvrnorm(250, mu_Y, Sigma)
#' cluster_Y2 <- MASS::mvrnorm(250, mu_Y - 10000, Sigma)
#'
#' X <- as.data.table(rbind(cluster_X1, cluster_X2))
#' setnames(X, vars)
#'
#' Y <- as.data.table(rbind(cluster_Y1, cluster_Y2))
#' setnames(Y, vars)
#'
#' # t-SNE
#' res_tsne <- compare_embedding(X, Y, vars, method = 'tsne', side_by_side = TRUE)
#' print(res_tsne$plot)
#' res_tsne <- compare_embedding(X, Y, vars, method = 'tsne', side_by_side = FALSE)
#' print(res_tsne$plot)
#'
#' # UMAP
#' res_umap <- compare_embedding(X, Y, vars, method = 'umap', side_by_side = FALSE)
#' print(res_umap$plot)
#'
#' # MDS (Sammon)
#' res_mds <- compare_embedding(X, Y, vars, method = 'mds', side_by_side = FALSE)
#' print(res_mds$plot)
compare_embedding <- function(X, Y, vars, method = 'tsne', perplexity = 30, n_neighbors = 15, side_by_side = FALSE) {
  X <- as.data.table(X)[, ..vars]
  Y <- as.data.table(Y)[, ..vars]

  get_embedding <- function(data, method) {
    if (method == 'tsne') {
      emb <- Rtsne::Rtsne(data, perplexity = perplexity)$Y
    } else if (method == 'umap') {
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
    plot <- gridExtra::grid.arrange(plot_X, plot_Y, ncol = 2)
  }

  list(embedding = list(X = emb_X, Y = emb_Y), plot = plot)
}

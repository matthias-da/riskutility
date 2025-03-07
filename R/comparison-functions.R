#' Tools for Comparing Original and Synthetic/Anonymized Datasets
#'
#' Provides a comprehensive set of tools to assess similarity between an original dataset (X)
#' and its synthetic or anonymized counterpart (Y). Includes visualization, statistical tests,
#' dimensionality reduction, multivariate analyses, feature importance, and predictive model comparisons.
#'
#' @section Visualization Functions:
#' \describe{
#'   \item{\code{compare_histograms()}}{Compare distributions using histograms (weighted and unweighted, faceted).}
#'   \item{\code{compare_boxplots()}}{Compare numeric variable distributions using weighted/unweighted boxplots.}
#' }
#'
#' @section Statistical Tests:
#' \describe{
#'   \item{\code{compare_ks_test()}}{Kolmogorov-Smirnov test for comparing distributions.}
#'   \item{\code{compare_chisq_gof()}}{Chi-square Goodness-of-Fit test with multi-dimensional frequency tables (weighted/unweighted).}
#' }
#'
#' @section Distributional Comparisons:
#' \describe{
#'   \item{\code{compare_wasserstein()}}{Calculate Wasserstein distances between numeric or nominal variables.}
#'   \item{\code{compare_multivariate_distribution()}}{Compare joint distributions using Mahalanobis distance or mutual information.}
#' }
#'
#' @section Summary Statistics:
#' \describe{
#'   \item{\code{compare_means_frequencies()}}{Compare means, medians, robust statistics, skewness, kurtosis for numeric data, and frequencies for categorical data.}
#'   \item{\code{compare_correlation()}}{Compare correlations: Pearson, Spearman, robust correlations for numeric, categorical, and mixed data types.}
#'   \item{\code{multivariate_summary()}}{Calculate joint summary statistics for numeric and categorical data.}
#' }
#'
#' @section Dimensionality Reduction:
#' \describe{
#'   \item{\code{compare_pca()}}{Principal Component Analysis comparison with biplot options.}
#'   \item{\code{compare_embedding()}}{t-SNE, UMAP, or Sammon's Mapping (MDS) visualizations for embedding comparisons.}
#' }
#'
#' @section Machine Learning and Predictive Modeling:
#' \describe{
#'   \item{\code{compare_model_performance()}}{Evaluate predictive model performance (accuracy, precision, recall, F1-score, AUC, R-squared, RMSE) using cross-validation.}
#'   \item{\code{compare_feature_importance()}}{Evaluate stability of feature importance measures from Random Forests, Decision Trees, Gradient Boosting, Permutation Importance, and SHAP values.}
#' }
#'
#' @section Dependencies:
#' The package leverages several R packages including \code{ggplot2}, \code{data.table}, \code{caret}, \code{vip}, \code{simPop}, \code{survey}, \code{Rtsne}, \code{uwot}, \code{MASS}, \code{robustbase}, and \code{psych}.
#'
#' @author
#' Your Name <your.email@example.com>
#'
#' @seealso
#' \code{\link[ggplot2]{ggplot}}, \code{\link[data.table]{data.table}}, \code{\link[caret]{train}}, \code{\link[vip]{vi_permute}}
#'
#' @examples
#' # Example for histogram comparison
#' compare_histograms(X, Y, num_var = "income", cat_vars = c("gender"),
#'                    weight_X = "weight", weight_Y = "weight")
#'
#' # Example for PCA comparison
#' set.seed(123)
#' X <- data.frame(
#'   income = rnorm(500, mean = 50000, sd = 10000),
#'   age = rnorm(500, mean = 40, sd = 10),
#'   expenditure = rnorm(500, mean = 2000, sd = 500)
#' )
#' Y <- data.frame(
#'   income = rnorm(1000, mean = 48000, sd = 12000),
#'   age = rnorm(1000, mean = 42, sd = 11),
#'   expenditure = rnorm(1000, mean = 2200, sd = 600)
#' )
#' res_pca <- compare_pca(X, Y, vars = c("income", "age", "expenditure"), biplot = TRUE)
#' print(res_pca$plot)
#'
#' # Example for embeddings with realistic clusters
#' library(MASS)
#' Sigma <- matrix(c(1,0.6,0.3,0.2,0.1,
#'                   0.6,1,0.4,0.3,0.2,
#'                   0.3,0.4,1,0.5,0.3,
#'                   0.2,0.3,0.5,1,0.4,
#'                   0.1,0.2,0.3,0.4,1), ncol=5)
#' cluster1 <- MASS::mvrnorm(250, c(50,40,2000,5,3), Sigma*1000)
#' cluster2 <- MASS::mvrnorm(250, c(70,60,3000,8,6), Sigma*1000)
#' X <- data.table(rbind(cluster1, cluster2))
#' names(X) <- c("income","age","expenditure","score1","score2")
#'
#' cluster1_y <- MASS::mvrnorm(500, c(48,38,2100,4.5,2.5), Sigma*1200)
#' cluster2_y <- MASS::mvrnorm(500, c(68,58,2900,7.5,5.5), Sigma*1200)
#' Y <- data.table(rbind(cluster1_y, cluster2_y))
#' names(Y) <- c("income","age","expenditure","score1","score2")
#'
#' res_embed <- compare_embedding(X, Y, vars = names(X), method = 'umap')
#' print(res_embed$plot)
#'
#' @docType package
#' @name ComparisonTools
NULL

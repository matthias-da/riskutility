#' Compare Correlation Matrices between Two Datasets
#'
#' This function computes and compares the correlation matrices for a specified set of variables between two
#' datasets (X and Y), such as an original dataset and an anonymized/synthetic version. For continuous variables,
#' the user can choose Pearson, Spearman, or robust correlation methods (using either MM‐estimation or the MCD estimator).
#' For categorical variables, the function computes Cramér’s V. For variable pairs with one continuous and one nominal variable,
#' if the nominal variable is binary, the point‐biserial correlation is computed; otherwise, the correlation ratio (eta squared)
#' is used and transformed into a correlation‐like measure.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the anonymized/synthetic dataset.
#' @param vars A character vector specifying the variables to include in the correlation analysis.
#' @param method A character string specifying the correlation method for continuous variables. Options are
#'   \code{"pearson"}, \code{"spearman"}, \code{"robust_mcd"}, or \code{"robust_mm"}. For categorical variable pairs, the function
#'   uses Cramér's V.
#' @param mixed_method A character string specifying the method for continuous vs nominal pairs. Options are
#'   \code{"point_biserial"} (if the nominal variable is binary) or \code{"eta_squared"} (for nominal variables with more than two levels).
#'   The default is \code{"point_biserial"}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{corr_X}{The correlation matrix computed from dataset X.}
#'   \item{corr_Y}{The correlation matrix computed from dataset Y.}
#'   \item{diff}{The absolute difference matrix between the two correlation matrices.}
#' }
#'
#' @details The function determines the type of each variable (continuous or categorical) based on its class.
#' Continuous variables are correlated using the specified method. For robust methods, the MCD estimator from the
#' \code{robustbase} package or MM-estimation via \code{MASS::cov.rob} is used. For categorical variables, Cramér’s V
#' is computed from the chi-square statistic. For mixed pairs, if one variable is continuous and the other nominal,
#' point-biserial correlation is used when the nominal variable has two levels, and the correlation ratio (eta squared)
#' is computed otherwise.
#'
#' @importFrom data.table as.data.table
#' @importFrom robustbase covMcd
#' @importFrom MASS cov.rob
#' @export
#'
#' @examples
#' set.seed(123)
#' X <- data.frame(
#'   income = rnorm(500, mean = 50000, sd = 10000),
#'   age = rnorm(500, mean = 40, sd = 10),
#'   gender = factor(sample(c("Male", "Female"), 500, replace = TRUE)),
#'   region = factor(sample(c("North", "South", "East", "West"), 500, replace = TRUE))
#' )
#'
#' Y <- data.frame(
#'   income = rnorm(1000, mean = 48000, sd = 12000),
#'   age = rnorm(1000, mean = 42, sd = 11),
#'   gender = factor(sample(c("Male", "Female"), 1000, replace = TRUE)),
#'   region = factor(sample(c("North", "South", "East", "West"), 1000, replace = TRUE))
#' )
#'
#' # Continuous correlation using Pearson
#' res1 <- compare_correlation_matrices(X, Y, vars = c("income", "age"), method = "pearson")
#' print(res1$corr_X)
#' print(res1$corr_Y)
#' print(res1$diff)
#'
#' # Categorical correlation using Cramér's V (computed automatically)
#' res2 <- compare_correlation_matrices(X, Y, vars = c("gender", "region"), method = "pearson")
#' print(res2$corr_X)
#' print(res2$corr_Y)
#'
#' # Mixed example: income (continuous) vs gender (nominal); point-biserial or eta-squared is computed.
#' res3 <- compare_correlation_matrices(X, Y, vars = c("income", "gender"), method = "pearson", mixed_method = "point_biserial")
#' print(res3$corr_X)
#' print(res3$corr_Y)
compare_correlation_matrices <- function(X, Y, vars, method = "pearson", mixed_method = "point_biserial") {

  library(data.table)
  library(robustbase)
  library(MASS)

  # Convert X and Y to data.table if needed
  X <- as.data.table(X)
  Y <- as.data.table(Y)

  # Determine variable types for each variable in vars
  var_types <- sapply(X[, vars, with = FALSE], function(v) {
    if (is.numeric(v)) "continuous" else "categorical"
  })

  # Initialize correlation matrices
  corr_X <- matrix(NA, nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
  corr_Y <- matrix(NA, nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))

  # Helper functions
  # Continuous correlations
  compute_continuous_corr <- function(x, y, method) {
    cor(x, y, method = method, use = "pairwise.complete.obs")
  }

  compute_robust_mcd_corr <- function(x, y) {
    cov_mcd <- covMcd(cbind(x, y))
    cor(cov_mcd$cov)[1, 2]
  }

  compute_robust_mm_corr <- function(x, y) {
    cov_mm <- cov.rob(cbind(x, y), method = "MM")
    cor(cov_mm$cov)[1, 2]
  }

  # Categorical correlation (Cramér's V)
  compute_cramers_v <- function(x, y) {
    tbl <- table(x, y)
    chi2 <- suppressWarnings(chisq.test(tbl)$statistic)
    n <- sum(tbl)
    min_dim <- min(nrow(tbl) - 1, ncol(tbl) - 1)
    sqrt(as.numeric(chi2) / (n * min_dim))
  }

  # Mixed correlation: point-biserial correlation
  compute_point_biserial <- function(x, binary) {
    binary <- as.numeric(binary) - min(as.numeric(binary))
    cor(x, binary, method = "pearson", use = "pairwise.complete.obs")
  }

  # Mixed correlation: correlation ratio (eta squared converted to correlation)
  compute_eta_ratio <- function(x, factor_var) {
    aov_res <- aov(x ~ factor(factor_var))
    ss_between <- sum((tapply(x, factor(factor_var), mean) - mean(x))^2 * table(factor(factor_var)))
    ss_total <- sum((x - mean(x))^2)
    eta2 <- ss_between / ss_total
    sqrt(eta2)
  }

  # Loop over pairs of variables
  for (i in seq_along(vars)) {
    for (j in seq_along(vars)) {
      if (i <= j) {
        # Determine the types for each variable
        type_i <- var_types[i]
        type_j <- var_types[j]
        x_i <- X[[vars[i]]]
        x_j <- X[[vars[j]]]
        y_i <- Y[[vars[i]]]
        y_j <- Y[[vars[j]]]

        # Compute correlation for dataset X
        corr_val_X <- NA
        if (type_i == "continuous" && type_j == "continuous") {
          if (method %in% c("pearson", "spearman")) {
            corr_val_X <- compute_continuous_corr(x_i, x_j, method)
          } else if (method == "robust_mcd") {
            corr_val_X <- compute_robust_mcd_corr(x_i, x_j)
          } else if (method == "robust_mm") {
            corr_val_X <- compute_robust_mm_corr(x_i, x_j)
          }
        } else if (type_i == "categorical" && type_j == "categorical") {
          corr_val_X <- compute_cramers_v(x_i, x_j)
        } else if ( (type_i == "continuous" && type_j == "categorical") ||
                    (type_i == "categorical" && type_j == "continuous") ) {
          # Identify which is continuous and which is categorical
          if (type_i == "continuous") {
            if (length(unique(x_j)) == 2 && mixed_method == "point_biserial") {
              corr_val_X <- compute_point_biserial(x_i, x_j)
            } else {
              corr_val_X <- compute_eta_ratio(x_i, x_j)
            }
          } else {
            if (length(unique(x_i)) == 2 && mixed_method == "point_biserial") {
              corr_val_X <- compute_point_biserial(x_j, x_i)
            } else {
              corr_val_X <- compute_eta_ratio(x_j, x_i)
            }
          }
        }
        corr_X[i, j] <- corr_val_X
        corr_X[j, i] <- corr_val_X

        # Compute correlation for dataset Y similarly
        corr_val_Y <- NA
        if (type_i == "continuous" && type_j == "continuous") {
          if (method %in% c("pearson", "spearman")) {
            corr_val_Y <- compute_continuous_corr(y_i, y_j, method)
          } else if (method == "robust_mcd") {
            corr_val_Y <- compute_robust_mcd_corr(y_i, y_j)
          } else if (method == "robust_mm") {
            corr_val_Y <- compute_robust_mm_corr(y_i, y_j)
          }
        } else if (type_i == "categorical" && type_j == "categorical") {
          corr_val_Y <- compute_cramers_v(y_i, y_j)
        } else if ( (type_i == "continuous" && type_j == "categorical") ||
                    (type_i == "categorical" && type_j == "continuous") ) {
          if (type_i == "continuous") {
            if (length(unique(y_j)) == 2 && mixed_method == "point_biserial") {
              corr_val_Y <- compute_point_biserial(y_i, y_j)
            } else {
              corr_val_Y <- compute_eta_ratio(y_i, y_j)
            }
          } else {
            if (length(unique(y_i)) == 2 && mixed_method == "point_biserial") {
              corr_val_Y <- compute_point_biserial(y_j, y_i)
            } else {
              corr_val_Y <- compute_eta_ratio(y_j, y_i)
            }
          }
        }
        corr_Y[i, j] <- corr_val_Y
        corr_Y[j, i] <- corr_val_Y
      }
    }
  }

  diff <- abs(corr_X - corr_Y)

  return(list(corr_X = corr_X, corr_Y = corr_Y, diff = diff))
}

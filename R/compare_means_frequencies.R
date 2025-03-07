#' Compare Means and Frequencies between Two Datasets
#'
#' This function compares the central tendencies of continuous variables and the frequency
#' distributions of categorical variables between two datasets (e.g., an original dataset and
#' an anonymized/synthetic dataset). For continuous variables, the function computes a set of
#' summary statistics for each variable, including arithmetic mean, median, Huber mean, standard deviation,
#' interquartile range (IQR), median absolute deviation (MAD), skewness, and kurtosis. Weighted estimates
#' are computed if sampling weights are provided. In addition, conditional summaries can be produced if grouping
#' variables are specified. For categorical variables, frequency and relative frequency tables are computed,
#' using weights if available.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the anonymized/synthetic dataset.
#' @param cont_vars A character vector specifying the continuous variables to compare.
#' @param cat_vars A character vector specifying the categorical variables for frequency comparison.
#' @param group_vars Optional. A character vector specifying the grouping variables for computing conditional summaries. Default is NULL.
#' @param weight_X Optional. A character string specifying the sampling weight variable in X.
#' @param weight_Y Optional. A character string specifying the sampling weight variable in Y.
#' @param stats A character vector of summary measures to compute for continuous variables. Options include "mean", "median", "sd", "IQR", "mad", "huber", "skewness", and "kurtosis".
#'              Default is c("mean", "median", "sd", "IQR", "mad", "huber", "skewness", "kurtosis").
#'
#' @return A list containing three elements:
#' \describe{
#'   \item{continuous_summary}{A data.table with overall summary statistics for each continuous variable in X and Y.}
#'   \item{conditional_summary}{A data.table with conditional summary statistics computed by the grouping variables (if provided); otherwise, NULL.}
#'   \item{categorical_summary}{A data.table with frequency and relative frequency tables for each categorical variable in X and Y.}
#' }
#'
#' @details For continuous variables, if sampling weights are provided the function computes weighted estimates.
#' The arithmetic mean and standard deviation are computed using \code{weighted.mean} and a weighted variance formula.
#' The median is computed using \code{Hmisc::wtd.quantile} when weights are provided.
#' The Huber mean is computed via \code{MASS::huber} (weighted Huber is not implemented). The IQR and MAD are computed
#' using base functions. For skewness and kurtosis, if weights are provided, the function computes them using custom
#' formulas; otherwise, it computes the unweighted versions.
#'
#' For categorical variables, frequencies are computed as the sum of weights (if provided) or as raw counts,
#' with relative frequencies calculated accordingly.
#'
#' @importFrom data.table as.data.table rbindlist
#' @importFrom Hmisc wtd.quantile
#' @importFrom MASS huber
#' @export
#'
#' @examples
#' set.seed(123)
#' X <- data.frame(
#'   age = sample(20:80, 500, replace = TRUE),
#'   income = rnorm(500, mean = 50000, sd = 10000),
#'   gender = sample(c("Male", "Female"), 500, replace = TRUE),
#'   region = sample(c("North", "South", "East", "West"), 500, replace = TRUE),
#'   weight = runif(500, 0.5, 1.5)
#' )
#'
#' Y <- data.frame(
#'   age = sample(20:80, 1000, replace = TRUE),
#'   income = rnorm(1000, mean = 48000, sd = 12000),
#'   gender = sample(c("Male", "Female"), 1000, replace = TRUE),
#'   region = sample(c("North", "South", "East", "West"), 1000, replace = TRUE),
#'   weight = runif(1000, 0.5, 1.5)
#' )
#'
#' result <- compare_means_frequencies(X, Y,
#'              cont_vars = c("age", "income"),
#'              cat_vars = c("gender", "region"),
#'              group_vars = c("region"),
#'              weight_X = "weight", weight_Y = "weight",
#'              stats = c("mean", "median", "sd", "IQR", "mad", "huber", "skewness", "kurtosis"))
#'
#' print(result$continuous_summary)
#' print(result$conditional_summary)
#' print(result$categorical_summary)
compare_means_frequencies <- function(X, Y, cont_vars, cat_vars, group_vars = NULL,
                                      weight_X = NULL, weight_Y = NULL,
                                      stats = c("mean", "median", "sd", "IQR", "mad", "huber", "skewness", "kurtosis")) {

  library(data.table)
  library(Hmisc)
  library(MASS)

  # Convert to data.table
  X <- as.data.table(X)
  Y <- as.data.table(Y)

  ## Helper functions for weighted skewness and kurtosis
  weighted_skewness <- function(x, w) {
    mu <- weighted.mean(x, w, na.rm = TRUE)
    sigma <- sqrt(sum(w * (x - mu)^2, na.rm = TRUE) / (sum(w) - 1))
    if(sigma == 0) return(NA)
    sum(w * (x - mu)^3, na.rm = TRUE) / (sum(w) * sigma^3)
  }

  weighted_kurtosis <- function(x, w) {
    mu <- weighted.mean(x, w, na.rm = TRUE)
    sigma <- sqrt(sum(w * (x - mu)^2, na.rm = TRUE) / (sum(w) - 1))
    if(sigma == 0) return(NA)
    sum(w * (x - mu)^4, na.rm = TRUE) / (sum(w) * sigma^4) - 3
  }

  ## Helper function: compute summary statistics for a vector, with optional weights.
  compute_stats <- function(x, w = NULL, stats) {
    if (is.null(w)) w <- rep(1, length(x))
    result <- list()
    # Arithmetic mean
    if ("mean" %in% stats) {
      result$mean <- weighted.mean(x, w, na.rm = TRUE)
    }
    # Median
    if ("median" %in% stats) {
      result$median <- as.numeric(wtd.quantile(x, weights = w, probs = 0.5, na.rm = TRUE))
    }
    # Huber mean (unweighted only)
    if ("huber" %in% stats) {
      result$huber <- huber(x)$mu
    }
    # Standard deviation
    if ("sd" %in% stats) {
      mu <- weighted.mean(x, w, na.rm = TRUE)
      result$sd <- sqrt(sum(w * (x - mu)^2, na.rm = TRUE) / (sum(w) - 1))
    }
    # IQR
    if ("IQR" %in% stats) {
      result$IQR <- IQR(x, na.rm = TRUE)
    }
    # MAD
    if ("mad" %in% stats) {
      result$mad <- mad(x, na.rm = TRUE)
    }
    # Skewness
    if ("skewness" %in% stats) {
      result$skewness <- weighted_skewness(x, w)
    }
    # Kurtosis
    if ("kurtosis" %in% stats) {
      result$kurtosis <- weighted_kurtosis(x, w)
    }
    return(result)
  }

  ## Overall continuous summaries for each variable in cont_vars for dataset X and Y
  cont_summary_list <- list()
  for (var in cont_vars) {
    # For dataset X
    if (!is.null(weight_X) && weight_X %in% names(X)) {
      stats_X <- compute_stats(X[[var]], X[[weight_X]], stats)
    } else {
      stats_X <- compute_stats(X[[var]], NULL, stats)
    }
    dt_X <- data.table(variable = var, dataset = "X", t(as.data.table(stats_X)))

    # For dataset Y
    if (!is.null(weight_Y) && weight_Y %in% names(Y)) {
      stats_Y <- compute_stats(Y[[var]], Y[[weight_Y]], stats)
    } else {
      stats_Y <- compute_stats(Y[[var]], NULL, stats)
    }
    dt_Y <- data.table(variable = var, dataset = "Y", t(as.data.table(stats_Y)))

    cont_summary_list[[var]] <- rbind(dt_X, dt_Y)
  }
  continuous_summary <- rbindlist(cont_summary_list)

  ## Conditional continuous summaries (if group_vars provided)
  if (!is.null(group_vars)) {
    cond_summary_list <- list()
    for (var in cont_vars) {
      # For dataset X conditional on group_vars
      if (!is.null(weight_X) && weight_X %in% names(X)) {
        cond_X <- X[, .(temp = list(compute_stats(get(var), get(weight_X), stats))), by = group_vars]
      } else {
        cond_X <- X[, .(temp = list(compute_stats(get(var), NULL, stats))), by = group_vars]
      }
      cond_X[, dataset := "X"]
      cond_X[, variable := var]

      # For dataset Y conditional on group_vars
      if (!is.null(weight_Y) && weight_Y %in% names(Y)) {
        cond_Y <- Y[, .(temp = list(compute_stats(get(var), get(weight_Y), stats))), by = group_vars]
      } else {
        cond_Y <- Y[, .(temp = list(compute_stats(get(var), NULL, stats))), by = group_vars]
      }
      cond_Y[, dataset := "Y"]
      cond_Y[, variable := var]

      # Updated unlist_stats: bind the list in each row
      unlist_stats <- function(dt) {
        stats_matrix <- do.call(rbind, dt$temp)
        dt <- cbind(dt[, !("temp"), with = FALSE], as.data.table(stats_matrix))
        return(dt)
      }

      cond_X <- unlist_stats(cond_X)
      cond_Y <- unlist_stats(cond_Y)

      cond_summary_list[[var]] <- rbind(cond_X, cond_Y, fill = TRUE)
    }
    conditional_summary <- rbindlist(cond_summary_list, fill = TRUE)
  } else {
    conditional_summary <- NULL
  }

  ## Categorical variables: frequencies and relative frequencies
  categorical_list <- list()
  for (var in cat_vars) {
    # For dataset X
    if (!is.null(weight_X) && weight_X %in% names(X)) {
      freq_X <- X[, .(freq = sum(get(weight_X), na.rm = TRUE)), by = var]
    } else {
      freq_X <- X[, .N, by = var]
      setnames(freq_X, "N", "freq")
    }
    freq_X[, dataset := "X"]
    freq_X[, variable := var]
    freq_X[, rel_freq := freq / sum(freq)]

    # For dataset Y
    if (!is.null(weight_Y) && weight_Y %in% names(Y)) {
      freq_Y <- Y[, .(freq = sum(get(weight_Y), na.rm = TRUE)), by = var]
    } else {
      freq_Y <- Y[, .N, by = var]
      setnames(freq_Y, "N", "freq")
    }
    freq_Y[, dataset := "Y"]
    freq_Y[, variable := var]
    freq_Y[, rel_freq := freq / sum(freq)]

    categorical_list[[var]] <- rbind(freq_X, freq_Y)
  }
  categorical_summary <- rbindlist(categorical_list, fill = TRUE)

  return(list(continuous_summary = continuous_summary,
              conditional_summary = conditional_summary,
              categorical_summary = categorical_summary))
}

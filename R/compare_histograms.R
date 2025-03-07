#' Compare Histograms of a Numeric Variable in Two Datasets
#'
#' This function compares the distribution of a numeric variable between
#' two datasets (e.g., an original dataset and an anonymized/synthetic version)
#' using histograms. It allows conditional faceting on categorical variables
#' and supports different visualization styles, including overlapping, side-by-side,
#' stacked, and separate histograms.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the anonymized/synthetic dataset.
#' @param num_var A character string specifying the numeric variable to compare.
#' @param cat_vars A character vector of categorical variables used for faceting.
#' @param weight_X Optional. A character string specifying the sampling weight variable in X.
#' @param weight_Y Optional. A character string specifying the sampling weight variable in Y.
#' @param facet_type Character string. The faceting type: `"wrap"` (default) for `facet_wrap()` or `"grid"` for `facet_grid()`.
#' @param bins Integer. Number of bins for the histogram (default: `30`).
#' @param transparency Numeric. Transparency level for overlapping histograms (default: `0.4`).
#' @param plot_type Character string. Type of plot:
#'   - `"overlap"` (default): Overlapping histograms with transparency.
#'   - `"side"`: Side-by-side histograms.
#'   - `"stack"`: Stacked histograms.
#'   - `"separate"`: Separate histograms for each dataset-category combination.
#'
#' @importFrom data.table rbindlist
#' @return A ggplot2 histogram visualization comparing the distributions.
#' @export
#'
#' @examples
#' set.seed(123)
#' X <- data.frame(
#'   age = sample(20:80, 500, replace = TRUE),
#'   gender = sample(c("Male", "Female"), 500, replace = TRUE),
#'   income = rnorm(500, mean = 50000, sd = 10000),
#'   weight = runif(500, 0.5, 1.5)
#' )
#'
#' Y <- data.frame(
#'   age = sample(20:80, 1000, replace = TRUE),
#'   gender = sample(c("Male", "Female"), 1000, replace = TRUE),
#'   income = rnorm(1000, mean = 48000, sd = 12000),
#'   weight = runif(1000, 0.5, 1.5)
#' )
#'
#' # Overlapping histograms with transparency
#' compare_histograms(X, Y, num_var = "income", cat_vars = c("gender"),
#'                    weight_X = "weight", weight_Y = "weight",
#'                    facet_type = "wrap", plot_type = "overlap")
#' # Side-by-side histograms
#' compare_histograms(X, Y, num_var = "income", cat_vars = c("gender"),
#'                    weight_X = "weight", weight_Y = "weight",
#'                    facet_type = "wrap", plot_type = "side")
#'
#' # Stacked histograms
#' compare_histograms(X, Y, num_var = "income", cat_vars = c("gender"),
#'                    weight_X = "weight", weight_Y = "weight",
#'                    facet_type = "wrap", plot_type = "stack")
#'
#' # Separate histograms
#' compare_histograms(X, Y, num_var = "income", cat_vars = c("gender"),
#'                    weight_X = "weight", weight_Y = "weight",
#'                    facet_type = "wrap", plot_type = "separate")

compare_histograms <- function(X, Y, num_var, cat_vars,
                               weight_X = NULL, weight_Y = NULL,
                               facet_type = "wrap", bins = 30,
                               transparency = 0.4, plot_type = "overlap") {

  # Convert to data.table for efficient processing
  X <- as.data.table(X)
  Y <- as.data.table(Y)

  # Check if the numeric variable exists in both datasets
  if (!(num_var %in% names(X) && num_var %in% names(Y))) {
    stop("The specified numeric variable must be in both X and Y.")
  }

  # Check if categorical variables exist in both datasets
  if (!all(cat_vars %in% names(X)) || !all(cat_vars %in% names(Y))) {
    stop("All categorical variables must be present in both X and Y.")
  }

  # Assign dataset labels
  X[, .dataset := "Original (X)"]
  Y[, .dataset := "Anonymized/Synthetic (Y)"]

  # Handle weights (if provided, apply; else use unweighted counts)
  if (!is.null(weight_X) && weight_X %in% names(X)) {
    X[, weight := get(weight_X)]
  } else {
    X[, weight := 1]
  }

  if (!is.null(weight_Y) && weight_Y %in% names(Y)) {
    Y[, weight := get(weight_Y)]
  } else {
    Y[, weight := 1]
  }

  # Combine data for plotting
  combined_data <- rbindlist(list(X, Y), use.names = TRUE, fill = TRUE)

  # Define ggplot object
  p <- ggplot(combined_data, aes(x = get(num_var), weight = weight, fill = .dataset)) +
    scale_fill_manual(values = c("blue", "red"), labels = c("Original (X)", "Anonymized/Synthetic (Y)")) +
    theme_minimal() +
    theme(legend.position = "top") +
    labs(x = num_var, y = "Density", fill = "Dataset",
         title = paste("Comparison of Histograms for", num_var))

  # Choose plotting method
  if (plot_type == "overlap") {
    p <- p + geom_histogram(aes(y = after_stat(density)), bins = bins, alpha = transparency, position = "identity")
  } else if (plot_type == "side") {
    p <- p + geom_histogram(aes(y = after_stat(density)), bins = bins, alpha = 0.8, position = "dodge")
  } else if (plot_type == "stack") {
    p <- p + geom_histogram(aes(y = after_stat(density)), bins = bins, alpha = 0.8, position = "stack")
  } else if (plot_type == "separate") {
    # Assign a new facet variable combining dataset type and categorical variables
    combined_data[, facet_group := paste0(.dataset, " - ", get(cat_vars[1]))]

    p <- ggplot(combined_data, aes(x = get(num_var), weight = weight, fill = .dataset)) +
      geom_histogram(aes(y = after_stat(density)), bins = bins, alpha = 0.7, position = "identity") +
      facet_wrap(~facet_group, scales = "free_x") +
      scale_fill_manual(values = c("blue", "red")) +
      theme_minimal() +
      theme(legend.position = "top") +
      labs(x = num_var, y = "Density", fill = "Dataset",
           title = paste("Separate Histograms for", num_var, "by", cat_vars[1]))

    print(p)
    return(invisible(p))
  } else {
    stop("Invalid plot_type. Choose 'overlap', 'side', 'stack', or 'separate'.")
  }

  # Add faceting based on categorical variables (only if not using "separate" mode)
  if (plot_type != "separate") {
    if (length(cat_vars) == 1) {
      if (facet_type == "wrap") {
        p <- p + facet_wrap(as.formula(paste("~", cat_vars)), scales = "free_x")
      } else {
        p <- p + facet_grid(as.formula(paste(". ~", cat_vars)), scales = "free_x")
      }
    } else if (length(cat_vars) == 2) {
      p <- p + facet_grid(as.formula(paste(cat_vars[1], "~", cat_vars[2])), scales = "free_x")
    }
  }

  # Print plot
  print(p)
}

#' Compare continuous distributions conditionally
#'
#' This function calculates and compares the empirical cumulative distribution
#' functions (ECDF) of a sample (\code{X}) and a population (\code{Y}) for
#' specified variables. It can handle weighted samples and provides options for
#' conditional CDFs, and approximation.
#'
#' @param X A data frame or data table.
#' @param Y A data frame or data table to compare with \code{X}.
#' @param variables A character vector specifying the variables for which the
#' comparisons are made. Must be present in both \code{X} and \code{Y}. The
#' class of variables determines the default kind of calculation. See details.
#' @param kind To be set when changing the default kind of calculation. See
#' details.
#' @param conditional A character vector specifying conditioning variables, or
#' \code{NULL} if no conditioning is desired. Must be of length 1 if specified.
#' @param weights A character vector of length 1 specifying the weights variable
#' in \code{X} and/or \code{Y}, or \code{NULL} if no weights are used.
#' @param approx A logical vector of length 2 indicating whether to approximate
#' the CDFs in \code{X} and \code{Y}, respectively. \code{TRUE} is recommended
#' for large data.
#' @param n_approx A numeric value specifying the number of points to use for
#' the approximation.
#' @param bounds A logical value indicating whether to use bounds in the CDF
#' calculation.
#' @param ... additional arguments passed to methods
#' @family comparison
#' @export
#' @return A list with class \code{"compare"} containing the following
#' components:
#' \describe{
#'   \item{\code{formula}}{A formula representing the CDF comparison.}
#'   \item{\code{ecdf}}{A data frame with the ECDF values.}
#'   \item{\code{kind}}{A character string indicating the type of comparison
#'   ("numeric").}
#' }
#'
#' @details
#' The following default calculations are taken:
#'
#' If \code{variables} links to one or two numeric variables: the (weighted)
#' ECDF is calculated for both data sets. With \code{kind} one can change to
#' \code{"density"}, \code{"density_bayes"}, \code{"density_ratio"} and
#' \code{"density_ratio_bayes"}.
#'
#' If \code{variables} links to one, two or three categorical variables:
#' barcharts and mosaic plots are made.
#'
#' If \code{variables} links to one numeric and one categorical variable:
#' boxplots statistics are calculated.
#'
#' If weights are provided, the weighted versions are
#' computed. The function supports conditional estimates, which are computed by
#' conditioning on the specified variable.
#'
#' If the \code{approx} argument is set to \code{TRUE}, the ECDFs are
#' approximated using the specified number of points (\code{n_approx}).
#'
#' @importFrom data.table data.table as.data.table
#'
#' @examples
#' # Example data 1: Complex sample and population data
#' S <- data.frame(
#'   age = sample(20:80, 100, replace = TRUE),
#'   income = rnorm(100, mean = 50000, sd = 10000),
#'   weights = runif(100, 0.5, 2)
#' )
#' U <- data.frame(
#'   age = sample(20:80, 10000, replace = TRUE),
#'   income = rnorm(10000, mean = 50000, sd = 10000)
#' )
#'
#' # Compare ECDFs for age and income
#' result <- compare_distributions_cont(S, U, variables = c("age", "income"), weights = "weights")
#'
#' # View the result
#' print(result$formula)
#' head(result$ecdf)
#' plot(result)
#' plot(result, which = "density")
#' plot(result, which = "density_bayes")
#'
#' # Example data 2: Complex sample and complex sample data
#' X <- data.frame(
#'   age = sample(20:80, 100, replace = TRUE),
#'   income = rnorm(100, mean = 50000, sd = 10000),
#'   weights = runif(100, 0.5, 2)
#' )
#' Y <- data.frame(
#'   age = sample(20:80, 100, replace = TRUE),
#'   income = rnorm(10000, mean = 50000, sd = 10000),
#'   weights = runif(100, 0.5, 2)
#' )
#'
#' # Compare ECDFs for age and income
#' result <- compare_distributions_cont(X, Y, variables = c("age", "income"), weights = "weights")
#'
#' # View the result
#' print(result$formula)
#' head(result$ecdf)
#' plot(result)
#' plot(result, which = "density")
#' plot(result, which = "density_bayes")
#' plot(result, which = "density_ratio")
#'
#' # Example data 3: Sample and sample data
#' X <- data.frame(
#'   age = sample(20:80, 100, replace = TRUE),
#'   income = rnorm(100, mean = 50000, sd = 10000)
#' )
#' Y <- data.frame(
#'   age = sample(20:80, 100, replace = TRUE),
#'   income = rnorm(100, mean = 50000, sd = 10000)
#' )
#'
#' # Compare ECDFs for age and income
#' result <- compare_distributions_cont(X, Y, variables = c("age", "income"), n_approx = 100)
#'
#' # View the result
#' print(result$formula)
#' head(result$ecdf)
#' plot(result)
#' plot(result, which = "density")
#' plot(result, which = "density_bayes")

compare_distributions_cont <- function(X, ...) {
  UseMethod("compare_distributions_cont")
}

#' @rdname compare_distributions_cont
#' @export
compare_distributions_cont.synth_pair <- function(X, ...) {
  compare_distributions_cont.default(X = X$original, Y = X$synthetic, ...)
}

#' @rdname compare_distributions_cont
#' @export
compare_distributions_cont.default <- function(X, Y, variables = NULL, kind = NULL, conditional = NULL,
                          weights = NULL, approx = c(FALSE, TRUE),
                          n_approx = 10000, bounds = TRUE, ...)
{
  # Function to check and convert to data.table if necessary
  convert_to_data_table <- function(X) {
    if (!inherits(X, "data.table")) {
      X <- as.data.table(X)
    }
    return(X)
  }

  # Convert to data.table if necessary
  X <- convert_to_data_table(X)
  Y <- convert_to_data_table(Y)

  n <- nrow(X)
  N <- nrow(Y)

  validate_weights <- function(weights) {
    # Check if weights is NULL
    if (is.null(weights)) {
      return(TRUE)
    }

    # Check if weights is a character vector of length 1
    if (is.character(weights) && length(weights) == 1) {
      return(TRUE)
    }

    # If neither condition is met, throw an error
    stop("weights must be either NULL or a character vector of length 1.")
  }

  validate_weights(weights)

  # for later use:
  n_support <- n_approx

  if(is.null(weights)){
    X$weights <- rep(1, n)
    weights.samp <- weights <- "weights"
  } else{
    weights.samp <- weights
  }

  if (!(weights %in% colnames(X))) {
    stop("The variables names specified in argument 'weight' must be available
         in the sample.\n")
  }
  if (weights %in% colnames(Y)) {
    weights.pop <- "weights"
  } else{
    weights.pop <- NULL
    Y$weights <- rep(1, nrow(Y))
  }

  if (!is.character(variables) || length(variables) == 0) {
    stop("'variables' must be a character vector of positive length!\n")
  }
  if (!(all(variables %in% colnames(Y)) & (all(variables %in% colnames(X))))) {
    stop("The variables names specified in argument 'variables' must be available
         both data sets.\n")
  }
  if (!is.null(conditional) && !is.character(conditional)) {
    stop("'conditional' must be a character vector or NULL!\n")
    if (length(conditional) != 1) {
      stop("argument 'conditional' must have length 1!\n")
    }
  }
  if (!(all(conditional %in% colnames(Y)) & (all(conditional %in% colnames(X))))) {
    stop("The variables names specified in argument 'conditional' must be
         available both data sets.")
  }
  if (!is.logical(approx) || length(approx) == 0)
    approx <- formals()$approx
  else approx <- sapply(rep(approx, length.out = 2), isTRUE)
  if (any(approx) && (!is.numeric(n_approx) || length(n_approx) == 0)) {
    stop("'n_approx' is not numeric or does not have positive length.")
  }
  else n_approx <- ifelse(approx, n_approx[1], NA)
  bounds <- isTRUE(bounds)

  if(nrow(Y) > nrow(X)){
    lab <- c("Sample", "Population")
  } else{
    lab <- c("Sample 1", "Sample 2")
  }

  if (!requireNamespace("simPop", quietly = TRUE)) {
    stop("Package 'simPop' is required for compare_distributions_cont(). Please install it.",
         call. = FALSE)
  }
  tmp <- simPop::getCdf(variables, weights.samp, conditional, X, approx = approx[1],
                n = n_approx[1], name = lab[1])
  values <- tmp$values
  values$.x <- as.numeric(as.character(values$.x))
  app <- t(tmp$approx)
  tmp <- simPop::getCdf(variables, weights.pop, conditional, Y, approx = approx[2],
                n = n_approx[2], name = lab[2])
  values$.x <- as.numeric(as.character(values$.x))
  values <- rbind(values, tmp$values)

  app <- rbind(app, tmp$approx)
  form <- ".y~.x"
  if (length(variables) > 1)
    conditional <- c(".var", conditional)
  if (!is.null(conditional)) {
    conditional <- paste(conditional, collapse = " + ")
    form <- paste(form, conditional, sep = " | ")
  }

  # converting to factors
  values$.name <- factor(values$.name)
  values$.var <- factor(values$.var)

  ## Density

  # Compute density for each combination of sample and variable
  # Compute density for each combination of .name and .var
  # Convert to data.table

  dt <- as.data.table(values)

  ## Compute Bayes densities
  den_1d_num_workhorse <- function(x, y,
                                   weights_x = NULL, weights_y = NULL,
                                   bayesspace = TRUE,
                                   stepsize = 1000) {

    # Normalize weights (if provided), otherwise assume uniform weights
    normalize_weights <- function(w) {
      if (is.null(w)) return(rep(1 / length(x), length(x)))  # Uniform weights
      w <- w / sum(w, na.rm = TRUE)  # Normalize weights to sum to 1
      return(w)
    }

    weights_x <- normalize_weights(weights_x)
    weights_y <- normalize_weights(weights_y)

    # Kernel Density Estimation with weights
    kde_x <- density(x, weights = weights_x, na.rm = TRUE)
    kde_y <- density(y, weights = weights_y, na.rm = TRUE)

    # Define a sequence of points where we want to compute the density ratio
    points <- seq(min(kde_x$x, kde_y$x), max(kde_x$x, kde_y$x), length.out = stepsize)

    # Interpolate densities at the sequence of points
    density_X <- density_X_orig <- approx(kde_x$x, kde_x$y, xout = points, rule = 2)$y
    density_Y <- density_Y_orig <- approx(kde_y$x, kde_y$y, xout = points, rule = 2)$y

    # Handle zero density values to avoid log(0) and division by zero issues
    epsilon <- 1e-10
    density_X[density_X <= 0] <- epsilon
    density_Y[density_Y <= 0] <- epsilon

    if (!bayesspace) {
      distance <- sqrt(sum((density_X - density_Y)^2))
      # Compute the density ratio
      density_ratio <- density_X / density_Y
      kl <- KLDiv(density_X, density_Y)
      jsd <- JSDiv(density_X, density_Y)
    } else {
      # Compute Bayes-space transformation
      gm <- function(x) exp(mean(log(x[x > 0]), na.rm = TRUE))
      density_X <- log(density_X / gm(density_X))
      density_Y <- log(density_Y / gm(density_Y))
      distance <- sqrt(sum((density_X - density_Y)^2, na.rm = TRUE))
      # Compute the density ratio in log-space
      density_ratio <- density_X - density_Y
      kl <- KLDiv_bayes(density_X, density_Y)
      jsd <- JSDiv_bayes(density_X, density_Y)
    }

    mean_ratio <- mean(density_ratio, na.rm = TRUE)
    sd_ratio <- sd(density_ratio, na.rm = TRUE)

    return(list("distance" = distance,
                "density_ratio" = density_ratio,
                "kl" = kl,
                "jsd" = jsd,
                "mean_ratio" = mean_ratio,
                "sd_ratio" = sd_ratio,
                "denX" = density_X,
                "denY" = density_Y,
                "points" = points,
                "bayesspace" = bayesspace))
  }

  # Compute density for each combination of .name and .var using .y
  ll <- 0
  res <- list()
  for(j in variables){
    ll <- ll + 1
    res[[ll]] <- den_1d_num_workhorse(X[[j]], Y[[j]], weights_x = X$weights, weights_y = Y$weights,
                         bayesspace = FALSE)
  }

  # Convert the list into a structured data.table
  extract_densities <- function(res, var_names) {
    # Create a data.table by iterating over the list elements
    dt_list <- lapply(seq_along(res), function(i) {
      data.table(
        .x = res[[i]]$points,          # Extracting points
        .y = c(res[[i]]$denX, res[[i]]$denY),  # Stacking denX and denY
        .name = rep(c("Sample 1", "Sample 2"), each = length(res[[i]]$points)),  # Label samples
        .var = factor(i, levels = seq_along(res), labels = var_names)  # Assign variable names
      )
    })

    # Combine into a single data.table
    rbindlist(dt_list)
  }
  # Convert the list into a structured data.table
  extract_density_ratios <- function(res, var_names) {
    # Create a data.table by iterating over the list elements
    dt_list <- lapply(seq_along(res), function(i) {
      data.table(
        .x = res[[i]]$points,          # Extracting points
        .y = c(res[[i]]$density_ratio, res[[i]]$density_ratio),
        .name = rep(c("Sample 1", "Sample 2"), each = length(res[[i]]$points)),  # Label samples
        .var = factor(i, levels = seq_along(res), labels = var_names)  # Assign variable names
      )
    })

    # Combine into a single data.table
    rbindlist(dt_list)
  }

  # Apply function to res
  densities <- extract_densities(res, variables)
  density_ratios <- extract_density_ratios(res, variables)


  ## for the Bayes space:
  # Compute density for each combination of .name and .var using .y
  ll <- 0
  res_bayes <- list()
  for(j in variables){
    ll <- ll + 1
    res_bayes[[ll]] <- den_1d_num_workhorse(X[[j]], Y[[j]], weights_x = X$weights, weights_y = Y$weights,
                                      bayesspace = TRUE)
  }


  # Apply function to res
  densities_bayes <- extract_densities(res_bayes, variables)
  density_ratios_bayes <- extract_density_ratios(res_bayes, variables)

  results <- list("formula" = form,
               "ecdf" = values,
               "density" = densities,
               "density_bayes" = densities_bayes,
               "density_ratios" = density_ratios,
               "density_ratios_bayes" = density_ratios_bayes,
               "kind" = "numeric",
               "datasetNames" = lab)
  class(results) <- "compare_distributions_cont"
  return(results)
}
NULL


#' Print method for compare_distributions_cont objects
#'
#' @param x an object of class "compare_distributions_cont"
#' @param ... additional arguments (ignored)
#' @export
print.compare_distributions_cont <- function(x, ...) {
  cat("Continuous Distribution Comparison\n")
  cat("==================================\n\n")
  cat("  Formula:", x$formula, "\n")
  cat("  Kind:   ", x$kind, "\n")
  cat("  Datasets:", paste(x$datasetNames, collapse = " vs. "), "\n")
  n_vars <- length(unique(x$ecdf$.var))
  cat("  Variables compared:", n_vars, "\n")
  cat("  ECDF observations:", nrow(x$ecdf), "\n\n")
  cat("  Available plot types: 'ecdf', 'density', 'density_bayes',\n")
  cat("    'density_ratio', 'density_ratio_bayes'\n")
  invisible(x)
}


#' Plot Method for Objects of Class "compare_distributions_cont"
#'
#' This function provides visualizations for objects of class \code{"compare_distributions_cont"}
#' produced by the \code{\link{compare_distributions_cont}} function. It generates plots of the
#' empirical cumulative distribution functions (ECDF) and can optionally
#' produce interactive plots.
#'
#' @param x An object of class \code{"compare_distributions_cont"}.
#' @param ... Additional arguments passed to the plotting function.
#' @param which A character string specifying the type of plot to generate.
#' Options include \code{"ecdf"}, \code{"density"}, \code{"density_bayes"},
#' \code{"density_ratio"}, and \code{"density_ratio_bayes"}.
#' Default is \code{"ecdf"}.
#' @param interactive A logical value indicating whether to produce an
#' interactive plot using \code{plotly}. Default is \code{FALSE}.
#'
#' @return Invisibly returns the ggplot2 object or the plotly object if
#' \code{interactive = TRUE}.
#'
#' @details
#' The \code{plot.compare_distributions_cont} function generates various types of plots for
#' objects of class \code{"compare_distributions_cont"}. The primary plot type is the ECDF plot,
#' which can be generated by setting \code{which} to \code{"ecdf"}.
#' Additional plot types (\code{"density"}, \code{"density_bayes"},
#' \code{"density_ratio"}, and \code{"density_ratio_bayes"}) are currently not
#' implemented and will raise an error if requested.
#'
#' The function supports both static and interactive plots.
#' If \code{interactive = TRUE}, the function will use \code{plotly} to
#' create an interactive version of the plot.
#'
#' @examples
#' S <- data.frame(
#'   age = sample(20:80, 100, replace = TRUE),
#'   income = rnorm(100, mean = 50000, sd = 10000),
#'   weights = runif(100, 0.5, 2)
#' )
#' U <- data.frame(
#'   age = sample(20:80, 10000, replace = TRUE),
#'   income = rnorm(10000, mean = 50000, sd = 10000)
#' )
#'
#' # Compare ECDFs for age and income
#' result <- compare_distributions_cont(S, U, variables = c("age", "income"), weights = "weights")
#'
#' # Plot the ECDF
#' plot(result)
#'
#' \dontrun{
#' # Interactive plot
#' plot(result, interactive = TRUE)
#'
#' ## approx. 20 seconds computation time
#' data("eusilc13puf", package = "simPop")
#' # Function to replace NAs in factor columns with a new level
#' replace_na_in_factor <- function(factor_col, new_level = "not unique") {
#'   # Convert the factor to character
#'   char_col <- as.character(factor_col)
#'   # Replace NA with the new level
#'   char_col[is.na(char_col)] <- new_level
#'   # Convert back to factor and include the new level
#'   factor_col <- factor(char_col, levels = unique(c(levels(factor_col), new_level)))
#'   return(factor_col)
#' }
#'
#' # Apply the function to the relevant columns
#' eusilc13puf$pb220a <- replace_na_in_factor(eusilc13puf$pb220a)
#' eusilc13puf$pl031 <- replace_na_in_factor(eusilc13puf$pl031, new_level = "child")
#' eusilc13puf[is.na(eusilc13puf)] <- 0
#' eusilc13puf$age <- as.numeric(as.character(eusilc13puf$age))
#'
#' inp <- simPop::specifyInput(data=eusilc13puf, hhid="db030",
#'                             hhsize="hsize", strata="db040", weight="rb050")
#' simPop <- simPop::simStructure(data = inp, method = "direct",
#'                               basicHHvars=c("age", "rb090", "hsize", "db040", "pb220a"))
#' simPop <- simPop::simCategorical(simPop, additional=c("pl031"),
#'                                  method = "multinom", nr_cpus = 1)
#' # multinomial model with random draws
#' simPop <- simPop::simContinuous(simPop, additional="pgrossIncome",
#'                                 regModel = ~rb090+hsize+pl031+pb220a+age,
#'                                 upper=200000, equidist=FALSE, nr_cpus=1)
#'
#' eusilc13puf_synth <- data.frame(simPop::pop(simPop))
#' eusilc13puf_synth$age <- as.numeric(as.character(eusilc13puf_synth$age))
#'
#'
#' c1 <- compare_distributions_cont(eusilc13puf, eusilc13puf_synth,
#'               variables = c("age", "pgrossIncome"),
#'               weights = "rb050",
#'               n_approx = 10000)
#' p1 <- plot(c1)
#' p1
#' p1$ggplot_object + theme_dark()
#' c2 <- compare_distributions_cont(eusilc13puf, eusilc13puf_synth,
#'               variables = c("age", "pgrossIncome"),
#'               weights = "rb050",
#'               conditional = "pb220a",
#'               n_approx = 10000)
#' plot(c2)
#' c3 <- compare_distributions_cont(eusilc13puf, eusilc13puf_synth,
#'               variables = c("pgrossIncome"),
#'               weights = "rb050",
#'               conditional = "rb090",
#'               n_approx = 10000)
#' plot(c3)
#' c4 <- compare_distributions_cont(eusilc13puf, eusilc13puf_synth,
#'                     variables = c("age"),
#'                     weights = "rb050",
#'                     conditional = "rb090",
#'                    n_approx = 10000)
#' plot(c4)
#' c5 <- compareSU_cdf(eusilc13puf, eusilc13puf_synth,
#'                     variables = c("age"),
#'                     weights = "rb050",
#'                     conditional = "db040",
#'                     n_approx = 10000)
#' plot(c5)
#' }
#'
#' @export
plot.compare_distributions_cont <- function(x, ..., which = "ecdf", interactive = FALSE){

  if (!which %in% c("ecdf", "density", "density_bayes", "density_ratios","density_ratios_bayes","density_ratio","density_ratio_bayes")) {
    stop("Invalid choice for 'which'. Choose from 'ecdf', 'density', or 'density_bayes'.")
  }
  if(which == "density_ratio") which <- "density_ratios"
  if(which == "density_ratio_bayes") which <- "density_ratios_bayes"

  # Function to convert formula to ggplot2 syntax with color and appropriate faceting
  convert_formula_to_ggplot <- function(data, x_var, y_var, y_label) {

    facet_cmd <- ""

    # Check if faceting is needed
    if (".var" %in% names(data) && length(unique(data$.var)) > 1) {
      facet_cmd <- "facet_wrap(~ .var, scales = 'free')"
    }

    # Create ggplot2 command with color
    ggplot_cmd <- paste0(
      "ggplot(data, aes(x = ", x_var, ", y = ", y_var, ", colour = .name)) + ",
      "geom_line() + theme_bw() + theme(legend.title=element_blank()) + ",
      "xlab('x') + ylab('", y_label, "')",
      if (facet_cmd != "") paste0(" + ", facet_cmd) else ""
    )

    # Evaluate and return the ggplot2 command
    plot <- eval(parse(text = ggplot_cmd))
    return(plot)
  }

  # Select data based on user input
  if (which == "ecdf") {
    plot_data <- x$ecdf
    y_var <- ".y"
    y_label <- "F(x) (ECDF)"
  } else if (which == "density") {
    plot_data <- x$density
    y_var <- ".y"
    y_label <- "Density"
  } else if (which == "density_bayes") {
    plot_data <- x$density_bayes #x$density
    y_var <- ".y" #".y_bayes"
    y_label <- "Bayes-Transformed Density"
  } else if (which == "density_ratios_bayes") {
    plot_data <- x$density_ratios #x$density
    y_var <- ".y" #".y_bayes"
    y_label <- "Bayes-Transformed Density ratio"
  } else if (which == "density_ratios") {
    plot_data <- x$density_ratios_bayes #x$density
    y_var <- ".y" #".y_bayes"
    y_label <- "Density ratio"
  }

  # Generate the plot
  density_plot <- convert_formula_to_ggplot(plot_data, ".x", y_var, y_label)

  if(which %in% c("density_ratios","density_ratios_bayes")){
    density_plot <- density_plot + geom_hline(yintercept = 0, linetype = "dashed", color = "red") + theme(legend.position="none")
  }
  # Display the plot interactively or statically
  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Package 'plotly' is required for interactive plots. Falling back to static plot.")
      print(density_plot)
    } else {
      print(plotly::ggplotly(density_plot))
    }
  } else {
    print(density_plot)
  }

  invisible(density_plot)
}




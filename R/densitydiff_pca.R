#' Density ratio and entropy of first PC's
#'
#' Computes the density ratio between the selected principal components
#' and provides entropy measures and more.
#'
#' @rdname densitydiff_pca
#' @author Matthias Templ
#' @param X a numeric vector
#' @param Y a numeric vector
#' @param bayesspace if TRUE, a Bayes space compositional density approach
#' is taken
#' @param stepsize size of sequence of points where the densities are evaluated
#' @param strata_x vector holding the information which observation is related
#' to which strata. If it is not set to NULL,
#' density estimation is applied on each strata.
#' @param strata_y vector holding the information which observation is related
#' to which strata. If not NULL, density estimation is applied on each strata.
#' It must have the same levels than strata_x and there must be observed values
#' for each level.
#' @importFrom stats princomp
#' @importFrom stats runif
#' @family comparison
#' @export
#' @examples
#' # Simple example with multivariate numeric data
#' set.seed(123)
#' X <- data.frame(
#'   var1 = rnorm(100, 50, 10),
#'   var2 = rnorm(100, 1000, 200),
#'   var3 = rnorm(100, 30, 5),
#'   var4 = rnorm(100, 500, 100)
#' )
#' Y <- data.frame(
#'   var1 = rnorm(100, 52, 11),
#'   var2 = rnorm(100, 980, 220),
#'   var3 = rnorm(100, 31, 6),
#'   var4 = rnorm(100, 510, 110)
#' )
#' d1 <- densitydiff_pca(X, Y, bayesspace = FALSE)
#' d1
densitydiff_pca <- function(X,
                            Y,
                            bayesspace = TRUE,
                            stepsize = 1000,
                            strata_x = NULL,
                            strata_y = NULL){

  check_matrices_or_dataframes_and_numeric <- function(X, Y) {
    # Check if X and Y are matrices or data frames
    if (!is.matrix(X) && !is.data.frame(X)) {
      stop("X is neither a matrix nor a data frame.")
    }
    if (!is.matrix(Y) && !is.data.frame(Y)) {
      stop("Y is neither a matrix nor a data frame.")
    }

    # Check if all entries in X are numeric
    if (!all(sapply(X, is.numeric))) {
      stop("Data frame or matrix X contains non-numeric entries.")
    }

    # Check if all entries in Y are numeric
    if (!all(sapply(Y, is.numeric))) {
      stop("Data frame or matrix Y contains non-numeric entries.")
    }
  }

  check_matrices_or_dataframes_and_numeric(X, Y)

  # if(!is.null(weight_x)){
  #   # Check if weight_x is a numeric vector
  #   if (!is.numeric(weight_x) || !is.vector(weight_x)) {
  #     stop("Error: 'weight_x' must be a numeric vector.")
  #   }
  #   # Check if weight_x is of the same length than x
  #   if(!identical(length(x), length(weights_x))) stop("Error: 'lenght of weight_x must be equal to x.")
  #   # Check if weight_y is of the same length than y
  #   if(!identical(length(y), length(weights_y))) stop("Error: 'lenght of weight_y must be equal to y.")
  # }
  #
  # if(!is.null(weight_y)){
  #   # Check if weight_y is a numeric vector
  #   if (!is.numeric(weight_y) || !is.vector(weight_y)) {
  #     stop("Error: 'weight_y' must be a numeric vector.")
  #   }
  # }

  # Check if bayesspace is a logical value
  if (!is.logical(bayesspace) || length(bayesspace) != 1) {
    stop("Error: 'bayesspace' must be a single logical value, TRUE or FALSE.")
  }

  # Check if stepsize is a positive integer
  if (!is.numeric(stepsize) || length(stepsize) != 1 || stepsize <= 0 || stepsize %% 1 != 0) {
    stop("Error: 'stepsize' must be a positive integer.")
  }

  # Check if strata_x is a factor or NULL
  if (!is.null(strata_x)) {
    if (!is.factor(strata_x) || length(strata_x) != nrow(X)) {
      stop("Error: 'strata_x' must be a factor with the same length as 'X'.")
    }
  }

  # Check if strata_y is a factor or NULL
  if (!is.null(strata_y)) {
    if (!is.factor(strata_y) || length(strata_y) != nrow(Y)) {
      stop("Error: 'strata_y' must be a factor with the same length as 'Y'.")
    }

    if (!is.null(strata_x) && !all(levels(strata_y) %in% levels(strata_x))) {
      stop("Error: 'strata_y' must have the same levels as 'strata_x'.")
    }
  }

  # Check if strata_x and strata_y are both NULL or both non-NULL
  if (is.null(strata_x) != is.null(strata_y)) {
    stop("Error: Both 'strata_x' and 'strata_y' must be NULL or both must be
         non-NULL.")
  }

  # Check for empty levels in strata_y
  if (!is.null(strata_x)) {
    if (any(table(strata_x) == 0)) {
      stop("Error: 'strata_x' contains levels with no observations.")
    }
  }

  # Check for empty levels in strata_y
  if (!is.null(strata_y)) {
    if (any(table(strata_y) == 0)) {
      stop("Error: 'strata_y' contains levels with no observations.")
    }
  }

  den_pca_workhorse <- function(X,
                                Y,
                                choices = 1L:2L,
                                bayesspace = TRUE,
                                stepsize = 1000){
    # Estimate the selected principal components of X and Y
    scoresX <- princomp(X)$scores[, choices]
    scoresY <- princomp(Y)$scores[, choices]
    k <- length(choices)
    kl <- jsd <- mean_ratio <- sd_ratio <- numeric(k)
    for(i in 1:k){
      # Kernel Density Estimation
      kde_x <- density(scoresX[, i])
      kde_y <- density(scoresY[, i])

      # Define a sequence of points where we want to compute the density ratio
      points <- seq(min(kde_x$x, kde_y$x), max(kde_x$x, kde_y$x),
                    length.out = stepsize)

      # Interpolate densities at the sequence of points
      density_X <- density_X_orig <- approx(kde_x$x, kde_x$y, xout = points)$y
      density_Y <- density_Y_orig <- approx(kde_y$x, kde_y$y, xout = points)$y



      if(!bayesspace){
        distance <- sqrt(sum((density_X - density_Y)^2))
        # Compute the density ratio
        density_ratio <- density_X / density_Y
        #kl <- sum(density_X * log(density_ratio))
        kl[i] <- KLDiv(density_X, density_Y)
        jsd[i] <- JSDiv(density_X, density_Y)
      } else{
        # Geometric mean (inline to avoid robCompositions dependency)
        gm <- function(x) exp(mean(log(x[x > 0]), na.rm = TRUE))
        density_X <- log(density_X / gm(density_X))
        density_Y <- log(density_Y / gm(density_Y))
        distance <- sqrt(sum((density_X - density_Y)^2, na.rm=TRUE))
        # Compute the density ratio
        density_ratio <- density_X - density_Y
        # kl <- length(density_X_orig) / 2 * log(mean(density_X_orig / density_Y_orig, na.rm = TRUE) * mean(density_Y_orig / density_X_orig, na.rm = TRUE))
        kl[i] <- KLDiv_bayes(density_X, density_Y)
        jsd[i] <- JSDiv_bayes(density_X, density_Y)
      }

      mean_ratio[i] <- mean(density_ratio, na.rm = TRUE)
      sd_ratio[i] <- sd(density_ratio, na.rm = TRUE)
    }

    return(list(#"distance" = distance,
                #"density_ratio" = density_ratio,
                "kl" = kl,
                "jsd" = jsd,
                "mean_ratio" = mean_ratio,
                "sd_ratio" = sd_ratio#,
                #"denX" = density_X,
                #"denY" = density_Y,
                #"points" = points,
                #"bayesspace" = bayesspace,
                #"strata_x" = strata_x,
                #"strata_y" = strata_y
                )
    )
  }

  if(is.null(strata_x)){
    result <- den_pca_workhorse(X, Y, choices = 1:2, bayesspace, stepsize)
  } else {
    result <- list()
    j <- 0
    for(i in levels(strata_x)){
      j <- j + 1
      result[[j]] <- den_pca_workhorse(X[strata_x == i], Y[strata_y == i], bayesspace, stepsize)
    }
  }
  class(result) <- "denpca"
  return(result)
}

# Register global variables to avoid NOTE
utils::globalVariables(c("x", "y"))

NULL


#' Print method for denpca objects
#'
#' @param x an object of class "denpca"
#' @param ... additional arguments (ignored)
#' @export
print.denpca <- function(x, ...) {
  cat("Density Ratio Comparison (PCA-based)\n")
  cat("====================================\n\n")
  if (is.list(x) && !is.null(x$kl)) {
    n_pc <- length(x$kl)
    cat("  Principal components:", n_pc, "\n")
    for (i in seq_len(n_pc)) {
      cat(sprintf("  PC%d: KL = %.4f, JSD = %.4f, mean ratio = %.4f\n",
                  i, x$kl[i], x$jsd[i], x$mean_ratio[i]))
    }
  } else {
    cat("  Stratified result with", length(x), "strata\n")
    for (i in seq_along(x)) {
      cat(sprintf("  Stratum %d: KL = [%s], JSD = [%s]\n",
                  i,
                  paste(sprintf("%.4f", x[[i]]$kl), collapse = ", "),
                  paste(sprintf("%.4f", x[[i]]$jsd), collapse = ", ")))
    }
  }
  cat("\n")
  invisible(x)
}


#' Plot method for denpca objects
#'
#' This function plots objects of class \code{denpca}.
#'
#' @param x An object of class \code{denpca}.
#' @param y Not used.
#' @param ... Additional arguments passed to the plotting functions.
#' @param which Which plot to show: 1 for density, 2 for density ratio.
#' @return A plot.
#' @rdname plot.denpca
#' @method plot denpca
#' @export
plot.denpca <- function(x, y, ..., which = 1){
  # TBD
  is_list_of_lists <- function(obj) {
    # Check if obj is a list
    if (!is.list(obj)) {
      return(FALSE)
    }

    # Check if every element of obj is a list
    for (element in obj) {
      if (!is.list(element)) {
        return(FALSE)
      }
    }

    return(TRUE)
  }

  show <- rep(FALSE, 2)
  show[which] <- TRUE

  if(is_list_of_lists(x)){
    if(x[[1]]$bayesspace){
      ylab1 <- "densities (Bayes space)"
      ylab2 <- "density ratio (Bayes space)"
    } else {
      lab1 <- "densities"
      lab2 <- "density ratio"
    }
    extractFromList <- function(x, what){
      df <- data.frame(sapply(x, function(x) x[[what]]))
      colnames(df) <- levels(x[[1]]$strata_x)
      df <- reshape2::melt(df)
      df$what <- what
      colnames(df)[1] <- "strata"
      return(df)
    }
    appendMe <- function(dfNames) {
      do.call(rbind, lapply(dfNames, function(x) {
        cbind(get(x), source = x)
      }))
    }

    denX <- extractFromList(x, what = "denX")
    denY <- extractFromList(x, what = "denY")
    points <- extractFromList(x, what = "points")
    denX$x <- points$value
    denY$x <- points$value
    #    den <- appendMe(c("denX", "denY"))
    den <- rbind(denX, denY)
    denRatio <- extractFromList(x, what = "density_ratio")
    denRatio$x <- points$value

    if(show[1L]){
      print(ggplot(den, aes(x = x, y = value, colour = what)) +
              geom_line() +
              facet_wrap(~ strata) +
              ylab(ylab2) +
              theme_minimal() +
              theme(legend.position="none") +
              geom_hline(yintercept = 0, color = "grey", linetype = "dashed"))
    }
    if(show[2L]){
      print(ggplot(denRatio, aes(x = x, y = value, colour = what)) +
              geom_line() +
              facet_wrap(~ strata) +
              ylab(ylab2) +
              theme_minimal() +
              theme(legend.position="none") +
              geom_hline(yintercept = 0, color = "grey", linetype = "dashed"))
    }

  } else {
    if(x$bayesspace){
      ylab1 <- "densities (Bayes space)"
      ylab2 <- "density ratio (Bayes space)"
    } else {
      ylab1 <- "density"
      ylab2 <- "density ratio"
    }
    if(sum(show) > 1){
      par(mfrow = c(1,2))
    }
    if(show[1L]){
      plot(x = x$points, y = x$denX, type = "l", ylab = ylab1, xlab = "")
      lines(x = x$points, y = x$denY, col = "red")
      abline(h = 1, lty = 2, col = "gray")
    }
    if(show[2L]){
      plot(x = x$points,
           y = x$density_ratio,
           type = "l",
           ylab = ylab2, xlab = "")
      abline(h = 1, lty = 2, col = "gray")
    }

  }
}






# plot(d1$points, d1$density_ratio)
# plot(d1e$points, d1e$density_ratio)
#

# d2 <- densitydiff_1d_num(x, y, strata_x = strata_x, strata_y = strata_y)
#plot(d1)
# plot(d1, which = 1:2)
# plot(d1e, which = 1:2)


#plot(d2, which = 1:2)

#
# # Plot the results
# par(mfrow = c(2,2))
# plot(kde_x)
# lines(kde_y, col = "blue")
# plot(points, density_ratio, type = "l", col = "red",
#      xlab = "Points", ylab = "Density Ratio",
#      main = "Density Ratio Between Two Distributions")
# abline(h = 1, lty = 2, col = "gray")
#
# plot(points, density_X_clr, type = "l")
# lines(points, density_Y_clr, col = "blue")
# abline(h = 0, lty = 2, col = "gray")
# plot(points, density_ratio_clr, type = "l", col = "red",
#      xlab = "Points", ylab = "Density Ratio",
#      main = "Density Ratio Between Two Distributions")
# abline(h = 1, lty = 2, col = "gray")


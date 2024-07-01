#' Propensity score utility measure
#'
#' How well can you predict if an observation is from which data set?
#'
#' @param X data frame
#' @param Y data frame with the same structure as X
#' @param form formula. If NULL all variables are used as predictors
#' @author Matthias Templ
#' @return Propensity score measures
#' @references
#' Templ, M. Statistical Disclosure Control for Microdata: Methods and Applications in R.
#' \emph{Springer International Publishing}, 287 pages, 2017. ISBN 978-3-319-50272-4. \doi{10.1007/978-3-319-50272-4}
#' data(eusilc13puf, package = "simPop")
#' eusilc13puf$age <- as.numeric(as.character(eusilc13puf$age))
#' keyvars <- c("age", "rb090", "db040", "pl031", "pb220a")
#' require(sdcMicro)
#' sdc <- createSdcObj(eusilc13puf,
#'          keyVars = keyvars,
#'          numVars = "pgrossIncome",
#'          w = "rb050",
#'          hhId = "db030")
#' sdc <- globalRecode(sdc,
#'          column = "age",
#'          breaks = c(1,9,19,29,39,49,59,69,100),
#'          labels = 1:8)
#' sdc <- localSuppression(sdc)
#' sdc <- microaggregation(sdc)
#' eusilc13puf_anon <- extractManipData(sdc)
#'
#' propscore(eusilc13puf, eusilc13puf_anon, na = "remove")
#'
#'
#' #' \dontrun{
#' ## approx. 20 seconds computation time
#' require(simPop)
#' inp <- specifyInput(data=eusilc13puf, hhid="db030", hhsize="hsize", strata="db040", weight="rb050")
#' simPop <- simStructure(data = inp, method = "direct",
#'   basicHHvars=c("age", "rb090", "hsize", "db040"))
#' simPop <- simCategorical(simPop, additional=c("pl031", "pb220a"), method = "multinom", nr_cpus = 1)
#' # multinomial model with random draws
#' simPop <- simContinuous(simPop, additional="pgrossIncome",
#'               regModel = ~rb090+hsize+pl031+pb220a,
#'               upper=200000, equidist=FALSE, nr_cpus=1)
#' eusilc13puf_synth <- data.frame(pop(simPop))
#'
#' propscore(eusilc13puf, eusilc13puf_synth, na = "remove")
#' propscore(eusilc13puf, eusilc13puf_synth, na = "remove", cluster = "db030")
#'
#' }
propscore <- function(X, Y,
               form = NULL,
               method = "rf",
               adjust_size = TRUE,
               cluster = NULL,
               na = "impute"){
  # Check if X and Y are data frames
  if (!is.data.frame(X)) {
    stop("X must be a data frame.")
  }
  if (!is.data.frame(Y)) {
    stop("Y must be a data frame.")
  }

  if (is.null(form)) {
    # If no formula is provided, check if X and Y have the same structure
    if (!all(names(X) == names(Y))) {
      stop("X and Y must have the same structure.")
    }
    if (!all(sapply(X, class) == sapply(Y, class))) {
      stop("X and Y must have the same column classes.")
    }
    form <- paste("group", paste(colnames(X), collapse = " + "), collapse = " ~ ")
  } else {
    # Check if form is of class 'formula'
    if (!inherits(form, "formula")) {
      stop("form must be a valid formula.")
    }

    # Extract variable names from the formula
    formula_vars <- all.vars(form)[-1]

    # Check if variables in the formula appear in X (and Y)
    missing_vars <- setdiff(formula_vars, names(X))
    if (length(missing_vars) > 0) {
      stop(paste("The following variables from the formula are missing in X (and Y):", paste(missing_vars, collapse = ", ")))
    }

    # Check if the columns in X and Y for the formula variables have the same structure
    for (var in formula_vars) {
      if (!identical(class(X[[var]]), class(Y[[var]]))) {
        stop(paste("The variable", var, "must have the same class in both X and Y."))
      }
    }
  }

  # Check the cluster argument
  if (!is.null(cluster)) {
    if (!is.character(cluster) || length(cluster) != 1) {
      stop("cluster must be a single string representing a variable name in X.")
    }
    if (!(cluster %in% names(X))) {
      stop(paste("The cluster variable", cluster, "is not found in X."))
    }
  }

  # Check the na argument
  valid_na_values <- c("impute", "remove", "stop")
  if (!(na %in% valid_na_values)) {
    stop(paste("na must be one of", paste(valid_na_values, collapse = ", "), "."))
  }

  vars <- attributes(terms.formula(form))$term.labels
  if(any(is.na(X[, vars]))){
    if(na == "stop") stop("data contains missing values")
    if(na == "impute"){
      X[, vars] <- VIM::kNN(X[, vars], imp_var = FALSE)
      Y[, vars] <- VIM::kNN(Y[, vars], imp_var = FALSE)
    }
    if(na == "remove"){
      X <- X[complete.cases(X[, vars]), ]
      Y <- Y[complete.cases(Y[, vars]), ]
    }
  }

  # Function to sample data based on the cluster variable
  sample_by_cluster <- function(X, Y, cluster) {
    X_clusters <- unique(X[[cluster]])
    Y_clusters <- unique(Y[[cluster]])

    if (nrow(X) > nrow(Y)) {
      # Sample households from X to match the number of persons in Y
      sampled_clusters <- sample(X_clusters, length(Y_clusters))
      X_sampled <- X[X[[cluster]] %in% sampled_clusters, ]
      return(list(X = X_sampled, Y = Y))
    } else {
      # Sample households from Y to match the number of persons in X
      sampled_clusters <- sample(Y_clusters, length(X_clusters))
      Y_sampled <- Y[Y[[cluster]] %in% sampled_clusters, ]
      return(list(X = X, Y = Y_sampled))
    }
  }

  # Sample the data based on the cluster variable if cluster is provided
  if (!is.null(cluster)) {
    sampled_data <- sample_by_cluster(X, Y, cluster)
    X <- sampled_data$X
    Y <- sampled_data$Y
  }

  mi_x <- n_x <- nrow(X)
  mi_y <- n_y <- nrow(Y)

  cr <- n_y / (n_x + n_y)

  if(adjust_size & !(identical(n_x, n_y))){
#    ma <- max(c(n_x, n_y))
    if(n_x > n_y){
      index0 <- sample(1:n_x, size = n_y, replace = FALSE)
      X <- X[index0, ]
    } else {
      index0 <- sample(1:n_y, size = n_x, replace = FALSE)
      Y <- Y[index0, ]
    }
    cr <- nrow(Y) / (nrow(X) + nrow(Y))
    mi_x <- mi_y <- min(c(n_x, n_y))
  }

    Z <- rbind(X[, vars], Y[, vars])
    Z$group <- factor(rep(c(0,1), times = c(mi_x, mi_y)))
    if(method == "logreg"){
      res <- glm(form, data = Z, family = binomial())
      p <- data.frame("prediction" = predict(res, type = "response"),
                      "group" = rep(c("real", "synth"), times = c(mi_x, mi_y)))
    } else{
      res <- randomForest::randomForest(form, data = Z)
      # p <- data.frame("prediction" = predict(res, data = Z, predict.all = TRUE)$pred,
      #                 "group" = rep(c("real", "synth"), times = c(mi_x, mi_y)))
      p <- data.frame("prediction" = res$votes[,2],
                      "group" = rep(c("real", "synth"), times = c(mi_x, mi_y)))
    }

    ps <- 1 / (2 * mi_x) * sum((p$prediction - cr)^2)
    t1 <- as.numeric(table(p$prediction < cr))
    ps_ratio <- t1[1] / t1[2]

    # Kernel Density Estimation, difference in distributions
    kde_x <- density(p$prediction[1:mi_x])
    kde_y <- density(p$prediction[(mi_x+1):length(ps)])
    # Define a sequence of points where we want to compute the density ratio
    points <- seq(min(kde_x$x, kde_y$x), max(kde_x$x, kde_y$x), length.out = 1000)
    # Interpolate densities at the sequence of points
    density_X <- density_X_orig <- approx(kde_x$x, kde_x$y, xout = points)$y
    density_Y <- density_Y_orig <- approx(kde_y$x, kde_y$y, xout = points)$y
    distance <- sqrt(sum((density_X - density_Y)^2, na.rm = TRUE))
    # Compute the density ratio
    density_ratio <- density_X / density_Y
    kl <- sum(density_X * log(density_ratio), na.rm = TRUE)
    density_X <- log(density_X / robCompositions::gm(density_X))
    density_Y <- log(density_Y / robCompositions::gm(density_Y))
    distance <- sqrt(sum((density_X - density_Y)^2, na.rm=TRUE))
    # Compute the density ratio
    density_ratio_bayes <- density_X - density_Y
    kl_bayes <- length(density_X_orig) / 2 * log(mean(density_X_orig / density_Y_orig, na.rm = TRUE) * mean(density_Y_orig / density_X_orig, na.rm = TRUE))
    mean_ratio <- mean(density_ratio, na.rm = TRUE)
    sd_ratio <- sd(density_ratio, na.rm = TRUE)
    mean_ratio_bayes <- mean(density_ratio_bayes, na.rm = TRUE)
    sd_ratio_bayes <- sd(density_ratio_bayes, na.rm = TRUE)


    results <- list("predictions" = p,
                "ps_ratio" = ps_ratio,
                "ps_score" = ps,
                "cr" = cr,
                "mean_ps_x" = mean(p$prediction[1:mi_x], na.rm = TRUE),
                "mean_ps_y" = mean(p$prediction[(mi_x + 1):length(p$prediction)], na.rm = TRUE),
                "density_ratio" = density_ratio,
                "density_ratio_bayes" = density_ratio_bayes,
                "kl" = kl,
                "kl_bayes" = kl_bayes,
                "mean_ratio" = mean_ratio,
                "sd_ratio" = sd_ratio,
                "mean_ratio_bayes" = mean_ratio_bayes,
                "sd_ratio_bayes" = sd_ratio_bayes)
    class(results) <- "propscore"
    return(results)

}

print.propscore <- function(x, ...){
  cat("mean propensity scores for x: ", x$mean_ps_x)
  cat("\nmean propensity scores for y: ", x$mean_ps_y)
  cat("\npropensity score statistics: ", x$ps_score, "\n")
}

summary.propscore <- function(object, ...){
  cat("TBD")
}


plot.propscore <- function(x, y, ..., which = 1){

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



# data(eusilc13puf, package = "simPop")
# eusilc13puf$age <- as.numeric(as.character(eusilc13puf$age))
# keyvars <- c("age", "rb090", "db040", "pl031", "pb220a")
# require(sdcMicro)
# sdc <- createSdcObj(eusilc13puf,
#          keyVars = keyvars,
#          numVars = "pgrossIncome",
#          w = "rb050",
#          hhId = "db030")
# sdc <- globalRecode(sdc,
#          column = "age",
#          breaks = c(1,9,19,29,39,49,59,69,100),
#          labels = 1:8)
# sdc <- localSuppression(sdc)
# sdc <- microaggregation(sdc)
# eusilc13puf_anon <- extractManipData(sdc)
#
# propscore(eusilc13puf, eusilc13puf_anon, na = "remove")
#
#
# ## approx. 20 seconds computation time
# require(simPop)
# inp <- specifyInput(data=eusilc13puf, hhid="db030", hhsize="hsize", strata="db040", weight="rb050")
# simPop <- simStructure(data = inp, method = "direct",
#   basicHHvars=c("age", "rb090", "hsize", "db040"))
# simPop <- simCategorical(simPop, additional=c("pl031", "pb220a"), method = "multinom", nr_cpus = 1)
# # multinomial model with random draws
# simPop <- simContinuous(simPop, additional="pgrossIncome",
#               regModel = ~rb090+hsize+pl031+pb220a,
#               upper=200000, equidist=FALSE, nr_cpus=1)
# eusilc13puf_synth <- data.frame(pop(simPop))
#
# propscore(eusilc13puf, eusilc13puf_synth, na = "remove")
# propscore(eusilc13puf, eusilc13puf_synth, na = "remove", method = "rf")
# propscore(eusilc13puf, eusilc13puf_synth, na = "remove", adjust_size = FALSE)
# propscore(eusilc13puf, eusilc13puf_synth, na = "remove", adjust_size = FALSE, method = "rf")
# propscore(eusilc13puf, eusilc13puf_synth, na = "remove", method = "rf")
# propscore(eusilc13puf, eusilc13puf_synth, na = "remove", cluster = "db030")
#



#' Comparison of confidence intervals
#'
#' Given a sample, how much higher is the change compared to the sampling error.
#'
#' @author Matthias Templ
#' @param x a vector
#' @param y a vector
#' @param estimator a function with a vector of lenght 1 as output
#' @param R number of bootstrap replicates
#' @importFrom stats median
#' @importFrom stats quantile
#' @family comparison
#' @export
#' @rdname ci_overlap
#' @return A list with the following elements:
#' \describe{
#'   \item{m_x}{The point estimate of the function provided in function
#'   argument estimator for the vector x.}
#'   \item{m_y}{The point estimate of the function provided in function
#'   argument estimator for the vector x.}
#'   \item{ci_x}{The confidence interval of x.}
#'   \item{ci_y}{The confidence interval of y.}
#'   \item{overlap}{Overlap in confidence intervals.}
#'   \item{overlap_perc}{Overlap in confidence intervals in percentage}
#'   \item{ratioSE}{How much differ is the point estimate of y compared to the
#'   half of the confidence interval.}
#' }
#' @examples
#' x <- rnorm(100, 10)
#' y <- x + runif(100, 0, 1)
#'
#' ci_overlap(x, y, R = 1000)
#'
#' x <- factor(sample(1:2, size = 100, replace = TRUE))
#' y <- factor(sample(1:2, size = 100, replace = TRUE))
#'
#' my_func <- function(x){
#'   table(x)[1]
#' }
#'
#' ci_overlap(x, y, estimator = my_func)
#'
ci_overlap <- function(x, y, estimator = median, R = 1000){
  # calculate the estimator on x
  m_x <- estimator(x)
  if(length(m_x) > 1) stop("estimator must return a vector of lenght 1")
  # calculate the confidence interval
  nr <- length(x)
  ci_x <- quantile(replicate(R, estimator(sample(x, size = nr, replace = TRUE))),
                 c(0.025, 0.975))
  # calculate the estimator on y
  m_y <- estimator(y)
  # calculate the confidence interval
  nr_y <- length(y)
  ci_y <- quantile(replicate(R, estimator(sample(y, size = nr_y, replace = TRUE))),
                 c(0.025, 0.975))
  # is it less than the sampling error?
  lessSE <- ifelse(m_x > ci_x[1] & m_x < ci_x[2], "smaller", "larger")
  ## overlap of ci's
  # Determine the overlapping interval
  lower_bound_overlap <- max(ci_x[1], ci_y[1])
  upper_bound_overlap <- min(ci_x[2], ci_y[2])

  # Check if there is an overlap
  if (lower_bound_overlap < upper_bound_overlap) {
    overlap_length <- upper_bound_overlap - lower_bound_overlap
    total_length <- max(ci_x[2], ci_y[2]) - min(ci_x[1], ci_y[1])
    overlap <- (overlap_length / total_length)
    percentage_overlap <- overlap * 100
  } else {
    overlap <- percentage_overlap <- 0
  }

  # ratio SE to est_y
  ratioSE <- abs(m_y - m_x) / abs((ci_x[2] - ci_x[1])/2)

  return(list("m_x" = m_x,
              "m_y" = m_y,
              "ci_x" = ci_x,
              "ci_y" = ci_y,
              "overlap" = overlap,
              "overlap_perc" = percentage_overlap,
              "ratioSE" = ratioSE))
}



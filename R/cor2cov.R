
#' Compute covariance matrix from correlation matrix
#'
#' R and vector of standard deviations S. covariance matrix
#' required as input for mvrnorm
#'
cor2cov <- function(R, S) {
  sweep(sweep(R, 1, S, "*"), 2, S, "*")
}

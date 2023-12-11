
#' Compute covariance matrix
#'
#' @param R correlation matrix
#' @param S vector of standard deviations covariance matrix
#' @return Required input for \code{mvrnorm}.
#' @keywords internal
#'
cor2cov <- function(R, S) {
  sweep(sweep(R, 1, S, "*"), 2, S, "*")
}


#' Compute covariance matrix
#'
#' @param cormat correlation matrix
#' @param S vector of standard deviations covariance matrix
#' @return Required input for \code{mvrnorm}.
#' @keywords internal
#'
cor2cov <- function(cormat, S) {
  # cormat * S %*% t(S)  # alternative
  sweep(sweep(cormat, 1, S, "*"), 2, S, "*")
}

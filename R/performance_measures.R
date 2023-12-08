
## Functions to evaluate performance measures

#' Objective function to minimize for standard method of moments MAIC
#'
#' @param beta Beta
#' @param X X
#' @keywords internal
#' 
Q <- function(beta, X) {
  sum(exp(X %*% beta))
}

#' Bias estimate
#' 
#' @param theta.hat Theta hat
#' @param theta Theta
#' @keywords internal
#' 
bias <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  sum(theta.hat)/nsim - theta
}

#' Monte Carlo SE of bias estimate
#' @param theta.hat theta hat
#' @return \eqn{sqrt(1/(nsim*(nsim-1))*tmp)}
#' @keywords internal
#' 
bias.mcse <- function(theta.hat) {
  nsim <- length(theta.hat)
  tmp <- sum((theta.hat - mean(theta.hat))^2)
  sqrt(1/(nsim*(nsim-1))*tmp)
}

#' Coverage estimate
#' 
#' @param low Low
#' @param upp Upper
#' @param theta Theta
#' @return \eqn{sum(in_range)/nsim}
#' @keywords internal
#' 
coverage <- function(low, upp, theta) {
  nsim <- length(low)
  theta_inside_range <- theta >= low & theta <= upp
  in_range <- ifelse(theta_inside_range, 1, 0)
  sum(in_range)/nsim
}

#' Monte Carlo SE of coverage estimate
#' 
#' @param coverage Coverage
#' @param nsim Number of simulations
#' @return \eqn{sqrt((coverage*(1 - coverage))/nsim)}
#' @keywords internal
#' 
coverage.mcse <- function(coverage, nsim) {
  sqrt((coverage*(1 - coverage))/nsim)
}

#' Mean squared error estimate
#'
#' @param theta.hat Theta hat
#' @param theta Theta
#' @return \eqn{sum((theta.hat - theta)^2)/nsim}
#' @keywords internal
#' 
mse <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  sum((theta.hat - theta)^2)/nsim
}

#' Monte Carlo SE of MSE estimate
#' 
#' @param theta.hat Theta hat
#' @param theta Theta
#' @return \eqn{sqrt(sum((tmp - mse.est)^2)/(nsim*(nsim-1)))}
#' @keywords internal
#' 
mse.mcse <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  tmp <- (theta.hat - theta)^2
  mse.est <- sum(tmp)/nsim
  sqrt(sum((tmp - mse.est)^2)/(nsim*(nsim-1)))
}

#' Mean absolute error estimate
#' 
#' @param theta.hat Theta hat
#' @param theta Theta
#' @return \eqn{sum(abs(theta.hat - theta))/nsim}
#' @keywords internal
#' 
mae <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  sum(abs(theta.hat - theta))/nsim
}

#' Monte Carlo SE of any continuous performance metric
#' 
#' @param pm pm
#' @return \eqn{sqrt(sum((pm - pm_mean)^2)/(nsim*(nsim-1)))}
#' @keywords internal
#' 
mcse.estimate <- function(pm) {
  nsim <- length(pm)
  pm_mean <- sum(pm)/nsim
  sqrt(sum((pm - pm_mean)^2)/(nsim*(nsim-1)))
}

#' Empirical standard error 
#' 
#' @param theta.hat Theta
#' @return \eqn{sqrt(tmp/(nsim-1))}
#' @keywords internal
#' 
empse <- function(theta.hat) {
  nsim <- length(theta.hat)
  tmp <- sum((theta.hat - mean(theta.hat))^2)
  sqrt(tmp/(nsim-1))
}

#' EmpSE MCSE
#' 
#' @param empse EMPSE
#' @param nsim Number of simulations
#' @return \eqn{empse/(sqrt(2*(nsim-1)))}
#' @keywords internal
#' 
empse.mcse <- function(empse, nsim) {
  empse/(sqrt(2*(nsim-1)))
} 

#' Variability ratio
#' 
#' @param theta.hat Theta hat
#' @param std.err Standard error
#' @return Ratio
#' @keywords internal
#' 
var.ratio <- function(theta.hat, std.err) {
  nsim <- length(theta.hat)
  num <- sum(std.err)/nsim
  denom <- sqrt(sum((theta.hat - mean(theta.hat))^2)/(nsim-1))
  num/denom
}

#' Variability ratio MCSE
#' 
#' Approximation of ratio variance based on independence of avg. se and emp.se
#' see \insertCite{wolter2007}{mimR}
#'
#' @param avg.se Average SE
#' @param emp.se Emp SE
#' @param var.avg.se Variance of average SE
#' @param var.emp.se Variance of Emp SE
#' 
#' @references
#' \insertRef{wolter2007}{mimR}
#' 
#' @keywords internal
#' 
var.ratio.mcse <- function(avg.se, emp.se, var.avg.se, var.emp.se) {

  sqrt((1/emp.se^2)*var.avg.se + (((avg.se^2)/(emp.se^4))*var.emp.se))
}

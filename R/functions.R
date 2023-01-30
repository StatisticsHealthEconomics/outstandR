# This file contains user-defined MAIC function for MAIC and functions to evaluate 
# the performance measures of interest

### MAIC functions

# Function to estimate MAIC weights 
maic <- function(X.EM) {
  # X.EM: centered S=1 effect modifiers
  # objective function to be minimized for standard method of moments MAIC
  Q <- function(alpha, X) {
    return(sum(exp(X %*% alpha)))
  }
  X.EM <- as.matrix(X.EM)
  N <- nrow(X.EM)
  K.EM <- ncol(X.EM)
  alpha <- rep(1,K.EM) # arbitrary starting point for the optimizer
  # objective function minimized using BFGS
  Q.min <- optim(fn=Q, X=X.EM, par=alpha, method="BFGS")
  hat.alpha <- Q.min$par # finite solution is the logistic regression parameters
  log.hat.w <- rep(0, N)
  for (k in 1:K.EM) {
    log.hat.w <- log.hat.w + hat.alpha[k]*X.EM[,k]
  }
  hat.w <- exp(log.hat.w) # estimated weights
  return(hat.w)
}

### Functions to evaluate performance measures
#
# bias estimate
bias <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  est <- sum(theta.hat)/nsim - theta
  return(est)
}

# Monte Carlo SE of bias estimate
bias.mcse <- function(theta.hat) {
  nsim <- length(theta.hat)
  tmp <- sum((theta.hat - mean(theta.hat))^2)
  mcse <- sqrt(1/(nsim*(nsim-1))*tmp)
  return(mcse)
}

# coverage estimate
coverage <- function(theta.hat.low, theta.hat.upp, theta) {
  nsim <- length(theta.hat.low)
  est <- sum(ifelse(theta>=theta.hat.low & theta<=theta.hat.upp,1,0))/nsim
  return(est)
}

# Monte Carlo SE of coverage estimate
coverage.mcse <- function(coverage, nsim) {
  mcse <- sqrt((coverage*(1-coverage))/nsim)
  return(mcse)
}

# MSE estimate
mse <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  est <- sum((theta.hat-theta)^2)/nsim
  return(est)
}

# Monte Carlo SE of MSE estimate
mse.mcse <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  tmp <- (theta.hat-theta)^2
  mse.est <- sum(tmp)/nsim
  mcse <- sqrt(sum((tmp - mse.est)^2)/(nsim*(nsim-1)))
  return(mcse)
}

# MAE estimate
mae <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  est <- sum(abs(theta.hat-theta))/nsim
  return(est)
}

# Monte Carlo SE of any continuous performance metric
mcse.estimate <- function(perf.measure) {
  nsim <- length(perf.measure)
  perf.measure.mean <- sum(perf.measure)/nsim
  mcse <- sqrt(sum((perf.measure-perf.measure.mean)^2)/(nsim*(nsim-1)))
  return(mcse)
}

# Empirical standard error 
empse <- function(theta.hat) {
  nsim <- length(theta.hat)
  tmp <- sum((theta.hat - mean(theta.hat))^2)
  est <- sqrt(tmp/(nsim-1))
  return(est)
}

# EmpSE MCSE
empse.mcse <- function(empse, nsim) {
  mcse <- empse/(sqrt(2*(nsim-1)))
  return(mcse)
} 

# Variability ratio
var.ratio <- function(theta.hat, std.err) {
  nsim <- length(theta.hat)
  num <- sum(std.err)/nsim
  denom <- sqrt(sum((theta.hat-mean(theta.hat))^2)/(nsim-1))
  est <- num/denom
  return(est)    
}

# Variability ratio MCSE
var.ratio.mcse <- function(avg.se, emp.se, var.avg.se, var.emp.se) {
  # approximation of ratio variance based on independence of avg. se and emp.se
  # see Wolter, K., 2007. Introduction to variance estimation. 
  mcse <- sqrt((1/emp.se^2)*var.avg.se + (((avg.se^2)/(emp.se^4))*var.emp.se))
  return(mcse)             
}

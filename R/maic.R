
#' Estimate MAIC weights 
#' 
#' @param X.EM centered S=1 effect modifiers
#' @return estimated weights
#' 
maic <- function(X.EM) {
  # objective function to be minimized for standard method of moments MAIC
  Q <- function(alpha, X) {
    sum(exp(X %*% alpha))
  }
  
  X.EM <- as.matrix(X.EM)
  N <- nrow(X.EM)
  K.EM <- ncol(X.EM)
  alpha <- rep(1, K.EM) # arbitrary starting point for the optimizer
  # objective function minimized using BFGS
  Q.min <- optim(fn=Q, X=X.EM, par=alpha, method="BFGS")
  hat.alpha <- Q.min$par # finite solution is the logistic regression parameters
  log.hat.w <- rep(0, N)
  for (k in seq_len(K.EM)) {
    log.hat.w <- log.hat.w + hat.alpha[k]*X.EM[,k]
  }
  
  exp(log.hat.w)
}


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


#' MAIC bootstrap
#' 
#' @return fitted treatment coefficient is marginal effect for A vs. C
#' 
maic.boot <- function(data, indices) {
  dat <- data[indices, ]      # AC bootstrap sample
  N <- nrow(dat)              # number of subjects in sample
  x.EM <- dat[,c("X1","X2")]  # AC effect modifiers
  # BC effect modifier means, assumed fixed
  theta <- BC.ALD[c("mean.X1", "mean.X2")]
  K.EM <- ncol(x.EM)          # number of effect modifiers
  # center the AC effect modifiers on the BC means
  x.EM$X1 <- x.EM$X1 - theta$mean.X1
  x.EM$X2 <- x.EM$X2 - theta$mean.X2
  # MAIC weight estimation using method of moments
  alpha <- rep(1,K.EM)        # arbitrary starting point for the optimizer
  # objective function minimized using BFGS
  Q.min <- optim(fn=Q, X.EM=as.matrix(x.EM), par=alpha, method="BFGS")
  # finite solution is the logistic regression parameters
  hat.alpha <- Q.min$par
  log.hat.w <- rep(0, N)
  for (k in seq_len(K.EM)) {
    log.hat.w <- log.hat.w + hat.alpha[k]*x.EM[,k]
  }
  hat.w <- exp(log.hat.w)      # estimated weights
  # fit weighted logistic regression model using glm
  outcome.fit <- glm(y ~ trt, family="quasibinomial",
                     weights=hat.w,
                     data=dat)
  
  coef(outcome.fit)["trt"]
}

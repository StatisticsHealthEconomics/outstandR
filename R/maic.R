
#' Estimate MAIC weights 
#' 
#' @param X.EM centred S=1 effect modifiers
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
  init <- rep(1, K.EM) # arbitrary starting point for the optimizer
  # objective function minimized using BFGS
  Q.min <- optim(fn=Q, X=X.EM, par=init, method="BFGS")
  
  # finite solution is the logistic regression parameters
  hat_alpha <- Q.min$par
  log.hat_w <- rep(0, N)
  
  for (k in seq_len(K.EM)) {
    log.hat_w <- log.hat_w + hat_alpha[k]*X.EM[, k]
  }
  
  exp(log.hat_w)
}


#' MAIC bootstrap
#' 
#' @param data original data
#' @param indices vector of indices which define the bootstrap sample
#' @return fitted treatment coefficient is marginal effect for A vs C
#' 
maic.boot <- function(data, indices) {
  dat <- data[indices, ]  # bootstrap sample
  X.EM <- dat[,c("X1","X2")]  # AC effect modifiers

  ##TODO: this seems odd. where is BC.ALD passed from?
  ##      the call uses AC.IPD data
  ##TODO: why is this centering used in maic.boot() and not maic()?
  # BC effect modifier means, assumed fixed
  theta <- BC.ALD[c("mean.X1", "mean.X2")]
  # centre AC effect modifiers on BC means
  X.EM$X1 <- X.EM$X1 - theta$mean.X1
  X.EM$X2 <- X.EM$X2 - theta$mean.X2
  
  hat_w <- maic(X.EM)
  
  # fit weighted logistic regression model
  outcome.fit <- glm(y ~ trt,
                     family = "quasibinomial",
                     weights = hat_w,
                     data = dat)
  
  coef(outcome.fit)["trt"]
}

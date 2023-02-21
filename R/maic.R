
#' Estimate MAIC weights 
#' 
#' @references J. E. Signorovitch et al,
#' Comparative Effectiveness Without Head-to-Head Trials:
#' A Method for Matching-Adjusted Indirect Comparisons Applied
#' to Psoriasis Treatment with Adalimumab or Etanercept,
#' Pharmacoeconomics 2010; 28 (10): 935-945
#' 
#' @param X.EM Centred S=1 effect modifiers; matrix or data frame
#' @return Estimated weights for each individual; vector
#' 
maic_weights <- function(X.EM) {
  # objective function to minimize for standard method of moments MAIC
  Q <- function(beta, X) {
    sum(exp(X %*% beta))
  }
  
  X.EM <- as.matrix(X.EM)
  
  N <- nrow(X.EM)    # number of individuals
  K <- ncol(X.EM)    # number of covariates
  init <- rep(1, K)  # arbitrary starting point for optimizer
  Q.min <- optim(fn=Q, X=X.EM, par=init, method="BFGS")
  
  # finite solution is the logistic regression parameters
  hat_beta <- Q.min$par
  log.hat_w <- rep(0, N)
  
  # linear eqn for logistic
  for (k in seq_len(K)) {
    log.hat_w <- log.hat_w + hat_beta[k]*X.EM[, k]
  }
  
  exp(log.hat_w)
}


#' MAIC bootstrap
#' 
#' @param data original data
#' @param indices vector of indices which define the bootstrap sample
#' @return fitted treatment coefficient is marginal effect for A vs C
#' 
maic.boot <- function(data, indices, formula) {
  dat <- data[indices, ]  # bootstrap sample
  
  # effect_modifier_names <- get_effect_modifiers(formula)
  # X.EM <- dat[, effect_modifier_names]
  
  X.EM <- dat[, c("X1","X2")]  # AC effect modifiers
  
  ##TODO: this seems odd. where is BC.ALD passed from?
  ##      the call uses AC.IPD data
  ##TODO: why is this centering used in maic.boot() and not maic()?
  # BC effect modifier means, assumed fixed
  theta <- dat[c("mean.X1", "mean.X2")]
  # centre AC effect modifiers on BC means
  X.EM$X1 <- X.EM$X1 - theta$mean.X1
  X.EM$X2 <- X.EM$X2 - theta$mean.X2
  
  # X.EM <- centre_effect_modifiers(X.EM, formula)
  
  hat_w <- maic_weights(X.EM)
  
  # fit weighted logistic regression model
  outcome.fit <- glm(formula,
                     family = "quasibinomial",
                     weights = hat_w,
                     data = dat)
  
  coef(outcome.fit)["trt"]
}



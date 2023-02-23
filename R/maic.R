
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
#' @param indices vector of indices, same length as original,
#'   which define the bootstrap sample
#' @return fitted treatment coefficient is marginal effect for A vs C
#' 
maic.boot <- function(data, indices, formula, dat_ALD) {
  dat <- data[indices, ]  # bootstrap sample
  
  effect_modifier_names <- get_effect_modifiers(formula)
  X.EM <- dat[, effect_modifier_names]
  
  ##TODO: why is this centering used in maic.boot() and not maic()?
  
  browser()
  # BC effect modifier means, assumed fixed
  mean_names <- get_mean_names(dat_ALD, effect_modifier_names)

  # centre AC effect modifiers on BC means
  X.EM <- X.EM - dat_ALD[, mean_names]
 
  hat_w <- maic_weights(X.EM)
  
  treat_nm <- get_treatment_name(formula)
  formula_treat <- glue::glue("{formula[[2]]} ~ {treat_nm}")
  
  # fit weighted logistic regression model
  fit <- glm(formula_treat,
             family = "quasibinomial",
             weights = hat_w,
             data = dat)
  
  coef(fit)[treat_nm]
}



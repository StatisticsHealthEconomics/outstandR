
#' Estimate MAIC weights 
#' 
#' Matching-adjusted indirect comparison weights.
#' Method is taken from
#' \insertCite{Signorovitch2010}{mimR}.
#' 
#' @param X_EM Centred \eqn{S=1} effect modifiers; matrix or data frame
#' @return Estimated weights for each individual; vector
#' 
#' @references
#' \insertRef{Signorovitch2010}{mimR}
#' 
maic_weights <- function(X_EM) {
  X_EM <- as.matrix(X_EM)
  
  N <- nrow(X_EM)    # number of individuals
  K <- ncol(X_EM)    # number of covariates
  init <- rep(1, K)  # arbitrary starting point for optimizer
  Q.min <- optim(fn=Q, X=X_EM, par=init, method="BFGS")
  
  # finite solution is the logistic regression parameters
  hat_beta <- Q.min$par
  log.hat_w <- rep(0, N)
  
  # linear equation for logistic
  for (k in seq_len(K)) {
    log.hat_w <- log.hat_w + hat_beta[k]*X_EM[, k]
  }
  
  exp(log.hat_w)
}


#' MAIC bootstrap sample
#' 
#' Matching-adjusted indirect comparison bootstrap sampling.
#' 
#' @template args-ipd
#' @param indices Vector of indices, same length as original,
#'   which define the bootstrap sample
#' @param formula Linear regression formula
#' @template args-ald
#' @return Fitted treatment coefficient is marginal effect for _A_ vs _C_
#' @seealso [IPD_stats.maic()]
#' 
maic.boot <- function(ipd, indices, formula, ald) {
  
  dat <- ipd[indices, ]  # bootstrap sample
  
  effect_modifier_names <- get_effect_modifiers(formula)
  X_EM <- dat[, effect_modifier_names]
  
  ##TODO: why is this centering used in maic.boot() and not maic()?
  
  # BC effect modifier means, assumed fixed
  mean_names <- get_mean_names(ald, effect_modifier_names)

  # centre AC effect modifiers on BC means
  dat_ALD_means <- ald[, mean_names][rep(1, nrow(X_EM)), ]
  X_EM <- X_EM - dat_ALD_means
 
  hat_w <- maic_weights(X_EM)
  
  treat_nm <- get_treatment_name(formula)
  formula_treat <- glue::glue("{formula[[2]]} ~ {treat_nm}")

  # fit weighted logistic regression model
  fit <- glm(formula_treat,
             family = "quasibinomial",
             weights = hat_w,
             data = cbind(dat, hat_w = hat_w))
  
  coef(fit)[treat_nm]
}



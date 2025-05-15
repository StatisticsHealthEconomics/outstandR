
#' Estimate MAIC weights 
#' 
#' Matching-adjusted indirect comparison weights.
#' Method is taken from
#' \insertCite{Signorovitch2010}{outstandR}.
#' 
#' @param X_EM Centred \eqn{S=1} effect modifiers; matrix or data frame
#' @return Estimated weights for each individual; vector
#' @importFrom stats optim
#' @references
#' \insertRef{Signorovitch2010}{outstandR}
#' @keywords internal
#' 
maic_weights <- function(X_EM) {
  X_EM <- as.matrix(X_EM)
  
  N <- nrow(X_EM)    # number of individuals
  K <- ncol(X_EM)    # number of covariates
  init <- rep(1, K)  # arbitrary starting point for optimizer
  
  ##TODO: what about scaling X_EM?
  ##      because large values return error
  
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

#' Objective function to minimize for standard method of moments MAIC
#'
#' @param beta Beta coefficient
#' @param X Covariate value
#' @keywords internal
#' 
Q <- function(beta, X) {
  sum(exp(X %*% beta))
}

#' MAIC bootstrap sample
#' 
#' Matching-adjusted indirect comparison bootstrap sampling.
#' 
#' @eval study_data_args(include_ipd = TRUE, include_ald = TRUE)
#' @param indices Vector of indices, same length as original,
#'   which define the bootstrap sample
#' @eval reg_args(include_formula = TRUE, include_family = TRUE)
#' @param hat_w MAIC weights; default `NULL` which calls [maic_weights()]
#' 
#' @return Vector of fitted probabilities for treatments _A_ and _C_
#' @seealso [IPD_stats.maic()]
#' @keywords internal
#' 
maic.boot <- function(ipd, indices = 1:nrow(ipd),
                      formula, family, ald,
                      trt_var,
                      hat_w = NULL) {
  
  dat <- ipd[indices, ]  # bootstrap sample
  n_ipd <- length(indices)
  
  effect_modifier_names <- get_eff_mod_names(formula)
  
  X_EM <- dat[, effect_modifier_names]
  
  # centre AC effect modifiers on BC means
  dat_ALD_means <- ald |> 
    dplyr::filter(variable %in% effect_modifier_names,
                  statistic == "mean") |> 
    tidyr::pivot_wider(names_from = variable) |> 
    dplyr::select(all_of(effect_modifier_names)) |> 
    tidyr::uncount(weights = n_ipd)
  
  centred_EM <- X_EM - dat_ALD_means
  
  if (is.null(hat_w)) {
    hat_w <- maic_weights(centred_EM)
  }
  
  formula_treat <- glue::glue("{formula[[2]]} ~ {trt_var}")
  
  # so can use non-integer weights
  if (family$family == "binomial") {
    family <- quasibinomial()
  }
  
  # fit weighted logistic regression model
  fit <- glm(formula = formula_treat,
             family = family,
             weights = hat_w / mean(hat_w),
             data = cbind(dat, hat_w = hat_w))
  
  # extract model coefficients
  coef_fit <- coef(fit)
  
  # probabilities using inverse link
  linkinv <- family$linkinv
  
  pC <- unname(linkinv(coef_fit[1]))                # probability for control group
  pA <- unname(linkinv(coef_fit[1] + coef_fit[2]))  # probability for treatment group
  
  c(pC = pC, pA = pA)
}


#' @export
#' @importFrom boot boot
#' 
calc_maic <- function(strategy,
                      ipd, ald) {
  args_list <- 
    list(R = strategy$R,
         formula = strategy$formula,
         family = strategy$family,
         trt_var = strategy$trt_var,
         data = ipd,
         ald = ald)
  
  maic_boot <- do.call(boot::boot, c(statistic = maic.boot, args_list))
  
  list(mean_A = maic_boot$t[, 2],
       mean_C = maic_boot$t[, 1])  
}

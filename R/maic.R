
#' Estimate MAIC weights 
#' 
#' Matching-adjusted indirect comparison weights.
#' Method is taken from
#' \insertCite{Signorovitch2010}{outstandR}.
#' 
#' @param X_EM Centred \eqn{S=1} effect modifiers IPD covariates; matrix or data frame
#' @return Estimated weights for each individual; vector
#' 
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
  
  # find betas
  Q.min <- optim(fn = Q, X = X_EM, par = init, method = "BFGS")
  
  # check for convergence issues
  if (Q.min$convergence != 0) {
    warning(
      paste("optim did not converge (code:", Q.min$convergence, "). Message:", Q.min$message))
  }
  
  # finite solution is the logistic regression parameters
  hat_beta <- Q.min$par
  
  # Calculate the log weights
  # log(w_i) = sum(beta_k * X_ik)
  log.hat_w <- as.vector(X_EM %*% hat_beta)
  
  exp(log.hat_w)
}

#' Objective function to minimize for standard method of moments MAIC
#'
#' @param beta Beta coefficient to find
#' @param X Covariate value matrix, centred
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
  n_trts <- length(unique(dat[[trt_var]]))
  
  # ensure bootstrap sample contains more than one treatment level
  if (n_trts < 2) {
    warning("Bootstrap sample contains less than two treatment levels. Returning NA.")
    return(c(pC = NA, pA = NA))
  }
  
  effect_modifier_names <- get_eff_mod_names(formula, trt_var)
  
  if (length(effect_modifier_names) > 0) {

    X_EM_prepared <- matrix(NA, nrow = n_ipd, ncol = length(effect_modifier_names))
    colnames(X_EM_prepared) <- effect_modifier_names
    
    # determine covariate types
    for (em_name in effect_modifier_names) {
      ipd_col <- dat[[em_name]]
      
      # Attempt to get 'mean' from ALD. If not found, try 'prop'
      # This assumes ALD contains either 'mean' or 'proportion' for each effect modifier
      ald_mean_val <- ald |>
        dplyr::filter(variable == em_name, statistic == "mean") |>
        dplyr::pull(value)
      
      ald_prop_val <- ald |>
        dplyr::filter(variable == em_name, statistic == "prop") |>
        dplyr::pull(value)
      
      if (length(ald_mean_val) > 0) {
        # continuous if 'mean' in ALD
        X_EM_prepared[, em_name] <- ipd_col - ald_mean_val
        # scaling continuous variables by their standard deviation can improve optimizer performance
        col_sd <- sd(X_EM_prepared[, em_name])
        
        if (col_sd > 0) {
          X_EM_prepared[, em_name] <- X_EM_prepared[, em_name] / col_sd
        }
      } else if (length(ald_prop_val) > 0) {
        # binary if 'prop' in ALD
        # assumes binary variables are coded 0/1 in IPD
        X_EM_prepared[, em_name] <- ipd_col - ald_prop_val
      } else {
        stop(paste("Neither 'mean' nor 'prop' found in ALD for covariate:", em_name))
      }
    }
    
    # calculate MAIC weights if not provided
    if (is.null(hat_w)) {
      hat_w <- maic_weights(X_EM = X_EM_prepared)
    }
  } else {
    # if no covariates, all weights are 1 (unadjusted comparison)
    hat_w <- rep(1, n_ipd)
  }
  
  formula_treat <- glue::glue("{formula[[2]]} ~ {trt_var}")
  
  # so can use non-integer weights
  if (family$family == "binomial") {
    family <- quasibinomial()
  }
  
  # fit weighted regression model
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
                      analysis_params) {
  args_list <- 
    list(R = strategy$R,
         formula = strategy$formula,
         family = strategy$family,
         trt_var = strategy$trt_var,
         data = analysis_params$ipd,
         ald = analysis_params$ald)
  
  maic_boot <- do.call(boot::boot, c(statistic = maic.boot, args_list))
  
  list(mean_A = maic_boot$t[, 2],
       mean_C = maic_boot$t[, 1])  
}

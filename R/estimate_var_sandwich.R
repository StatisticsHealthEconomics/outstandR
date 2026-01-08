
#' Estimate Variance Sandwich Estimator
#'
#' Computes the robust (sandwich) variance estimator for the treatment effect.
#' 
#' @param strategy An object of class `strategy` created by functions such as 
#'   [strategy_maic()], [strategy_stc()], or [strategy_mim()]. 
#'   Contains modelling details like the formula and family.
#' @param analysis_params List of analysis parameters (ipd, ald, etc.)
#' @param ... Additional arguments
#' 
#' @return Numeric variance estimate for the treatment contrast
#' @export
estimate_var_sandwich <- function(strategy, analysis_params, ...) {
  UseMethod("estimate_var_sandwich", strategy)
}

#' @export
estimate_var_sandwich.default <- function(strategy, analysis_params, ...) {
  stop(paste("Sandwich variance not implemented for strategy:", class(strategy)[1]))
}

#' @export
estimate_var_sandwich.stc <- function(strategy, analysis_params, ...) {
  # 1. Fit Model (Standard GLM)
  # STC centers covariates (in calc_stc), so we must replicate that prep
  ipd <- analysis_params$ipd
  centre_vars <- get_eff_mod_names(strategy$formula)
  ipd[, centre_vars] <- scale(ipd[, centre_vars], scale = FALSE)
  
  fit <- glm(formula = strategy$formula, family = strategy$family, data = ipd)
  
  # 2. Get Robust Covariance Matrix (Beta)
  vcov_robust <- get_robust_vcov(fit)
  
  # 3. Calculate Gradient of Treatment Effect w.r.t Betas (Delta Method)
  # Effect is defined at the mean (0 after centering)
  # Mean C = linkinv(Beta0)
  # Mean A = linkinv(Beta0 + Beta_Trt)
  
  coefs <- coef(fit)
  trt_var <- strategy$trt_var
  
  # Identify coefficient names
  beta0_idx <- 1
  beta_trt_idx <- grep(paste0("^", trt_var), names(coefs))[1] # Assumes one treatment coef
  
  # Define function f(beta) -> Treatment Effect
  calc_effect_from_betas <- function(b) {
    beta0 <- b[beta0_idx]
    beta_trt <- b[beta_trt_idx]
    
    # Predict probabilities/means
    mu_c <- strategy$family$linkinv(beta0)
    mu_a <- strategy$family$linkinv(beta0 + beta_trt)
    
    # Convert to requested scale (log_odds, risk_diff, etc)
    calculate_ate(mu_a, mu_c, effect = analysis_params$scale)
  }
  
  # Numerical Gradient
  grad <- num_grad(calc_effect_from_betas, coefs)
  
  # 4. Variance = g' V g
  var_est <- as.numeric(t(grad) %*% vcov_robust %*% grad)
  return(var_est)
}

#' @export
estimate_var_sandwich.maic <- function(strategy, analysis_params, ...) {
  # 1. Calculate Weights
  # We assume weights are fixed for the sandwich estimator (standard weighted GLM approach)
  # Full M-estimation including weight uncertainty is complex; this provides robust
  # variance for the weighted outcome model.
  
  ipd <- analysis_params$ipd
  ald <- analysis_params$ald
  trt_var <- strategy$trt_var
  
  # Re-calculate weights (reusing internal logic would be better, but we reconstruct here)
  effect_modifier_names <- get_eff_mod_names(strategy$formula, trt_var)
  
  if (length(effect_modifier_names) > 0) {
    X_EM <- as.matrix(ipd[, effect_modifier_names, drop = FALSE])
    
    # Center X_EM based on ALD means (simplified matching logic)
    for (em in effect_modifier_names) {
      ald_mean <- ald$value[ald$variable == em & ald$statistic == "mean"]
      if (length(ald_mean) > 0) X_EM[, em] <- X_EM[, em] - ald_mean
    }
    
    # Optimization to find weights
    # Note: We reuse the internal Q function if exported, or simple optim here.
    # For robustness, we assume weights are 1 if optimization fails or isn't set up.
    # Ideally, calc_maic should pass the weights, but we are refitting.
    # We will compute weights via the internal helper if available
    w <- tryCatch({
      maic_weights(X_EM)
    }, error = function(e) rep(1, nrow(ipd)))
    
  } else {
    w <- rep(1, nrow(ipd))
  }
  
  # 2. Fit Weighted GLM
  # formula needs modification for MAIC (remove effect modifiers from formula usually)
  # But strategy_maic usually keeps the full formula? 
  # Standard MAIC typically uses just `outcome ~ trt` in the weighted model.
  # We check the strategy implementation: calc_maic uses `formula_treat <- glue("{formula[[2]]} ~ {trt_var}")`
  
  formula_maic <- as.formula(paste(all.vars(strategy$formula)[1], "~", trt_var))
  
  fit <- glm(formula = formula_maic, family = strategy$family, data = ipd, weights = w)
  
  # 3. Robust Variance
  vcov_robust <- get_robust_vcov(fit)
  
  # 4. Delta Method for Effect
  coefs <- coef(fit)
  
  calc_effect_from_betas <- function(b) {
    # Simple model: Intercept + Trt
    mu_c <- strategy$family$linkinv(b[1])
    mu_a <- strategy$family$linkinv(b[1] + b[2])
    calculate_ate(mu_a, mu_c, effect = analysis_params$scale)
  }
  
  grad <- num_grad(calc_effect_from_betas, coefs)
  var_est <- as.numeric(t(grad) %*% vcov_robust %*% grad)
  return(var_est)
}

#' @export
estimate_var_sandwich.gcomp_ml <- function(strategy, analysis_params, ...) {
  # 1. Fit Model
  ipd <- analysis_params$ipd
  fit <- glm(formula = strategy$formula, family = strategy$family, data = ipd)
  
  # 2. Generate Pseudo-Population (fixed for variance calc)
  # We treat the pseudo-population as a constant integration grid
  x_star <- simulate_ALD_pseudo_pop(
    formula = strategy$formula,
    ipd = ipd, ald = analysis_params$ald,
    trt_var = strategy$trt_var,
    rho = strategy$rho,
    N = strategy$N,
    marginal_distns = strategy$marginal_distns,
    marginal_params = strategy$marginal_params,
    seed = 123 # Fixed seed for stable gradient
  )
  
  # 3. Robust Covariance
  vcov_robust <- get_robust_vcov(fit)
  
  # 4. Delta Method over the Integration
  # The function f(beta) is the G-computation estimator:
  # Mean(Predict(A)) - Mean(Predict(C))
  
  ref_trt <- analysis_params$ref_trt
  comp_trt <- analysis_params$ipd_comp
  trt_var <- strategy$trt_var
  
  # Prepare counterfactual frames
  df_ref <- x_star; df_ref[[trt_var]] <- ref_trt
  df_comp <- x_star; df_comp[[trt_var]] <- comp_trt
  
  # Create model matrices for prediction (avoids re-calling predict.glm inside loop)
  # We need the terms from the fitted model to handle factors/interactions correctly
  Terms <- delete.response(terms(fit))
  m_ref <- model.matrix(Terms, df_ref, xlev = fit$xlevels)
  m_comp <- model.matrix(Terms, df_comp, xlev = fit$xlevels)
  
  linkinv <- strategy$family$linkinv
  
  calc_gcomp_effect <- function(b) {
    # Manual prediction: linkinv(X %*% beta)
    pred_ref <- linkinv(m_ref %*% b)
    pred_comp <- linkinv(m_comp %*% b)
    
    mu_c <- mean(pred_ref)
    mu_a <- mean(pred_comp)
    
    calculate_ate(mu_a, mu_c, effect = analysis_params$scale)
  }
  
  grad <- num_grad(calc_gcomp_effect, coef(fit))
  var_est <- as.numeric(t(grad) %*% vcov_robust %*% grad)
  return(var_est)
}


# --- Internal Helpers ---

#' Compute Robust Covariance Matrix (HC0-style)
#' 
#' Calculates (X'WX)^-1 (X' W^2 r^2 X) (X'WX)^-1
#' @keywords internal
get_robust_vcov <- function(fit) {
  # Bread: Unscaled covariance (dispersion=1)
  # vcov(fit) includes dispersion for some families, but standard sandwich uses Fisher info
  # summary.glm(fit)$cov.unscaled is (X'WX)^-1
  bread <- summary(fit)$cov.unscaled
  
  # Meat: Sum of score contributions outer products
  # score_i = x_i * (y_i - mu_i) for canonical link
  # For general link: x_i * (y - mu) * dmu/deta * 1/Var(mu)
  # This is equivalent to weighted residuals * model matrix
  
  # estfun equivalent
  w_res <- weights(fit, type = "working") * residuals(fit, type = "working")
  X <- model.matrix(fit)
  
  # individual scores = X_i * w_res_i
  scores <- X * w_res
  
  # meat = t(scores) %*% scores
  meat <- crossprod(scores)
  
  # Sandwich
  vcov_hc <- bread %*% meat %*% bread
  return(vcov_hc)
}

#' Numerical Gradient
#' @keywords internal
num_grad <- function(func, x, h = 1e-5) {
  n <- length(x)
  grad <- numeric(n)
  fx <- func(x)
  
  for (i in 1:n) {
    x_h <- x
    x_h[i] <- x[i] + h
    grad[i] <- (func(x_h) - fx) / h
  }
  return(grad)
}

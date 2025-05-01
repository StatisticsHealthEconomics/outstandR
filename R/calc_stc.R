
#' Calculate simulated treatment comparison statistics
#' @importFrom stats glm
#' @export
#' 
calc_stc <- function(strategy, ipd, ...) {
  
  treat_nm <- strategy$trt_var
  
  # centre covariates
  centre_vars <- get_eff_mod_names(strategy$formula)
  
  ipd[, centre_vars] <- scale(ipd[, centre_vars], scale = FALSE)
  
  fit <- glm(formula = strategy$formula,
             family = strategy$family,
             data = ipd)
  
  # extract model coefficients
  coef_fit <- coef(fit)
  
  # safer than treat_nm in case of factor level append
  treat_coef_name <- grep(paste0("^", treat_nm), names(coef_fit), value = TRUE)
  
  # probability for control and treatment group
  # estimating treatment effect at means because of centring
  linkinv <- strategy$family$linkinv
  mean_C <- linkinv(coef_fit[1])
  mean_A <- linkinv(coef_fit[1] + coef_fit[treat_coef_name])
  
  list(mean_A = mean_A,
       mean_C = mean_C)
}

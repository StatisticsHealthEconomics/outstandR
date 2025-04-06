
#' Calculate simulated treatment comparison statistics
#' @importFrom stats glm
#' @export
#' 
calc_stc <- function(strategy, ipd, ...) {
  
  # centre covariates
  centre_vars <- get_eff_mod_names(strategy$formula)
  
  ipd[, centre_vars] <- scale(ipd[, centre_vars], scale = FALSE)
  
  fit <- glm(formula = strategy$formula,
             family = strategy$family,
             data = ipd)
  
  # extract model coefficients
  coef_fit <- coef(fit)
  
  treat_nm <- get_treatment_name(strategy$formula)
  
  # safer than treat_nm in case of factor level append
  treat_coef_name <- grep(paste0("^", treat_nm), names(coef_fit), value = TRUE)
  
  # probability for control and treatment group
  # estimating treatment effect at means because of centring
  mean_C <- strategy$family$linkinv(coef_fit[1])
  mean_A <- strategy$family$linkinv(coef_fit[1] + coef_fit[treat_coef_name])
  
  list(mean_A = mean_A,
       mean_C = mean_C)
}

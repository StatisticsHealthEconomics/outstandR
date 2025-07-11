
#' Calculate simulated treatment comparison statistics
#' @return A list:
#' \describe{
#'   \item{`mean_A`}{Mean for comparator treatment group "A".}
#'   \item{`mean_C`}{Mean for reference treatment group "C".}
#' }
#' @importFrom stats glm
#' @export
#' 
calc_stc <- function(strategy, analysis_params, ...) {
  
  ipd <- analysis_params$ipd
  trt_var <- strategy$trt_var
  
  # centre covariates
  centre_vars <- get_eff_mod_names(strategy$formula)
  
  ipd[, centre_vars] <- scale(ipd[, centre_vars], scale = FALSE)
  
  fit <- glm(formula = strategy$formula,
             family = strategy$family,
             data = ipd)
  
  # extract model coefficients
  coef_fit <- coef(fit)
  
  # safer than trt_var in case of factor level append
  treat_coef_name <-
    grep(pattern = paste0("^", trt_var, "[^:]*$"),
         names(coef_fit), value = TRUE)
  
  # probability for control and treatment group
  # estimating treatment effect at means because of centring
  linkinv <- strategy$family$linkinv
  mean_C <- linkinv(coef_fit[1])
  mean_A <- linkinv(coef_fit[1] + coef_fit[treat_coef_name])
  
  list(mean_A = mean_A,
       mean_C = mean_C)
}

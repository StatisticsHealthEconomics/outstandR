
#' Calculate simulated treatment comparison statistics
#'
#' @param strategy An object of class `strategy` created by functions such as
#'   [strategy_maic()], [strategy_stc()], or [strategy_mim()].
#'   Contains modelling details like the formula and family.
#' @param analysis_params List of analysis parameters. Must contain `ipd`
#'   (individual patient data).
#' @param ... Additional arguments.
#'
#' @return A list containing:
#' * `means`: A list containing:
#'     * `A`: Mean for comparator treatment group "A".
#'     * `C`: Mean for reference treatment group "C".
#' * `model`: The fitted [stats::glm()] object.
#'
#' @importFrom stats glm coef
#' @keywords internal
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
  # estimating treatment effect at means because of centering
  linkinv <- strategy$family$linkinv
  mean_C <- linkinv(coef_fit[1])
  mean_A <- linkinv(coef_fit[1] + coef_fit[treat_coef_name])
  
  list(
    means = list(
      A = mean_A,
      C = mean_C),
    model = list(
      fit = fit))
}

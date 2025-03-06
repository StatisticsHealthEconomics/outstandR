
#' Calculate simulated treatment comparison statistics
#' 
calc_stc <- function(strategy, ipd, ...) {
  
  # centre covariates
  term.labels <- attr(terms(strategy$formula), "term.labels")
  centre_vars <- gsub("trt:", "", term.labels[grepl(":", term.labels)])
  
  ipd[, centre_vars] <- scale(ipd[, centre_vars], scale = FALSE)
  
  fit <- glm(formula = strategy$formula,
             family = strategy$family,
             data = ipd)
  
  # extract model coefficients
  coef_fit <- coef(fit)
  
  treat_nm <- get_treatment_name(strategy$formula)
  
  # probability for control and treatment group
  mean_A <- strategy$family$linkinv(coef_fit[1])
  mean_C <- strategy$family$linkinv(coef_fit[1] + coef_fit[treat_nm])
  
  list(mean_A = mean_A,
       mean_C = mean_C)
}

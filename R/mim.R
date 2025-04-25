
#' Multiple imputation marginalization (MIM)
#' 
#' @param strategy Strategy object
#' @template args-ipd
#' @template args-ald
#' @param ... Argument to pass to Stan model
#' 
#' @importFrom rstanarm posterior_predict stan_glm
#' @keywords internal
#' 
calc_mim <- function(strategy,
                     ipd, ald, ...) {
  
  formula <- strategy$formula
  family <- strategy$family
  rho <- strategy$rho
  N <- strategy$N
  
  x_star <- simulate_ALD_pseudo_pop(formula, ipd, ald, rho, N)
  
  ## SYNTHESIS STAGE ##
  
  # first-stage logistic regression model fitted to index RCT using MCMC (Stan)
  outcome.model <- stan_glm(
    formula = formula,
    data = ipd,
    family = family,
    algorithm = "sampling", ...)
  
  # create augmented target dataset
  target.t1 <- target.t0 <- x_star
  target.t1$trt <- 1
  target.t0$trt <- 0
  
  aug.target <- rbind(target.t0, target.t1)
  
  # complete syntheses by drawing binary outcomes
  # from their posterior predictive distribution
  y_star <-
    rstanarm::posterior_predict(
      outcome.model, newdata = aug.target)
  
  ## ANALYSIS STAGE ##
  
  M <- nrow(y_star)
  
  # fit second-stage regression to each synthesis using maximum-likelihood estimation
  reg2.fits <- lapply(1:M, function(m) {
    glm(y_star[m, ] ~ trt, data = aug.target, family = family)
  })
  
  # treatment effect point estimates in each synthesis
  coef_fit <- do.call(rbind, lapply(reg2.fits, function(fit) coef(fit)))
  
  ##TODO: how to transform this to the prob scale?
  # point estimates for the variance in each synthesis
  hats.v <- unlist(lapply(reg2.fits,
                          function(fit)
                            vcov(fit)["trt", "trt"]))
  
  mean_C <- family$linkinv(coef_fit[, 1])                  # probability for control
  mean_A <- family$linkinv(coef_fit[, 1] + coef_fit[, 2])  # probability for treatment
  
  tibble::lst(mean_A, mean_C,
              hats.v, M)
}

#' Wald-type interval estimates
#' 
#' Constructed using t-distribution with nu degrees of freedom
#'
#' @param M Number of syntheses used in analysis stage (high for low Monte Carlo error)
#' @param bar.v "within" variance (average of variance point estimates)
#' @param b "between" variance (sample variance of point estimates)
#'
#' @keywords internal
#' 
wald_type_interval <- function(M, bar.v, b) {
  (M - 1) * (1 + bar.v / ((1 + 1 / M) * b)) ^ 2
}

#' Variance estimate by pooling
#' 
#' Use combining rules to estimate
#' 2003
#' 
#' @param M Number of syntheses used in analysis stage (high for low Monte Carlo error)
#' @param bar.v "within" variance (average of variance point estimates)
#' @param b "between" variance (sample variance of point estimates)
#' 
#' @keywords internal
#' 
var_by_pooling <- function(M, bar.v, b) {
  (1 + (1 / M)) * b - bar.v
}


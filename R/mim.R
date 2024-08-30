
#' Multiple imputation marginalization (MIM)
#' 
#' @param formuala
#' @template args-ipd
#' @template args-ald
#' @param M Number of syntheses used in analysis stage (high for low Monte Carlo error)
#' @importFrom rstanarm posterior_predict stan_glm
#' @keywords internal
#' 
mim <- function(formula,
                ipd, ald,
                M = 1000,
                n.chains = 2,
                warmup = 1000,
                iters = 4000) {
  
  if (!inherits(formula, "formula"))
    stop("formula argument must be of formula class.")
  
  x_star <- simulate_ALD_pseudo_pop(formula, ipd, ald)
  
  ## SYNTHESIS STAGE ##
  
  # first-stage logistic regression model fitted to index RCT using MCMC (Stan)
  outcome.model <- stan_glm(
    formula,
    data = ipd,
    family = binomial,
    algorithm = "sampling",
    iter = iters,
    warmup = warmup,
    chains = n.chains,
    thin = (n.chains * (iters - warmup)) / M)
  
  # create augmented target dataset
  target.t1 <- target.t0 <- x_star
  target.t1$trt <- 1
  target.t0$trt <- 0
  
  aug.target <- rbind(target.t0, target.t1)
  
  # complete syntheses by drawing binary outcomes
  # from their posterior predictive distribution
  y_star <- rstanarm::posterior_predict(outcome.model, newdata = aug.target)
  
  ## ANALYSIS STAGE ##
  
  # fit second-stage regression to each synthesis using maximum-likelihood estimation
  reg2.fits <- lapply(1:M, function(m)
    glm(y_star[m, ] ~ trt, data = aug.target, family = binomial))
  
  # treatment effect point estimates in each synthesis
  hats.delta <- unlist(lapply(reg2.fits,
                              function(fit)
                                coef(fit)["trt"][[1]]))
  
  # point estimates for the variance in each synthesis
  hats.v <- unlist(lapply(reg2.fits,
                          function(fit)
                            vcov(fit)["trt", "trt"]))
  
  tibble::lst(hats.delta, hats.v, M)
}

#' Wald-type interval estimates
#' Constructed using t-distribution with nu degrees of freedom
#' @keywords internal
#' 
wald_type_interval <- function(M, bar.v, b) {
  (M - 1) * (1 + bar.v / ((1 + 1 / M) * b)) ^ 2
}

#' Variance estimate by pooling
#' Use combining rules to estimate
#' @keywords internal
#' 
var_by_pooling <- function(M, bar.v, b) {
  (1 + (1 / M)) * b - bar.v
}


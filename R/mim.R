
#' Multiple imputation marginalization (MIM)
#' 
#' @param ipd.index,ipd.target Index RCT and target covariate datasets
#' @param M Number of syntheses used in analysis stage (high for low Monte Carlo error)
#'
mim <- function(ipd.index,
                ipd.target,
                M = 1,
                n.chains = 2,
                warmup = 10,
                iters = 1000) {
  
  ## SYNTHESIS STAGE ##
  
  # first-stage logistic regression model fitted to index RCT using MCMC (Stan)
  outcome.model <- stan_glm(
    y ~ (trt * x1) + (trt * x2),
    data = ipd.index,
    family = binomial,
    algorithm = "sampling",
    iter = iters,
    warmup = warmup,
    chains = n.chains,
    # thin to use M independent samples in analysis stage
    thin = (n.chains * (iters - warmup)) / M
  )
  
  # create augmented target dataset
  target.t1 <- target.t0 <- ipd.target
  target.t1$trt <- 1
  target.t0$trt <- 0
  
  aug.target <- rbind(target.t0, target.t1)
  
  # complete syntheses by drawing binary outcomes from their posterior predictive distribution
  y.star <- posterior_predict(outcome.model, newdata = aug.target)
  
  ## ANALYSIS STAGE ##
  
  # fit second-stage regression to each synthesis using maximum-likelihood estimation
  reg2.fits <- lapply(1:M, function(m)
    glm(y.star[m, ] ~ trt, data = aug.target, family = binomial))
  
  # treatment effect point estimates in each synthesis
  hats.delta <- unlist(lapply(reg2.fits, function(fit)
    coef(fit)["trt"][[1]]))
  
  # point estimates for the variance in each synthesis
  hats.v <- unlist(lapply(reg2.fits, function(fit)
    vcov(fit)["trt", "trt"]))
  
  # quantities originally defined by Rubin (1987) for multiple imputation
  bar.delta <- mean(hats.delta)  # average of treatment effect point estimates
  bar.v <- mean(hats.v)          # "within" variance (average of variance point estimates)
  b <- var(hats.delta)           # "between" variance (sample variance of point estimates)
  
  # pooling: average of point estimates is marginal log odds ratio
  hat.Delta <- bar.delta
  
  # pooling: use combining rules to estimate the variance
  hat.var.Delta <- (1 + (1 / M)) * b - bar.v
  
  # Wald-type interval estimates constructed using t-distribution with nu degrees of freedom
  nu <- (M - 1) * (1 + bar.v / ((1 + 1 / M) * b)) ^ 2
  
  lci.Delta <- hat.Delta + qt(0.025, df = nu) * sqrt(hat.var.Delta)
  uci.Delta <- hat.Delta + qt(0.975, df = nu) * sqrt(hat.var.Delta)
  
  hat.Delta
}

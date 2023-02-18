
#'
gcomp_stan <- function(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"),
                       dat = AC.IPD) {
  # outcome logistic regression fitted to IPD using MCMC (Stan)
  outcome.model <- stan_glm(formula,
                            data = dat,
                            family = binomial,
                            algorithm = "sampling",
                            iter = 4000, warmup = 2000, chains = 2)
  # counterfactual datasets
  data.trtA <- data.trtC <- x_star
  
  # intervene on treatment while keeping set covariates fixed
  data.trtA$trt <- 1 # everyone receives treatment A
  data.trtC$trt <- 0 # all observations receive C
  
  # draw binary responses from posterior predictive distribution
  list(
    y.star.A = posterior_predict(outcome.model, newdata=data.trtA),
    y.star.C = posterior_predict(outcome.model, newdata=data.trtC))
}


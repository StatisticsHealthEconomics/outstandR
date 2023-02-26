
#'
gcomp_stan <- function(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"),
                       dat = AC.IPD) {
  cov_names <- get_covariate_names(formula)
  treat_names <- get_treatment_names(formula)
  
  rho <- cor(AC.IPD[, cov_names])
  
  #  covariate simulation for BC trial using copula package
  cop <-
    normalCopula(param = c(rho[1,2],rho[1,3],rho[1,4],
                                    rho[2,3],rho[2,4],
                                             rho[3,4]),
                 dim = 4,
                 dispstr = "un") # AC IPD pairwise correlations
  
  # sample covariates from approximate joint distribution using copula
  mvd <- mvdc(copula = cop,
              margins = c("norm", "norm", # Gaussian marginals
                          "norm", "norm"),
              # BC covariate means and standard deviations
              paramMargins = list(list(mean=BC.ALD$mean.X1, sd=BC.ALD$sd.X1),
                                  list(mean=BC.ALD$mean.X2, sd=BC.ALD$sd.X2),
                                  list(mean=BC.ALD$mean.X3, sd=BC.ALD$sd.X3),
                                  list(mean=BC.ALD$mean.X4, sd=BC.ALD$sd.X4)))
  # simulated BC pseudo-population of size 1000
  x_star <- as.data.frame(rMvdc(1000, mvd))
  colnames(x_star) <- cov_names
  
  # outcome logistic regression fitted to IPD using MCMC (Stan)
  outcome.model <- stan_glm(formula,
                            data = dat,
                            family = binomial,
                            algorithm = "sampling",
                            iter = 4000, warmup = 2000, chains = 2)
  # counterfactual datasets
  data.trtA <- data.trtC <- x_star
  
  # intervene on treatment while keeping set covariates fixed
  data.trtA[[treat_names]] <- 1 # everyone receives treatment A
  data.trtC[[treat_names]] <- 0 # all observations receive C
  
  # draw binary responses from posterior predictive distribution
  list(
    y.star.A = posterior_predict(outcome.model, newdata=data.trtA),
    y.star.C = posterior_predict(outcome.model, newdata=data.trtC))
}


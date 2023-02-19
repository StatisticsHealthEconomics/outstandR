
#' G-computation maximum likelihood bootstrap
#' 
#' @return mean difference in expected log-odds
#' 
gcomp_ml.boot <- function(data, indices,
                          formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")) {
  dat <- data[indices, ]
  gcomp_ml_log_odds_ratio(formula, dat) 
}


#' marginal A vs. C log-odds ratio (mean difference in expected log-odds)
#' estimated by transforming from probability to linear predictor scale
#
gcomp_ml_log_odds_ratio <- function(formula, dat) {
  rho <- cor(AC.IPD[, c("X1","X2","X3","X4")])
  
  # covariate simulation for BC trial using copula package
  # AC IPD pairwise correlations
  cop <-
    normalCopula(param = c(rho[1,2], rho[1,3], rho[1,4],
                                     rho[2,3], rho[2,4],
                                               rho[3,4]),
                 dim = 4,
                 dispstr = "un")
  
  # sample covariates from approximate joint distribution using copula
  mvd <- mvdc(copula = cop,
              margins = c("norm", "norm",  # Gaussian marginals
                          "norm", "norm"),
              # BC covariate means and standard deviations
              paramMargins = list(list(mean=BC.ALD$mean.X1, sd=BC.ALD$sd.X1),
                                  list(mean=BC.ALD$mean.X2, sd=BC.ALD$sd.X2),
                                  list(mean=BC.ALD$mean.X3, sd=BC.ALD$sd.X3),
                                  list(mean=BC.ALD$mean.X4, sd=BC.ALD$sd.X4)))
  # simulated BC pseudo-population of size 1000
  x_star <- as.data.frame(rMvdc(1000, mvd))
  colnames(x_star) <- c("X1", "X2", "X3", "X4")
  
  # outcome logistic regression fitted to IPD using maximum likelihood
  fit <- glm(formula, data = dat,
             family = binomial)
  
  # counterfactual datasets
  data.trtA <- data.trtC <- x_star
  
  # intervene on treatment while keeping set covariates fixed
  data.trtA$trt <- 1 # dataset where everyone receives treatment A
  data.trtC$trt <- 0 # dataset where all observations receive C
  
  # predict counterfactual event probs, conditional on treatment/covariates
  hat.mu.A.i <- predict(fit, type="response", newdata=data.trtA)
  hat.mu.C.i <- predict(fit, type="response", newdata=data.trtC)
  
  hat.mu.A <- mean(hat.mu.A.i) # (marginal) mean probability prediction under A
  hat.mu.C <- mean(hat.mu.C.i) # (marginal) mean probability prediction under C
  
  log(hat.mu.A/(1-hat.mu.A)) - log(hat.mu.C/(1-hat.mu.C))
  # qlogis(hat.mu.A) - qlogis(hat.mu.C)#'
}

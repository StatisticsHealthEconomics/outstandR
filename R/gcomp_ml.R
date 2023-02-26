
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
  
  covariate_names <- get_covariate_names(formula)
  mean_names <- get_mean_names(formula, dat)
  sd_names <- get_sd_names(formula, dat)
  treat_name <- get_treatment_name(formula)
  
  n_covariates <- length(covariate_names)
  rho <- cor(dat[, covariate_names])
  
  # covariate simulation for BC trial using copula package
  # AC IPD pairwise correlations
  t_rho <- t(rho)  # so extract along rows
  cop <-
    normalCopula(param = t_rho(lower.tri(t_rho)),
                 dim = n_covariates,
                 dispstr = "un")
  
  # BC covariate means & standard deviations
  mean_sd_margins <- list()
  for (i in seq_len(n_covariates)) {
    mean_sd_margins <- c(mean_sd_margins,
                         list(mean = dat[[mean_names[i]]],
                              sd = dat[[sd_names[i]]]))
  }
  
  # sample covariates from approximate joint distribution using copula
  mvd <- mvdc(copula = cop,
              margins = rep("norm", n_covariates),  # Gaussian marginals
              paramMargins = mean_sd_margins)
  
  # simulated BC pseudo-population
  x_star <- as.data.frame(rMvdc(n = 1000, mvd))
  colnames(x_star) <- covariate_names
  
  # outcome logistic regression fitted to IPD using maximum likelihood
  fit <- glm(formula,
             data = dat,
             family = binomial)
  
  # counterfactual datasets
  data.trtA <- data.trtC <- x_star
  
  # intervene on treatment while keeping set covariates fixed
  data.trtA[[treat_name]] <- 1 # dataset where everyone receives treatment A
  data.trtC[[treat_name]] <- 0 # dataset where all observations receive C
  
  # predict counterfactual event probs, conditional on treatment/covariates
  hat.mu.A.i <- predict(fit, type="response", newdata=data.trtA)
  hat.mu.C.i <- predict(fit, type="response", newdata=data.trtC)
  
  hat.mu.A <- mean(hat.mu.A.i) # (marginal) mean probability prediction under A
  hat.mu.C <- mean(hat.mu.C.i) # (marginal) mean probability prediction under C
  
  log(hat.mu.A/(1-hat.mu.A)) - log(hat.mu.C/(1-hat.mu.C))
  # qlogis(hat.mu.A) - qlogis(hat.mu.C)#'
}

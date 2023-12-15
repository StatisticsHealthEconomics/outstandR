
#' G-computation using Stan
#'
#' Calculate draws of binary responses from posterior predictive distribution
#' from the Bayesian G-computation method and Hamiltonian Monte-Carlo.
#' 
#' @param formula Linear regression formula object; default \eqn{y = X_3 + X_4 + \beta_t X_1 + \beta_t X_2} 
#' @template args-ipd
#' @template args-ald
#'
#' @return A list of \eqn{y^*_A} and \eqn{y^*_C} posterior predictions
#' @importFrom copula normalCopula mvdc rMvdc
#' @importFrom rstanarm stan_glm posterior_predict
#' @export
#'
gcomp_stan <- function(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"),
                       ipd, ald) {
  
  if (class(formula) != "formula")
    stop("formula argument must be of formula class.")
  
  treat_names <- get_treatment_name(formula)
  cov_names <- get_covariate_names(formula)
  
  # remove treatment
  cov_names <- cov_names[cov_names != treat_names]
  n_covariates <- length(cov_names)
  rho <- cor(ipd[, cov_names])
  
  # covariate simulation for BC trial using copula package
  cop <-
    copula::normalCopula(param = c(rho[1,2],rho[1,3],rho[1,4],
                                   rho[2,3],rho[2,4],
                                   rho[3,4]),
                         dim = n_covariates,
                         dispstr = "un")  # AC IPD pairwise correlations
  
  # sample covariates from approximate joint distribution using copula
  mvd <- mvdc(copula = cop,
              margins = c("norm", "norm",  # Gaussian marginals
                          "norm", "norm"),
              # BC covariate means and standard deviations
              paramMargins = list(list(mean=ald$mean.X1, sd=ald$sd.X1),
                                  list(mean=ald$mean.X2, sd=ald$sd.X2),
                                  list(mean=ald$mean.X3, sd=ald$sd.X3),
                                  list(mean=ald$mean.X4, sd=ald$sd.X4)))
  # simulated BC pseudo-population of size 1000
  x_star <- as.data.frame(rMvdc(1000, mvd))
  
  colnames(x_star) <- cov_names
  
  # outcome logistic regression fitted to IPD using MCMC (Stan)
  outcome.model <- rstanarm::stan_glm(formula,
                                      data = ipd,
                                      family = binomial,
                                      algorithm = "sampling",
                                      iter = 4000, warmup = 2000, chains = 2)
  # counterfactual datasets
  data.trtA <- data.trtC <- x_star
  
  # intervene on treatment while keeping set covariates fixed
  data.trtA[[treat_names]] <- 1  # everyone receives treatment A
  data.trtC[[treat_names]] <- 0  # all observations receive treatment C
  
  # draw binary responses from posterior predictive distribution
  list(
    y.star.A = rstanarm::posterior_predict(outcome.model, newdata=data.trtA),
    y.star.C = rstanarm::posterior_predict(outcome.model, newdata=data.trtC))
}


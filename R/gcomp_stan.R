
#' G-computation using Stan
#'
#' Calculate draws of binary responses from posterior predictive distribution
#' from the Bayesian G-computation method using Hamiltonian Monte Carlo.
#' 
#' @param formula Linear regression `formula` object
#' @template args-ipd
#' @template args-ald
#'
#' @return A list of \eqn{y^*_A} and \eqn{y^*_C} posterior predictions
#' @importFrom copula normalCopula mvdc rMvdc
#' @importFrom rstanarm stan_glm posterior_predict
#' @keywords internal
#'
gcomp_stan <- function(formula = NULL,
                       ipd, ald) {
  
  if (!inherits(formula, "formula"))
    stop("formula argument must be of formula class.")
  
  x_star <- simulate_ALD_pseudo_pop(formula, ipd, ald)
  
  # outcome logistic regression fitted to IPD using MCMC (Stan)
  outcome.model <-
    rstanarm::stan_glm(formula,
                       data = ipd,
                       family = binomial,
                       algorithm = "sampling",
                       iter = 4000, warmup = 2000, chains = 2)
  
  # counterfactual datasets
  data.trtA <- data.trtC <- x_star
  
  treat_name <- get_treatment_name(formula)
  
  # intervene on treatment while keeping set covariates fixed
  data.trtA[[treat_name]] <- 1  # everyone receives treatment A
  data.trtC[[treat_name]] <- 0  # all observations receive treatment C
  
  # draw binary responses from posterior predictive distribution
  list(
    y.star.A = rstanarm::posterior_predict(outcome.model, newdata=data.trtA),
    y.star.C = rstanarm::posterior_predict(outcome.model, newdata=data.trtC))
}



#' G-computation using Stan
#'
#' Calculate draws of binary responses from posterior predictive distribution
#' from the Bayesian G-computation method using Hamiltonian Monte Carlo.
#' 
#' @param formula Linear regression `formula` object
#' @param family A `family` object
#' @template args-ipd
#' @template args-ald
#'
#' @return A list of \eqn{y^*_A} and \eqn{y^*_C} posterior predictions
#' @importFrom copula normalCopula mvdc rMvdc
#' @importFrom rstanarm stan_glm posterior_predict
#' @export
#'
calc_gcomp_stan <- function(strategy,
                            ipd, ald) {
  
  formula <- strategy$formula
  family <- strategy$family
  iter <- strategy$iter
  warmup <- strategy$warmup
  chain <- strategy$chians
  
  if (!inherits(formula, "formula"))
    stop("formula argument must be of formula class.")
  
  x_star <- simulate_ALD_pseudo_pop(formula, ipd, ald)
  
  # outcome logistic regression fitted to IPD using MCMC (Stan)
  outcome.model <-
    rstanarm::stan_glm(formula,
                       data = ipd,
                       family = family,
                       algorithm = "sampling",
                       iter = iter,
                       warmup = warmup,
                       chains = chains)
  
  # counterfactual datasets
  data.trtA <- data.trtC <- x_star
  
  treat_name <- get_treatment_name(formula)
  
  # intervene on treatment while keeping set covariates fixed
  data.trtA[[treat_name]] <- 1  # everyone receives treatment A
  data.trtC[[treat_name]] <- 0  # receive treatment C
  
  ##TODO: is this going to work for all of the different data types?
  # draw responses from posterior predictive distribution
  y.star.A <- rstanarm::posterior_predict(outcome.model, newdata = data.trtA)
  y.star.C <- rstanarm::posterior_predict(outcome.model, newdata = data.trtC)
  
  # posterior means for each treatment group
  list(
    `0` = rowMeans(y.star.A),
    `1` = rowMeans(y.star.C))
}


#' @export
#' @importFrom boot boot
#' 
calc_gcomp_ml <- function(strategy,
                          ipd, ald) {
  args_list <- 
    list(R = strategy$R,
         formula = strategy$formula,
         family = strategy$family,
         N = strategy$N,
         data = ipd,
         ald = ald)
  
  gcomp_boot <- do.call(boot::boot, c(statistic = gcomp_ml.boot, args_list))
  
  list(mean_A = gcomp_boot$t[, 2],
       mean_C = gcomp_boot$t[, 1])  
}

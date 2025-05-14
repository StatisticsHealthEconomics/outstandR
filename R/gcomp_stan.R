
#' G-computation using Stan
#'
#' Calculate draws of binary responses from posterior predictive distribution
#' from the Bayesian G-computation method using Hamiltonian Monte Carlo.
#' 
#' @param strategy A list specifying the model strategy, including:
#'   \describe{
#'     \item{formula}{A linear regression `formula` object.}
#'     \item{family}{A `family` object specifying the distribution and link function (e.g., `binomial`).}
#'     \item{iter}{Number of iterations for the MCMC sampling.}
#'     \item{warmup}{Number of warmup iterations for the MCMC sampling.}
#'     \item{chains}{Number of MCMC chains.}
#'   }
#' @eval study_data_args(include_ipd = TRUE, include_ald = TRUE)
#'
#' @return A list of \eqn{y^*_A} and \eqn{y^*_C} posterior predictions:
#' \describe{
#'   \item{\code{`0`}}{Posterior means for treatment group C.}
#'   \item{\code{`1`}}{Posterior means for treatment group A.}
#' }
#' @importFrom copula normalCopula mvdc rMvdc
#' @importFrom rstanarm stan_glm posterior_predict
#' @examples
#' \dontrun{
#' strategy <- list(
#'   formula = outcome ~ treatment + age,
#'   family = binomial(),
#'   iter = 2000,
#'   warmup = 500,
#'   chains = 4
#' )
#' ipd <- data.frame(treatment = c(0, 1), outcome = c(1, 0), age = c(30, 40))
#' ald <- data.frame()
#' calc_gcomp_stan(strategy, ipd, ald)
#' }
#' @export
#'
calc_gcomp_stan <- function(strategy,
                            ipd, ald, 
                            ref_trt = "C",
                            comp_trt = "A", ...) {
  
  formula <- strategy$formula
  family <- strategy$family
  rho <- strategy$rho
  N <- strategy$N
  trt_var <- strategy$trt_var
  
  x_star <- simulate_ALD_pseudo_pop(formula, ipd, ald, trt_var, rho, N)
  
  # outcome logistic regression fitted to IPD using MCMC (Stan)
  outcome.model <-
    rstanarm::stan_glm(formula,
                       data = ipd,
                       family = family,
                       algorithm = "sampling",
                       ...)
  
  # counterfactual datasets
  data.comp <- data.ref <- x_star
  
  # intervene on treatment while keeping set covariates fixed
  data.comp[[trt_var]] <- comp_trt  # all receive treatment A
  data.ref[[trt_var]] <- ref_trt    # all receive treatment C
  
  ##TODO: is this going to work for all of the different data types?
  # draw responses from posterior predictive distribution
  y.star.comp <- rstanarm::posterior_predict(outcome.model, newdata = data.comp)
  y.star.ref <- rstanarm::posterior_predict(outcome.model, newdata = data.ref)
  
  # posterior means for each treatment group
  list(
    mean_A = rowMeans(y.star.comp),
    mean_C = rowMeans(y.star.ref))
}


#' G-computation Maximum Likelihood Bootstrap
#'
#' Computes the mean difference in treatment effects using bootstrap resampling.
#'
#' @param strategy A list specifying the model strategy, including:
#'   \describe{
#'     \item{R}{Number of bootstrap replications.}
#'     \item{formula}{A linear regression `formula` object.}
#'     \item{family}{A `family` object specifying the distribution and link function (e.g., `binomial`).}
#'     \item{N}{Synthetic sample size for g-computation.}
#'   }
#' @param ipd Individual patient data.
#' @param ald Aggregate-level data.
#'
#' @return A list containing:
#' \describe{
#'   \item{mean_A}{Bootstrap estimates for treatment group A.}
#'   \item{mean_C}{Bootstrap estimates for treatment group C.}
#' }
#' @importFrom boot boot
#' @examples
#' \dontrun{
#' strategy <- list(
#'   R = 1000,
#'   formula = outcome ~ treatment + age,
#'   family = binomial(),
#'   trt_var = "treatment",
#'   N = 1000
#' )
#' ipd <- data.frame(treatment = c(0, 1), outcome = c(1, 0), age = c(30, 40))
#' ald <- data.frame()
#' calc_gcomp_ml(strategy, ipd, ald)
#' }
#' @export
#'
calc_gcomp_ml <- function(strategy,
                          ipd, ald) {
  args_list <- 
    list(R = strategy$R,
         formula = strategy$formula,
         family = strategy$family,
         trt_var = strategy$trt_var,
         ref_trt = "C",
         comp_trt = "A",
         rho = strategy$rho,
         N = strategy$N,
         data = ipd,
         ald = ald)
  
  gcomp_boot <- do.call(boot::boot, c(statistic = gcomp_ml.boot, args_list))
  
  list(mean_A = gcomp_boot$t[, 2],
       mean_C = gcomp_boot$t[, 1])  
}

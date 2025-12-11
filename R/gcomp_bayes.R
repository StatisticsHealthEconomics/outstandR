
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
#' @param analysis_params List of analysis parameters 
#' @param ... Additional arguments
#'
#' @return A list of \eqn{y^*_A} and \eqn{y^*_C} posterior predictions:
#' \describe{
#'   \item{\code{`0`}}{Posterior means for reference treatment group "C".}
#'   \item{\code{`1`}}{Posterior means for comparator treatment group "A".}
#' }
#' @importFrom copula normalCopula mvdc rMvdc
#' @importFrom rstanarm stan_glm posterior_predict
#' @examples
#' strategy <- list(
#'   formula = y ~ trt:X1,
#'   family = binomial(),
#'   rho = NA,
#'   N = 1000L,
#'   trt_var = "trt",
#'   iter = 2000,
#'   warmup = 500,
#'   chains = 4)
#' 
#' ipd <- data.frame(trt = sample(c("A", "C"), size = 100, replace = TRUE),
#'                   X1 = rnorm(100, 1, 1),
#'                   y = sample(c(1,0), size = 100, prob = c(0.7,0.3), replace = TRUE))
#' 
#' ald <- data.frame(trt = c(NA, NA, "B", "C", "B", "C"),
#'                   variable = c("X1", "X1", "y", "y", NA, NA),
#'                   statistic = c("mean", "sd", "sum", "sum", "N", "N"),
#'                   value = c(0.5, 0.1, 10, 12, 20, 25))
#' 
#' calc_gcomp_bayes(
#'   strategy,
#'   analysis_params = 
#'     list(ipd = ipd, ald = ald, 
#'          ref_trt = "C", ipd_comp = "A"))
#' 
#' @export
#'
calc_gcomp_bayes <- function(strategy,
                             analysis_params, ...) {
  
  formula <- strategy$formula
  family <- strategy$family
  rho <- strategy$rho
  N <- strategy$N
  trt_var <- strategy$trt_var
  ipd <- analysis_params$ipd 
  ald <- analysis_params$ald 
  ref_trt <- analysis_params$ref_trt
  comp_trt <- analysis_params$ipd_comp
  marginal_distns <- strategy$marginal_distns
  marginal_params <- strategy$marginal_params
  
  x_star <- simulate_ALD_pseudo_pop(formula, ipd, ald, trt_var, rho, N,
                                    marginal_distns = marginal_distns,
                                    marginal_params = marginal_params)
  
  # outcome logistic regression fitted to IPD using MCMC (Stan)
  outcome.model <-
    rstanarm::stan_glm(formula,
                       data = ipd,
                       family = family,
                       algorithm = "sampling",
                       ...)
  
  # counterfactual datasets
  data_comp <- data_ref <- x_star
  
  # intervene on treatment while keeping set covariates fixed
  data_comp[[trt_var]] <- comp_trt  # all receive comparator treatment
  data_ref[[trt_var]] <- ref_trt    # all receive reference treatment
  
  ##TODO: is this going to work for all of the different data types?
  
  # draw responses from posterior predictive distribution
  y.star.comp <- rstanarm::posterior_predict(outcome.model, newdata = data_comp)
  y.star.ref  <- rstanarm::posterior_predict(outcome.model, newdata = data_ref)
  
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
#' @param analysis_params List of analysis parameters
#' @return A list containing:
#' \describe{
#'   \item{mean_A}{Bootstrap estimates for comparator treatment group "A".}
#'   \item{mean_C}{Bootstrap estimates for reference treatment group "C".}
#' }
#' @importFrom boot boot
#' @examples
#' strategy <- list(
#'   formula = y ~ trt:X1,
#'   family = binomial(),
#'   rho = NA,
#'   N = 1000L,
#'   R = 1000L,
#'   marginal_distns = NA,
#'   marginal_params = NA,
#'   trt_var = "trt")
#' 
#' ipd <- data.frame(trt = sample(c("A", "C"), size = 100, replace = TRUE),
#'                   X1 = rnorm(100, 1, 1),
#'                   y = sample(c(1,0), size = 100, prob = c(0.7,0.3), replace = TRUE))
#' 
#' ald <- data.frame(trt = c(NA, NA, "B", "C", "B", "C"),
#'                   variable = c("X1", "X1", "y", "y", NA, NA),
#'                   statistic = c("mean", "sd", "sum", "sum", "N", "N"),
#'                   value = c(0.5, 0.1, 10, 12, 20, 25))
#' 
#' calc_gcomp_ml(
#'   strategy,
#'   analysis_params = 
#'     list(ipd = ipd, ald = ald, 
#'          ref_trt = "C", ipd_comp = "A"))
#'          
#' @export
#'
calc_gcomp_ml <- function(strategy,
                          analysis_params) {
  
  args_list <- 
    list(R = strategy$R,
         formula = strategy$formula,
         family = strategy$family,
         trt_var = strategy$trt_var,
         ref_trt = analysis_params$ref_trt,
         comp_trt = analysis_params$ipd_comp,
         rho = strategy$rho,
         N = strategy$N,
         marginal_distns = strategy$marginal_distns,
         marginal_params = strategy$marginal_params,
         data = analysis_params$ipd,
         ald = analysis_params$ald)
  
  gcomp_boot <- do.call(boot::boot, c(statistic = gcomp_ml.boot, args_list))
  
  list(mean_A = gcomp_boot$t[, 2],
       mean_C = gcomp_boot$t[, 1])  
}

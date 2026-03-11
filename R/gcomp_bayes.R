
#' Bayesian G-computation using Stan
#'
#' Calculate draws of binary responses from posterior predictive distribution
#' from the Bayesian G-computation method using Hamiltonian Monte Carlo.
#'
#' @param strategy A list specifying the model strategy, including:
#'   * `outcome_model`: A linear regression `formula` object.
#'   * `family`: A `family` object specifying the distribution and link function
#'     (e.g., `binomial`).
#'   * `iter`: Number of iterations for the MCMC sampling.
#'   * `warmup`: Number of warmup iterations for the MCMC sampling.
#'   * `chains`: Number of MCMC chains.
#' @param analysis_params List of analysis parameters. Must contain `ipd` and `ald`.
#' @param ... Additional arguments passed to [rstanarm::stan_glm()].
#'
#' @return A list containing:
#' * `means`: A list containing:
#'     * Posterior means for comparator treatment group.
#'     * Posterior means for reference treatment group.
#' * `model`: A list containing the `fit` object (from `stan_glm`), `rho`, `N`,
#'   and `stan_args`.
#'
#' @importFrom copula normalCopula mvdc rMvdc
#' @importFrom rstanarm stan_glm posterior_predict
#' @examples
#' strategy <- list(
#'   outcome_model = y ~ trt:X1,
#'   family = binomial(),
#'   rho = NA,
#'   N = 1000L,
#'   marginal_distns = NA,
#'   marginal_params = NA,
#'   trt_var = "trt",
#'   iter = 2000,
#'   warmup = 500,
#'   chains = 4)
#' 
#' ipd <- data.frame(
#'    trt = sample(c("A", "C"), size = 100, replace = TRUE),
#'    X1 = rnorm(100, 1, 1),
#'    y = sample(c(1,0), size = 100, prob = c(0.7, 0.3), replace = TRUE))
#' 
#' ald <- data.frame(
#'   trt = c(NA, NA, "B", "C", "B", "C"),
#'   variable = c("X1", "X1", "y", "y", NA, NA),
#'   statistic = c("mean", "sd", "sum", "sum", "N", "N"),
#'   value = c(0.5, 0.1, 10, 12, 20, 25))
#' 
#' calc_gcomp_bayes(
#'   strategy,
#'   analysis_params = list(
#'     ipd = ipd, ald = ald, 
#'     ref_trt = "C",
#'     ipd_comp = "A"))
#' 
#' @export
#'
calc_gcomp_bayes <- function(strategy,
                             analysis_params, ...) {
  
  verbose <- isTRUE(analysis_params$verbose)
  
  # extract seed, or default to a random integer if missing/NULL
  bayes_seed <- analysis_params$seed
  
  if (is.null(bayes_seed)) {
    bayes_seed <- sample.int(.Machine$integer.max, 1)
  }
  
  # If verbose is TRUE, Stan print every 500 iters. FALSE, silence it (0).
  stan_refresh <- if (verbose) 500 else 0
  
  default_stan_args <- list(
    algorithm = "sampling",
    chains = 2,
    iter = 2000,
    refresh = stan_refresh,
    seed = bayes_seed
  )
  
  if (verbose) {
    cli::cli_h2("G-Computation (Bayesian) Execution")
    cli::cli_alert_info("Compiling/Sampling Stan model...")
  }
  
  # merge with user-provided dots
  user_args <- list(...)
  stan_args <- modifyList(default_stan_args, user_args)
  
  formula <- strategy$outcome_model
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
  
  # outcome logistic regression fitted to IPD using MCMC (Stan)
  outcome_model <- do.call(rstanarm::stan_glm, c(
    list(formula = formula, data = ipd, family = family),
    stan_args
  ))
  
  # simulate ALD pseudo-population ONCE (N handles Monte Carlo Integration)
  x_star <- simulate_ALD_pseudo_pop(formula, ipd, ald, trt_var, rho, N,
                                    marginal_distns = marginal_distns,
                                    marginal_params = marginal_params)
  # counterfactual datasets
  data_comp <- data_ref <- x_star
  
  # Intervene safely, maintaining factor levels
  if (is.factor(ipd[[trt_var]])) {
    data_comp[[trt_var]] <- factor(comp_trt, levels = levels(ipd[[trt_var]]))
    data_ref[[trt_var]]  <- factor(ref_trt, levels = levels(ipd[[trt_var]]))
  } else {
    data_comp[[trt_var]] <- comp_trt  
    data_ref[[trt_var]]  <- ref_trt   
  }
  
  # # draw responses from posterior predictive distribution
  # y.star.comp <- rstanarm::posterior_predict(outcome_model, newdata = data_comp)
  # y.star.ref  <- rstanarm::posterior_predict(outcome_model, newdata = data_ref)
  
  # Draw EXPECTED responses (posterior_epred) over the fixed x_star grid
  y.star.comp <- rstanarm::posterior_epred(outcome_model, newdata = data_comp)
  y.star.ref  <- rstanarm::posterior_epred(outcome_model, newdata = data_ref)
  
  comp_draws <- rowMeans(y.star.comp)
  ref_draws  <- rowMeans(y.star.ref)
  
  # harmonized return structure with dynamic naming
  means_list <- 
    stats::setNames(list(comp_draws, ref_draws), 
                    c(comp_trt, ref_trt))
  point_est_list <- 
    stats::setNames(list(mean(comp_draws), mean(ref_draws)), 
                    c(comp_trt, ref_trt))
  
  list(
    means = means_list,
    point_estimates = point_est_list,
    model = list(
      fit = outcome_model,
      rho = rho,
      N = N,
      stan_args = stan_args)
  )
}

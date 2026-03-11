
#' G-computation Maximum Likelihood Bootstrap
#'
#' Computes the mean difference in treatment effects using bootstrap resampling.
#'
#' @param strategy A list specifying the model strategy, including:
#'   * `R`: Number of bootstrap replications.
#'   * `outcome_model`: A linear regression `formula` object for the outcome model.
#'   * `family`: A `family` object specifying the distribution and link function
#'     (e.g., `binomial`).
#'   * `N`: Synthetic sample size for g-computation.
#' @param analysis_params List of analysis parameters.
#'
#' @return A list containing:
#' * `means`: A list containing:
#'     * `A`: Bootstrap estimates for comparator treatment group "A".
#'     * `C`: Bootstrap estimates for reference treatment group "C".
#' * `model`: A list containing the `fit` object, `rho`, and `N`.
#'
#' @importFrom boot boot
#' @examples
#' strategy <- list(
#'   outcome_model = y ~ trt:X1,
#'   family = binomial(),
#'   rho = NA,
#'   N = 1000L,
#'   n_boot = 100L,
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
#'          ref_trt = "C", 
#'          ipd_comp = "A"))
#'          
#' @export
#'
calc_gcomp_ml <- function(strategy,
                          analysis_params) {
  
  verbose <- isTRUE(analysis_params$verbose)
  
  if (verbose) {
    cli::cli_h2("G-Computation (ML) Execution")
    cli::cli_alert_info("Fitting initial model...")
  }
  
  # Extract treatment names for dynamic naming
  ref_trt <- analysis_params$ref_trt
  comp_trt <- analysis_params$ipd_comp
  
  common_args <- list(
    outcome_model = strategy$outcome_model,
    family = strategy$family,
    trt_var = strategy$trt_var,
    ref_trt = ref_trt,
    comp_trt = comp_trt,
    rho = strategy$rho,
    N = strategy$N,
    marginal_distns = strategy$marginal_distns,
    marginal_params = strategy$marginal_params,
    ald = analysis_params$ald)
  
  # single run for fit
  args_orig <- c(common_args, list(ipd = analysis_params$ipd))
  original_run <- do.call(gcomp_ml_means, args_orig)
  
  args_boot <- c(
    common_args, list(
      data = analysis_params$ipd, 
      R = strategy$n_boot))
  
  if (verbose) {
    cli::cli_alert_info("Starting Bootstrap: {.val {strategy$n_boot}} replicates.")
    cli::cli_alert_info("Simulating pseudo-pop (N={strategy$N}) per replicate.")
    
    total_ops <- strategy$n_boot * strategy$N
    if (total_ops > 5e5) {
      cli::cli_alert_warning("Large computation detected ({.val {total_ops}} simulations). Grab a coffee.")
    }
  }
  
  # Run Bootstrap
  gcomp_boot <- do.call(boot::boot, c(statistic = gcomp_ml.boot, args_boot))
  
  # Dynamically assign names to the bootstrap distributions
  # (Assuming boot returns: col 1 = reference, col 2 = comparator based on original logic)
  means_list <- stats::setNames(
    list(gcomp_boot$t[, 2], gcomp_boot$t[, 1]), 
    c(comp_trt, ref_trt)
  )
  
  # Extract original sample point estimates from boot$t0 
  point_est_list <- stats::setNames(
    list(gcomp_boot$t0[2], gcomp_boot$t0[1]), 
    c(comp_trt, ref_trt)
  )
  
  list(
    means = means_list,
    point_estimates = point_est_list,
    model = list(
      fit = original_run$model,
      rho = strategy$rho,
      N = strategy$N,
      n_boot = strategy$n_boot)
  )
}

#' Bootstrap for G-computation via Maximum Likelihood
#'
#' For use with a bootstrapping function (e.g., [boot::boot()]). 
#' On each bootstrap sample of the data, it calculates
#' a relative treatment effect (e.g., log odds ratio, log relative risk, or
#' risk difference) using G-computation with maximum likelihood.
#'     
#' @param data A data frame containing the original individual participant data (IPD).
#' @param indices A vector of indices supplied by the bootstrapping function,
#'   used to resample `data`.
#' @eval reg_args(include_formula = TRUE, include_family = TRUE)
#' @param rho A named square matrix specifying the correlation between covariates
#'   for synthetic data generation. Defaults to `NA`, assuming independence.
#' @param N Synthetic sample size for G-computation
#' @param marginal_distns,marginal_params Marginal distributions and parameters
#' @param ald A data frame of aggregate-level data providing covariate distributions.
#'
#' @return A single numeric value representing the relative treatment effect
#' 
#' @seealso [strategy_gcomp_ml()]
#' 
#' @keywords internal
#' 
gcomp_ml.boot <- function(data, indices,
                          R, 
                          outcome_model = NULL,
                          family, trt_var,
                          ref_trt = NA,
                          comp_trt = NA,
                          rho = NA,
                          N = 1000, 
                          marginal_distns = NA,
                          marginal_params = NA,
                          ald) {
  dat <- data[indices, ]
  
  res <- gcomp_ml_means(outcome_model, family, dat, ald, trt_var, rho, N,
                        ref_trt, comp_trt,
                        marginal_distns, marginal_params)
  return(res$stats)
}


#' G-computation maximum likelihood mean outcomes by arm
#'
#' @eval reg_args(include_formula = TRUE, include_family = TRUE)
#' @eval study_data_args(include_ipd = TRUE, include_ald = TRUE)
#' @param rho A named square matrix of covariate correlations; default NA.
#' @param N Synthetic sample size for g-computation
#' @param marginal_distns,marginal_params Marginal distributions and parameters
#'
#' @return A list containing:
#'   * `stats`: Named vector of marginal mean probabilities
#'   * `model`: The fitted glm object
#'   
#' @seealso [strategy_gcomp_ml()] [gcomp_ml.boot()]
#' @importFrom copula normalCopula mvdc rMvdc
#' @importFrom stats predict glm
#' @keywords internal
#'
gcomp_ml_means <- function(outcome_model,
                           family,
                           ipd, ald,
                           trt_var,
                           rho = NA,
                           N = 1000,
                           ref_trt,
                           comp_trt,
                           marginal_distns = NA,
                           marginal_params = NA) {
  
  x_star <- simulate_ALD_pseudo_pop(formula = outcome_model,
                                    ipd = ipd, ald = ald,
                                    trt_var = trt_var, 
                                    rho = rho, N = N,
                                    marginal_distns = marginal_distns,
                                    marginal_params = marginal_params)
  
  # outcome logistic regression fitted to IPD using maximum likelihood
  fit <- glm(formula = outcome_model,
             family = family,
             data = ipd)
  
  # counterfactual datasets
  data.comp <- data.ref <- x_star
  
  # intervene on treatment while keeping set covariates fixed
  data.comp[[trt_var]] <- comp_trt  # all receive A
  data.ref[[trt_var]] <- ref_trt    # all receive C
  
  # predict counterfactual event probs, conditional on treatment/covariates
  hat.mu.comp <- predict(fit, type = "response", newdata = data.comp)
  hat.mu.ref <- predict(fit, type = "response", newdata = data.ref)
  
  list(
    # (marginal) mean probability prediction under A and C
    stats = c(`0` = mean(hat.mu.ref),
              `1` = mean(hat.mu.comp)),
    model = fit
  )
}




#' Bootstrap for G-computation via Maximum Likelihood
#'
#' This is a statistic function intended for use with a bootstrapping function
#' (e.g., [boot::boot()]). On each bootstrap sample of the data, it calculates
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
                          R, formula = NULL,
                          family, trt_var,
                          ref_trt = NA,
                          comp_trt = NA,
                          rho = NA,
                          N = 1000, 
                          marginal_distns = NA,
                          marginal_params = NA,
                          ald) {
  dat <- data[indices, ]
  
  res <- gcomp_ml_means(formula, family, dat, ald, trt_var, rho, N,
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
#' @seealso [strategy_gcomp_ml()], [gcomp_ml.boot()]
#' @importFrom copula normalCopula mvdc rMvdc
#' @importFrom stats predict glm
#' @keywords internal
#'
gcomp_ml_means <- function(formula,
                           family,
                           ipd, ald,
                           trt_var,
                           rho = NA,
                           N = 1000,
                           ref_trt,
                           comp_trt,
                           marginal_distns = NA,
                           marginal_params = NA) {
  
  x_star <- simulate_ALD_pseudo_pop(formula = formula,
                                    ipd = ipd, ald = ald,
                                    trt_var = trt_var, 
                                    rho = rho, N = N,
                                    marginal_distns = marginal_distns,
                                    marginal_params = marginal_params)
  
  # outcome logistic regression fitted to IPD using maximum likelihood
  fit <- glm(formula = formula,
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


#' Simulate Aggregate-Level Data Pseudo-Population
#'
#' Generates a synthetic cohort using a normal copula based on aggregate-level data.
#'
#' @eval reg_args(include_formula = TRUE, include_family = FALSE)
#' @eval study_data_args(include_ipd = TRUE, include_ald = TRUE)
#' @param rho A named square matrix of covariate correlations or single value; default NA takes from IPD.
#' @param N Sample size for the synthetic cohort. Default is 1000.
#' @param marginal_distns Marginal distributions names; vector default NA.
#'    Available distributions are given in stats::Distributions. See [copula::Mvdc()] for details
#' @param marginal_params Marginal distributions parameters;
#'    named list of lists, default NA. See [copula::Mvdc()] for details
#' @param seed Random seed
#' @param verbose Default `FALSE`
#' 
#' @return A data frame representing the synthetic pseudo-population.
#' @importFrom copula normalCopula mvdc
#' 
#' @keywords internal
#' 
simulate_ALD_pseudo_pop <- function(formula,
                                    ipd = NULL, ald = NULL,
                                    trt_var,
                                    rho = NA,
                                    N = 1000,
                                    marginal_distns = NA,
                                    marginal_params = NA,
                                    seed = NULL,
                                    verbose = FALSE) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  resolved <- 
    prepare_covariate_distns(
      formula, ald, trt_var, 
      marginal_distns, marginal_params, 
      verbose)
  
  marginal_distns <- resolved$distns
  marginal_params <- resolved$params
  
  covariate_names <- names(marginal_distns)
  n_covariates <- length(covariate_names)
  
  # --- Standard Simulation Logic ---
  
  # handle no covariates (intercept-only or treatment-only models)
  if (n_covariates == 0) {
    return(data.frame(row.names = seq_len(length.out = N)))
  }
  
  # don't require copula for single covariate
  if (n_covariates == 1) {
    # dynamically call appropriate random number generator
    rng_fun <- get(paste0("r", marginal_distns[1]))
    sim_vals <- do.call(rng_fun, c(list(n=N), marginal_params[[1]]))
    
    x_star <- matrix(sim_vals, ncol = 1,
                     dimnames = list(NULL, covariate_names))
    # return(as.data.frame(setNames(list(sim_vals), covariate_names)))  ##TODO
    return(as.data.frame(x_star))
  }
  
  # prepare correlation matrix
  if (!is.matrix(rho)) {
    if (is.na(rho)) {
      
      if (is.null(ipd)) {
        stop("'rho' must be provided when 'ipd' is not available.", call. = FALSE)
      }
      
      rho <- cor(ipd[, covariate_names], use = "pairwise.complete.obs")
    } else {
      rho <- matrix(rho, n_covariates, n_covariates,
                    dimnames = list(covariate_names, covariate_names))
      diag(rho) <- 1
    }
  }
  
  rho <- rho[covariate_names, covariate_names]
  
  # define Copula and Multivariate Distribution
  cor_params <- rho[lower.tri(rho, diag = FALSE)]
  
  cop <- copula::normalCopula(param = cor_params,
                              dim = n_covariates,
                              dispstr = "un")
  
  mvd <- copula::mvdc(copula = cop,
                      margins = marginal_distns,
                      paramMargins = marginal_params)
  
  # simulate data
  x_star <- as.data.frame(copula::rMvdc(n = N, mvd))
  colnames(x_star) <- covariate_names
  
  return(x_star)
}

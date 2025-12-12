
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
                                    seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  covariate_names <- get_covariate_names(formula)
  covariate_names <- covariate_names[covariate_names != trt_var]  # remove treatment
  n_covariates <- length(covariate_names)
  
  if (n_covariates == 0) {
    stop("No covariates found to simulate.", call. = FALSE)
  }
  
  # --- Logic for User-Defined Distributions with Auto-Parameterization ---
  
  # CASE 1: No distributions provided -> Full Auto-Detection (Norm/Binom only)
  # default behaviour
  if (all(is.na(marginal_distns))) {
    
    auto_distns <- character(n_covariates)
    auto_params <- vector("list", n_covariates)
    names(auto_params) <- names(auto_distns) <- covariate_names
    
    for (cov in covariate_names) {
      var_ald <- dplyr::filter(ald, .data$variable == cov)
      
      if (nrow(var_ald) == 0) {
        stop(paste("No ALD found for covariate: ", cov), call. = FALSE)
      }
      
      if ("prop" %in% var_ald$statistic) {
        auto_distns[cov] <- "binom"
        prob <- var_ald$value[var_ald$statistic == "prop"]
        auto_params[[cov]] <- list(size = 1, prob = prob)
        
      } else if (all(c("mean", "sd") %in% var_ald$statistic)) {
        auto_distns[cov] <- "norm"
        mean_val <- var_ald$value[var_ald$statistic == "mean"]
        sd_val <- var_ald$value[var_ald$statistic == "sd"]
        auto_params[[cov]] <- list(mean = mean_val, sd = sd_val)
        
      } else {
        stop(paste("For", cov, "provide 'prop' or 'mean'/'sd' in ALD."), call. = FALSE)
      }
    }
    marginal_distns <- auto_distns
    marginal_params <- auto_params
    
  } else {
    # CASE 2: Distributions provided -> Fill missing parameters
    message("Using user-supplied marginal distributions.")
    
    # Initialize params list if missing
    if (!is.list(marginal_params)) {
      marginal_params <- vector("list", n_covariates)
      names(marginal_params) <- covariate_names
    }
    
    # Handle unnamed distribution vectors (assume order matches formula)
    if (is.null(names(marginal_distns))) {
      if (length(marginal_distns) == n_covariates) {
        names(marginal_distns) <- covariate_names
      } else {
        stop("marginal_distns must be named or match number of covariates")
      }
    }
    
    for (cov in covariate_names) {
      # If parameters already exist for this covariate, skip calculation
      if (!is.null(marginal_params[[cov]])) next
      
      dist <- marginal_distns[cov]
      if (is.na(dist)) stop("Missing distribution for ", cov)
      
      var_ald <- dplyr::filter(ald, .data$variable == cov)
      if (nrow(var_ald) == 0) stop("No ALD found for ", cov)
      
      # Extract common stats
      m <- var_ald$value[var_ald$statistic == "mean"]
      s <- var_ald$value[var_ald$statistic == "sd"]
      p <- var_ald$value[var_ald$statistic == "prop"]
      
      # --- Parameter Transformations ---
      if (dist == "norm") {
        if (length(m) == 0 || length(s) == 0) stop("Need mean/sd for norm: ", cov)
        marginal_params[[cov]] <- list(mean = m, sd = s)
        
      } else if (dist == "binom") {
        if (length(p) == 0) stop("Need prop for binom: ", cov)
        marginal_params[[cov]] <- list(size = 1, prob = p)
        
      } else if (dist == "gamma") {
        # Method of Moments: shape = m^2/s^2, rate = m/s^2
        if (length(m) == 0 || length(s) == 0) stop("Need mean/sd for gamma: ", cov)
        rate_val <- m / (s^2)
        shape_val <- m^2 / (s^2)
        marginal_params[[cov]] <- list(shape = shape_val, rate = rate_val)
        
      } else if (dist == "lnorm") {
        # Method of Moments: Log-normal
        if (length(m) == 0 || length(s) == 0) stop("Need mean/sd for lnorm: ", cov)
        var_log <- log(1 + (s^2 / m^2))
        mean_log <- log(m) - 0.5 * var_log
        sd_log <- sqrt(var_log)
        marginal_params[[cov]] <- list(meanlog = mean_log, sdlog = sd_log)
        
      } else if (dist == "beta") {
        # Method of Moments: Beta
        if (length(m) == 0 || length(s) == 0) stop("Need mean/sd for beta: ", cov)
        if (s^2 >= m * (1 - m)) stop("Variance too high for Beta distribution: ", cov)
        
        term <- (m * (1 - m) / s^2) - 1
        shape1 <- m * term
        shape2 <- (1 - m) * term
        marginal_params[[cov]] <- list(shape1 = shape1, shape2 = shape2)
        
      } else {
        # Fallback for other distributions (e.g. t, weibull) not implemented
        stop(paste("Automatic parameter transformation not implemented for '", dist, 
                   "'. Please supply marginal_params for: ", cov, sep=""))
      }
    }
  }
  
  # --- Standard Simulation Logic ---
  
  # don't require copula for single covariate
  if (n_covariates <= 1) {
    # dynamically call appropriate random number generator
    rng_fun <- get(paste0("r", marginal_distns[1]))
    sim_vals <- do.call(rng_fun, c(list(n=N), marginal_params[[1]]))
    
    x_star <- matrix(sim_vals, ncol = 1,
                     dimnames = list(NULL, covariate_names))
    return(as.data.frame(x_star))
  }
  
  # prepare correlation matrix
  if (!is.matrix(rho)) {
    if (is.na(rho)) {
      
      if (is.null(ipd)) {
        stop("'rho' must be provided when 'ipd' is not available.")
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
  
  # marginal_params add size = 1 if 'binom' and size is missing
  for (i in seq_along(marginal_distns)) {
    if (marginal_distns[i] == "binom" &&
        is.null(marginal_params[[i]]$size)) {
      marginal_params[[i]]$size <- 1
    }
  }
  
  mvd <- copula::mvdc(copula = cop,
                      margins = marginal_distns,
                      paramMargins = marginal_params)
  
  # simulate data
  x_star <- as.data.frame(copula::rMvdc(n = N, mvd))
  colnames(x_star) <- covariate_names
  
  return(x_star)
}



#' Create counterfactual datasets for treatment comparison
#'
#' Creates two copies of a dataset with treatment variable set to different values.
#'
#' @param x_star Base dataset to create counterfactuals from
#' @param trt_var Name of the treatment variable
#' @param comp_trt Comparator treatment value
#' @param ref_trt Reference treatment value
#' @return A list with two elements:
#'   \describe{
#'     \item{comp}{Dataset with treatment set to comp_trt}
#'     \item{ref}{Dataset with treatment set to ref_trt}
#'   }
#' @keywords internal
#'
create_counterfactual_datasets <- function(x_star, trt_var, comp_trt, ref_trt) {
  data_comp <- data_ref <- x_star
  
  # intervene on treatment while keeping set covariates fixed
  data_comp[[trt_var]] <- comp_trt
  data_ref[[trt_var]] <- ref_trt
  
  list(comp = data_comp, ref = data_ref)
}


#' Extract treatment coefficient name from model coefficients
#'
#' Safely extracts the treatment coefficient name from model coefficients,
#' handling potential factor level appends.
#'
#' @param coef_names Names of model coefficients
#' @param trt_var Name of the treatment variable
#' @return The name of the treatment coefficient
#' @keywords internal
#'
extract_treatment_coef_name <- function(coef_names, trt_var) {
  grep(pattern = paste0("^", trt_var, "[^:]*$"),
       coef_names, value = TRUE)
}


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
#' @param ald A data frame of aggregate-level data providing covariate distributions.
#'
#' @return A single numeric value representing the relative treatment effect
#' 
#' @seealso [strategy_gcomp_ml()]
#' 
#' @examples
#' \dontrun{
#' data <- data.frame(trt = c("A", "C"),
#'                    y = c(1, 0))
#' gcomp_ml.boot(data, indices = 1:2, formula = y ~ trt,
#'               R = 100, family = binomial(), N = 1000, ald = NULL)
#' }
#' @keywords internal
#' 
gcomp_ml.boot <- function(data, indices,
                          R, formula = NULL,
                          family, trt_var,
                          ref_trt = NA,
                          comp_trt = NA,
                          rho = NA,
                          N = 1000, ald) {
  dat <- data[indices, ]
  gcomp_ml_means(formula, family, dat, ald, trt_var, rho, N,
                 ref_trt, comp_trt) 
}


#' G-computation maximum likelihood mean outcomes by arm
#'
#' @eval reg_args(include_formula = TRUE, include_family = TRUE)
#' @eval study_data_args(include_ipd = TRUE, include_ald = TRUE)
#' @param rho A named square matrix of covariate correlations; default NA.
#' @param N Synthetic sample size for g-computation
#'
#' @return A named vector containing the marginal mean probabilities under
#'   comparator "A" (`0`) and reference "C" (`1`) treatments.
#' @seealso [strategy_gcomp_ml()], [gcomp_ml.boot()]
#' @importFrom copula normalCopula mvdc rMvdc
#' @importFrom stats predict glm
#' @examples
#' \dontrun{
#' formula <- y ~ trt
#' family <- binomial()
#' ipd <- data.frame(trt = c("A", "C"),
#'                    y = c(1, 0))
#' ald <- data.frame()
#' gcomp_ml_means(formula, family, N = 1000, ipd = ipd, ald = ald)
#' }
#' @keywords internal
#'
gcomp_ml_means <- function(formula,
                           family,
                           ipd, ald,
                           trt_var,
                           rho = NA,
                           N = 1000,
                           ref_trt,
                           comp_trt) {
  
  x_star <- simulate_ALD_pseudo_pop(formula = formula,
                                    ipd = ipd, ald = ald,
                                    trt_var = trt_var, 
                                    rho = rho, N = N)
  
  # outcome logistic regression fitted to IPD using maximum likelihood
  fit <- glm(formula = formula,
             family = family,
             data = ipd)
  
  # counterfactual datasets
  counterfactuals <- create_counterfactual_datasets(x_star, trt_var, comp_trt, ref_trt)
  
  # predict counterfactual event probs, conditional on treatment/covariates
  hat.mu.comp <- predict(fit, type = "response", newdata = counterfactuals$comp)
  hat.mu.ref <- predict(fit, type = "response", newdata = counterfactuals$ref)
  
  # (marginal) mean probability prediction under A and C
  c(`0` = mean(hat.mu.ref),
    `1` = mean(hat.mu.comp))
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
#' @examples
#' \dontrun{
#' formula <- y ~ trt + age
#' ipd <- data.frame(tr = c("A", "C"), y = c(1, 0), age = c(30, 40))
#' ald <- data.frame()
#' simulate_ALD_pseudo_pop(formula, ipd, ald, trt_var = "trt", N = 1000)
#' }
#' @keywords internal
#' 
simulate_ALD_pseudo_pop <- function(formula,
                                    ipd = NULL, ald = NULL,
                                    trt_var,
                                    rho = NA,
                                    N = 1000,
                                    marginal_distns = NA,
                                    marginal_params = NA) {
  set.seed(1234)
  
  covariate_names <- get_covariate_names(formula)
  covariate_names <- covariate_names[covariate_names != trt_var]  # remove treatment
  n_covariates <- length(covariate_names)
  
  if (n_covariates == 0) {
    stop("No covariates found to simulate.")
  }
  
  if (is.character(marginal_distns) && is.list(marginal_params)) {
    message("user-supplied marginals.")
  } else {
    auto_distns <- character(n_covariates)
    auto_params <- vector("list", n_covariates)
    names(auto_params) <- names(auto_distns) <- covariate_names
    
    for (cov in covariate_names) {
      var_ald <- dplyr::filter(ald, .data$variable == cov)
      
      if (nrow(var_ald) == 0) {
        stop(paste("No ALD found for covariate: ", cov))
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
        stop(paste("For", cov, "provide 'prop' or 'mean'/'sd' in ALD."))
      }
    }
    
    marginal_distns <- auto_distns
    marginal_params <- auto_params
  }
  
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

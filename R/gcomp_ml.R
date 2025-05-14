
#' G-computation maximum likelihood bootstrap
#' 
#' Using bootstrap resampling, calculates the log odds ratio.
#'     
#' @param data Trial data 
#' @param indices Indices sampled from rows of `data`
#' @eval reg_args(include_formula = TRUE, include_family = TRUE)
#' @param rho A named square matrix of covariate correlations; default NA.
#' @param N Synthetic sample size for g-computation
#' @param ald Aggregate-level data for covariates.
#'
#' @return Mean difference in expected log-odds
#' @seealso [strategy_gcomp_ml()], [gcomp_ml_log_odds_ratio()]
#' @examples
#' \dontrun{
#' data <- data.frame(treatment = c(0, 1), outcome = c(1, 0))
#' gcomp_ml.boot(data, indices = 1:2, formula = outcome ~ treatment,
#'               R = 100, family = binomial(), N = 1000, ald = NULL)
#' }
#' @keywords internal
#' 
gcomp_ml.boot <- function(data, indices,
                          R, formula = NULL,
                          family, trt_var,
                          ref_trt = "C", comp_trt = "A",
                          rho = NA,
                          N = 1000, ald) {
  dat <- data[indices, ]
  gcomp_ml_means(formula, family, dat, ald, trt_var, rho, N,
                 ref_trt, comp_trt) 
}


#' G-computation Maximum Likelihood mean outcome
#' 
#' @section Log-Odds Ratio: 
#' Marginal _A_ vs _C_ log-odds ratio (mean difference in expected log-odds)
#' estimated by transforming from probability to linear predictor scale.
#'
#' \eqn{\log(\hat{\mu}_A/(1 - \hat{\mu}_A)) - \log(\hat{\mu}_C/(1 - \hat{\mu}_C))}
#'
#' @eval reg_args(include_formula = TRUE, include_family = TRUE)
#' @eval study_data_args(include_ipd = TRUE, include_ald = TRUE)
#' @param rho A named square matrix of covariate correlations; default NA.
#' @param N Synthetic sample size for g-computation
#'
#' @return A named vector containing the marginal mean probabilities under treatments A (`0`) and C (`1`).
#' @seealso [strategy_gcomp_ml()], [gcomp_ml.boot()]
#' @importFrom copula normalCopula mvdc rMvdc
#' @importFrom stats predict glm
#' @examples
#' \dontrun{
#' formula <- outcome ~ treatment
#' family <- binomial()
#' ipd <- data.frame(treatment = c(0, 1), outcome = c(1, 0))
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
                           ref_trt = "C",
                           comp_trt = "A") {
  
  x_star <- simulate_ALD_pseudo_pop(formula = formula,
                                    ipd = ipd, ald = ald,
                                    trt_var = trt_var, 
                                    rho = rho, N = N)
  
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
  hat.mu.comp <- predict(fit, type="response", newdata=data.comp)
  hat.mu.ref <- predict(fit, type="response", newdata=data.ref)
  
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
#' @param rho A named square matrix of covariate correlations; default NA.
#' @param N Sample size for the synthetic cohort. Default is 1000.
#' 
#' @return A data frame representing the synthetic pseudo-population.
#' @importFrom copula normalCopula mvdc
#' @examples
#' \dontrun{
#' formula <- outcome ~ treatment + age
#' ipd <- data.frame(treatment = c(0, 1), outcome = c(1, 0), age = c(30, 40))
#' ald <- data.frame()
#' simulate_ALD_pseudo_pop(formula, ipd, ald, trt_var = "treatment", N = 1000)
#' }
#' @keywords internal
#' 
simulate_ALD_pseudo_pop <- function(formula,
                                    ipd, ald,
                                    trt_var,
                                    rho = NA,
                                    N = 1000) {
  set.seed(1234)
  
  covariate_names <- get_covariate_names(formula)
  
  # remove treatment
  covariate_names <- covariate_names[covariate_names != trt_var]
  
  mean_names <- get_mean_names(ald, covariate_names)
  
  sd_names <- get_sd_names(ald, covariate_names)
  
  n_covariates <- length(covariate_names)
  
  # covariate simulation for BC trial using copula package
  
  # don't require copula for single covariate
  if (length(covariate_names) <= 1) {
    x_star <- 
      rnorm(n = N,
            mean = ald[[mean_names]],
            sd = ald[[sd_names]]) |> 
      matrix(ncol = 1, dimnames = list(NULL, covariate_names))
  
    return(x_star)
  }
  
  if (is.na(rho)) {
    # AC IPD pairwise correlations
    rho <- cor(ipd[, covariate_names])
  } else {
    # ensure in correct order
    rho <- rho[covariate_names, covariate_names]
  }
  # t_rho <- t(rho)  # extract along rows  ##TODO: isn't this symmetrical though?
  
  cor_ipd <- rho[lower.tri(rho, diag = FALSE)]
  # cor_ipd <- t_rho[lower.tri(t_rho, diag = FALSE)]
  
  cop <-
    copula::normalCopula(param = cor_ipd,
                         dim = n_covariates,
                         dispstr = "un")
  
  # aggregate BC covariate means & standard deviations
  mean_sd_margins <- list()
  
  for (i in covariate_names) {
    mean_sd_margins[[i]] <- list(mean = ald[[mean_names[i]]],
                                 sd = ald[[sd_names[i]]])
  }
  
  # sample covariates from approximate joint distribution using copula
  mvd <- copula::mvdc(
    copula = cop,
    margins = rep("norm", n_covariates),  # Gaussian marginals
    paramMargins = mean_sd_margins)

  # simulated BC pseudo-population
  x_star <- as.data.frame(copula::rMvdc(n = N, mvd))
  colnames(x_star) <- covariate_names
  
  ##TODO: this assume all covariate continuous
  ##      what about binary? threshold?
  x_star
}


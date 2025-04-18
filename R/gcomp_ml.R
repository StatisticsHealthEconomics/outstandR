
#' G-computation maximum likelihood bootstrap
#' 
#' Using bootstrap resampling, calculates the log odds ratio.
#'     
#' @param data Trial data 
#' @param indices Indices sampled from rows of `data`
#' @param formula Linear regression `formula` object
#' @param N Synthetic sample size for g-computation
#'
#' @return Mean difference in expected log-odds
#' @seealso [strategy_gcomp_ml()], [gcomp_ml_log_odds_ratio()]
#' @keywords internal
#' 
gcomp_ml.boot <- function(data, indices,
                          R, formula = NULL,
                          family, N = 1000, ald) {
  dat <- data[indices, ]
  gcomp_ml_means(formula, family, dat, ald) 
}


#' G-computation Maximum Likelihood mean outcome
#' 
#' @section Log-Odds Ratio: 
#' Marginal _A_ vs _C_ log-odds ratio (mean difference in expected log-odds)
#' estimated by transforming from probability to linear predictor scale.
#'
#' \eqn{\log(\hat{\mu}_A/(1 - \hat{\mu}_A)) - \log(\hat{\mu}_C/(1 - \hat{\mu}_C))}
#'
#' @param formula Linear regression `formula` object
#' @param family Family object
#' @template args-ipd
#' @template args-ald
#' @param N Synthetic sample size for g-computation
#'
#' @return Difference on relative treatment effect scale 
#' @seealso [strategy_gcomp_ml()], [gcomp_ml.boot()]
#' @importFrom copula normalCopula mvdc rMvdc
#' @importFrom stats predict glm
#' @keywords internal
#'
gcomp_ml_means <- function(formula,
                           family,
                           N = 1000,
                           ipd, ald) {
  
  if (!inherits(formula, "formula"))
    stop("formula argument must be of formula class.")
  
  x_star <- simulate_ALD_pseudo_pop(formula, ipd, ald, N)
  
  # outcome logistic regression fitted to IPD using maximum likelihood
  fit <- glm(formula = formula,
             family = family,
             data = ipd)
  
  # counterfactual datasets
  data.trtA <- data.trtC <- x_star
  
  treat_name <- get_treatment_name(formula)
  
  # # intervene on treatment while keeping set covariates fixed
  data.trtA[[treat_name]] <- 0  # all receive A
  data.trtC[[treat_name]] <- 1  # all receive C
  
  # predict counterfactual event probs, conditional on treatment/covariates
  hat.mu.A <- predict(fit, type="response", newdata=data.trtA)
  hat.mu.C <- predict(fit, type="response", newdata=data.trtC)
  
  # (marginal) mean probability prediction under A and C
  c(`0` = mean(hat.mu.A),
    `1` = mean(hat.mu.C))
}


#' Synthetic cohort using normal copula
#'
#' @param formula Linear regression `formula` object
#' @template args-ipd
#' @template args-ald
#' @param N Sample size
#' 
#' @importFrom copula normalCopula mvdc
#' @keywords internal
#' 
simulate_ALD_pseudo_pop <- function(formula,
                                    ipd, ald,
                                    N = 1000) {
  set.seed(1234)
  
  treat_name <- get_treatment_name(formula)
  covariate_names <- get_covariate_names(formula)
  
  # remove treatment
  covariate_names <- covariate_names[covariate_names != treat_name]
  
  mean_names <- get_mean_names(ald, covariate_names)
  
  sd_names <- get_sd_names(ald, covariate_names)
  
  n_covariates <- length(covariate_names)
  
  # covariate simulation for BC trial using copula package
  # AC IPD pairwise correlations
  rho <- cor(ipd[, covariate_names])
  t_rho <- t(rho)  # extract along rows
  
  cop <-
    copula::normalCopula(param = t_rho[lower.tri(t_rho)],
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


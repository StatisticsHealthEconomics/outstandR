
#' G-computation maximum likelihood bootstrap
#' 
#' Using bootstrap resampling, calculates the log odds ratio.
#'     
#' @param data Trial data 
#' @param indices Indices sampled from rows of `data`
#' @param formula Linear regression `formula` object
#'
#' @return Mean difference in expected log-odds
#' @seealso [strategy_gcomp_ml()], [gcomp_ml_log_odds_ratio()]
#' @keywords internal
#' 
gcomp_ml.boot <- function(data, indices,
                          R, formula = NULL,
                          family, ald, scale) {
  dat <- data[indices, ]
  gcomp_ml_ate(formula, family, dat, ald, scale) 
}


#' G-computation Maximum Likelihood ATE
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
#'
#' @return Difference of log-odds 
#' @seealso [strategy_gcomp_ml()], [gcomp_ml.boot()]
#' @importFrom copula normalCopula mvdc rMvdc
#' @importFrom stats predict glm
#' @keywords internal
#'
gcomp_ml_ate <- function(formula,
                         family,
                         ipd, ald,
                         scale) {
  
  if (!inherits(formula, "formula"))
    stop("formula argument must be of formula class.")
  
  x_star <- simulate_ALD_pseudo_pop(formula, ipd, ald)
  
  # outcome logistic regression fitted to IPD using maximum likelihood
  fit <- glm(formula,
             data = ipd,
             family = family)
  
  # counterfactual datasets
  data.trtA <- data.trtC <- x_star
  
  treat_name <- get_treatment_name(formula)
  
  # intervene on treatment while keeping set covariates fixed
  data.trtA[[treat_name]] <- 1  # all receive A
  data.trtC[[treat_name]] <- 0  # all receive C
  
  # predict counterfactual event probs, conditional on treatment/covariates
  hat.mu.A.i <- predict(fit, type="response", newdata=data.trtA)
  hat.mu.C.i <- predict(fit, type="response", newdata=data.trtC)
  
  hat.mu.A <- mean(hat.mu.A.i) # (marginal) mean probability prediction under A
  hat.mu.C <- mean(hat.mu.C.i) # (marginal) mean probability prediction under C
  
  calculate_ate(hat.mu.A, hat.mu.C, effect = scale)
}


#' Synthetic cohort using normal copula
#'
#' @param formula Linear regression `formula` object
#' @template args-ipd
#' @template args-ald
#' @importFrom copula normalCopula mvdc
#' @keywords internal
#' 
simulate_ALD_pseudo_pop <- function(formula, ipd, ald) {
  
  treat_name <- get_treatment_name(formula)
  covariate_names <- get_covariate_names(formula)
  
  # remove treatment
  covariate_names <- covariate_names[covariate_names != treat_name]
  
  ##TODO: this is a problem because we are only passing
  ##      the AC.IPD data and not BC.ALD
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
  
  # BC covariate means & standard deviations
  mean_sd_margins <- list()
  for (i in seq_len(n_covariates)) {
    mean_sd_margins[[i]] <- list(mean = ald[[mean_names[i]]],
                                 sd = ald[[sd_names[i]]])
  }
  
  # sample covariates from approximate joint distribution using copula
  mvd <- copula::mvdc(
    copula = cop,
    margins = rep("norm", n_covariates),  # Gaussian marginals
    paramMargins = mean_sd_margins)

  # simulated BC pseudo-population
  x_star <- as.data.frame(copula::rMvdc(n = 1000, mvd))
  colnames(x_star) <- covariate_names
  
  x_star
}


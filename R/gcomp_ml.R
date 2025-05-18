
#' G-computation maximum likelihood bootstrap
#' 
#' Using bootstrap resampling, calculates the relative treatment effect,
#' such as log odds ratio, log relative risk or risk difference.
#'     
#' @param data IPD trial data 
#' @param indices Indices sampled from rows of `data` for bootstrapping
#' @eval reg_args(include_formula = TRUE, include_family = TRUE)
#' @param rho A named square matrix of covariate correlations; default NA.
#' @param N Synthetic sample size for g-computation
#' @param ald Aggregate-level data for covariates.
#'
#' @return Relative treatment effect
#' @seealso [strategy_gcomp_ml()], [gcomp_ml_log_odds_ratio()]
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
                          ref_trt = "C", comp_trt = "A",
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
#' formula <- y ~ trt + age
#' ipd <- data.frame(tr = c("A", "C"), y = c(1, 0), age = c(30, 40))
#' ald <- data.frame()
#' simulate_ALD_pseudo_pop(formula, ipd, ald, trt_var = "trt", N = 1000)
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
  n_covariates <- length(covariate_names)
  
  ald_means <- dplyr::filter(ald, statistic == "mean", variable != "y")
  ald_sd <- dplyr::filter(ald, statistic == "sd", variable != "y")
  
  # same order as covariate names
  sd_values <- ald_sd$value[match(covariate_names, ald_sd$variable)]
  mean_values <- ald_means$value[match(covariate_names, ald_means$variable)]
  
  # covariate simulation for BC ALD trial using copula package
  
  # don't require copula for single covariate
  if (length(covariate_names) <= 1) {
    x_star <- 
      rnorm(n = N,
            mean = mean_values,
            sd = sd_values) |> 
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
  mean_sd_margins <- vector(mode = "list", length = n_covariates)
  
  for (i in 1:n_covariates) {
    mean_sd_margins[[i]] <- list(mean = mean_values[i],
                                 sd = sd_values[i])
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

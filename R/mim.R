
#' Multiple imputation marginalization (MIM)
#' 
#' @param strategy An object of class `strategy` created by functions such as 
#'   [strategy_maic()], [strategy_stc()], or [strategy_mim()]. 
#'   Contains modelling details like the formula and family.
#' @param analysis_params List of analysis parameters. Must contain `ipd` and `ald`.
#' @param ... Additional argument to pass to Stan model
#' 
#' @return A list containing:
#' * `means`: A list containing named vectors of posterior means (one per synthesis `n_imp`):
#'     * Comparator means.
#'     * Reference means.
#' * `model`: A list containing:
#'     * `fit`: The first-stage [rstanarm::stan_glm()] object.
#'     * `hats.v`: Vector of variance point estimates for each synthesis.
#'     * `n_imp`: Number of posterior prediction draws (syntheses).
#'     * `rho`, `N`, `stan_args`: Strategy and model parameters.
#'
#' @importFrom rstanarm posterior_predict stan_glm
#' @importFrom stats glm coef vcov as.formula setNames
#' @keywords internal
#' 
calc_mim <- function(strategy,
                     analysis_params,
                     ...) {
  
  default_stan_args <- list(
    algorithm = "sampling",
    chains = 2,
    iter = 2000,
    seed = 1234
  )
  
  stan_args <- modifyList(default_stan_args, list(...))
  
  ipd <- analysis_params$ipd
  ald <- analysis_params$ald 
  ref_trt <- analysis_params$ref_trt
  comp_trt <- analysis_params$ipd_comp
  formula <- strategy$outcome_model
  family <- strategy$family
  rho <- strategy$rho
  N <- strategy$N
  trt_var <- strategy$trt_var
  marginal_distns <- strategy$marginal_distns
  marginal_params <- strategy$marginal_params
  
  # Standardize number of imputations (syntheses) - default to 100 for MIM
  n_imp <- if (!is.null(strategy$n_imp)) strategy$n_imp else 100
  
  # SYNTHESIS STAGE ---
  
  # first-stage logistic regression model fitted to index RCT using MCMC (Stan)
  outcome_model <- do.call(rstanarm::stan_glm, c(
    list(formula = formula, data = ipd, family = family),
    stan_args
  ))
  
  # Resample ALD pseudo-population to capture covariate uncertainty
  x_star <- simulate_ALD_pseudo_pop(formula, ipd, ald, trt_var, rho, N,
                                    marginal_distns = marginal_distns,
                                    marginal_params = marginal_params)
  
  # create augmented target dataset
  target.comp <- target.ref <- x_star
  
  # Intervene safely, maintaining factor levels
  if (is.factor(ipd[[trt_var]])) {
    target.comp[[trt_var]] <- factor(comp_trt, levels = levels(ipd[[trt_var]]))
    target.ref[[trt_var]]  <- factor(ref_trt, levels = levels(ipd[[trt_var]]))
  } else {
    target.comp[[trt_var]] <- comp_trt  
    target.ref[[trt_var]]  <- ref_trt   
  }
  
  aug.target <- rbind(target.ref, target.comp)
  
  # set reference treatment as base level
  aug.target[[trt_var]] <- factor(aug.target[[trt_var]],
                                  levels = c(ref_trt, comp_trt))
  
  # complete syntheses by drawing binary outcomes
  # from their posterior predictive distribution
  # n_imp x iter
  y_star <- 
    rstanarm::posterior_predict(
      outcome_model, 
      newdata = aug.target,
      draws = n_imp)
  
  # ANALYSIS STAGE ---
  
  reg2.fits <- vector("list", n_imp)
  
  for(m in seq_len(n_imp)) {
    
    # fit second-stage regression to each synthesis using maximum-likelihood estimation
    data_m <- aug.target
    data_m$y <- y_star[m, ]
    # data_m$y <- as.numeric(y_star_m[1, ])
    
    # 4. Fit second-stage regression to the synthesis using ML estimation
    reg2.fits[[m]] <- stats::glm(stats::as.formula(paste("y ~", trt_var)), 
                                 data = data_m,
                                 family = family)
  }
  
  # treatment effect point estimates in each synthesis
  coef_fit <- do.call(rbind, lapply(reg2.fits, function(fit) stats::coef(fit)))
  
  # safer than trt_var in case of factor level append
  coef_names <- names(stats::coef(reg2.fits[[1]]))
  treat_coef_name <- 
    grep(pattern = paste0("^", trt_var, "[^:]*$"), coef_names, value = TRUE)
  
  # point estimates for the variance in each synthesis
  hats.v <- unlist(lapply(reg2.fits, function(fit)
    stats::vcov(fit)[treat_coef_name, treat_coef_name]))
  
  mean_ref <- family$linkinv(coef_fit[, 1])     # probability for reference
  mean_comp <- family$linkinv(coef_fit[, 1] + coef_fit[, treat_coef_name])  # probability for comparator
  
  # Harmonized Output Naming
  means_list <- stats::setNames(list(mean_comp, mean_ref), 
                                c(comp_trt, ref_trt))
  
  point_est_list <- stats::setNames(list(mean(mean_comp), mean(mean_ref)), 
                                    c(comp_trt, ref_trt))
  
  list(
    means = means_list,
    point_estimates = point_est_list,
    model = list(
      fit = outcome_model,
      hats.v = hats.v,
      n_imp = n_imp,
      rho = rho,
      N = N,
      stan_args = stan_args)
  )
}

#' Wald-type interval estimates
#' 
#' Constructed using t-distribution with nu degrees of freedom.
#'
#' @param n_imp Number of syntheses used in analysis stage (high for low Monte Carlo error)
#' @param bar.v "within" variance (average of variance point estimates)
#' @param b "between" variance (sample variance of point estimates)
#' @return Numeric value of Wald-type interval estimates.
#' @keywords internal
#' 
wald_type_interval <- function(n_imp, bar.v, b) {
  (n_imp - 1) * (1 + bar.v / ((1 + 1 / n_imp) * b)) ^ 2
}

#' Variance estimate by pooling
#' 
#' Use combining rules to estimate.
#' 
#' @param n_imp Number of syntheses used in analysis stage (high for low Monte Carlo error)
#' @param bar.v "within" variance (average of variance point estimates)
#' @param b "between" variance (sample variance of point estimates)
#' @return Numeric value of variance estimate using pooling.
#' @keywords internal
#' 
var_by_pooling <- function(n_imp, bar.v, b) {
  (1 + (1 / n_imp)) * b - bar.v
}



#' Multiple imputation marginalization (MIM)
#' 
#' @param strategy An object of class `strategy` created by functions such as 
#'   [strategy_maic()], [strategy_stc()], or [strategy_mim()]. 
#'   Contains modelling details like the formula and family.
#' @eval study_data_args(include_ipd = TRUE, include_ald = TRUE)
#' @param ref_trt Reference treatment. String.
#' @param comp_trt Comparator treatment. String.
#' @param ... Additional argument to pass to Stan model
#' 
#' @return Named list:
#' \describe{
#'    \item{mean_comp}{Comparator mean numeric value}
#'    \item{mean_ref}{Reference mean numeric value}
#'    \item{hats.v}{Point estimates for the variance in each synthesis}
#'    \item{M}{Number of posterior prediction draws}
#' }
#' 
#' @importFrom rstanarm posterior_predict stan_glm
#' @importFrom tibble tibble lst
#' @keywords internal
#' 
calc_mim <- function(strategy,
                     ipd, ald, 
                     ref_trt,
                     comp_trt, ...) {
  
  formula <- strategy$formula
  family <- strategy$family
  rho <- strategy$rho
  N <- strategy$N
  trt_var <- strategy$trt_var
  
  x_star <- simulate_ALD_pseudo_pop(formula, ipd, ald, trt_var, rho, N)
  
  # SYNTHESIS STAGE ---
  
  # first-stage logistic regression model fitted to index RCT using MCMC (Stan)
  outcome_model <- stan_glm(
    formula = formula,
    data = ipd,
    family = family,
    algorithm = "sampling", ...)
  
  # create augmented target dataset
  target.comp <- target.ref <- x_star
  target.comp[[trt_var]] <- comp_trt
  target.ref[[trt_var]] <- ref_trt
  
  aug.target <- rbind(target.ref, target.comp)
  
  # set reference treatment as base level
  aug.target[[trt_var]] <- factor(aug.target[[trt_var]],
                                  levels = c(ref_trt, comp_trt))
  
  # complete syntheses by drawing binary outcomes
  # from their posterior predictive distribution
  y_star <-
    rstanarm::posterior_predict(
      outcome_model, newdata = aug.target)
  
  # ANALYSIS STAGE ---
  
  M <- nrow(y_star)
  
  # fit second-stage regression to each synthesis using maximum-likelihood estimation
  reg2.fits <- lapply(1:M, function(m) {
    data_m <- aug.target
    data_m$y <- y_star[m, ]
    glm(as.formula(paste("y ~", trt_var)), data = data_m, family = family)
  })

  # treatment effect point estimates in each synthesis
  coef_fit <- do.call(rbind, lapply(reg2.fits, function(fit) coef(fit)))
  
  # safer than trt_var in case of factor level append
  coef_names <- names(coef(reg2.fits[[1]]))
  treat_coef_name <- grep(pattern = paste0("^", trt_var, "[^:]*$"), coef_names, value = TRUE)
  
  ##TODO: how to transform this to the prob scale?
  # point estimates for the variance in each synthesis
  hats.v <- unlist(lapply(reg2.fits,
                          function(fit)
                            vcov(fit)[treat_coef_name, treat_coef_name]))
  
  mean_ref <- family$linkinv(coef_fit[, 1])     # probability for reference
  mean_comp <- family$linkinv(coef_fit[, 1] + coef_fit[, treat_coef_name])  # probability for comparator
  
  list(
    means = list(
      A = mean_comp,
      C = mean_ref),
    model = list(
      fit = outcome_model,
      hats.v = hats.v,
      M = M,
      rho = rho,
      N = N,
      stan_args = list(...))
  )
}

#' Wald-type interval estimates
#' 
#' Constructed using t-distribution with nu degrees of freedom.
#'
#' @param M Number of syntheses used in analysis stage (high for low Monte Carlo error)
#' @param bar.v "within" variance (average of variance point estimates)
#' @param b "between" variance (sample variance of point estimates)
#' @return Numeric value of Wald-type interval estimates.
#' @keywords internal
#' 
wald_type_interval <- function(M, bar.v, b) {
  (M - 1) * (1 + bar.v / ((1 + 1 / M) * b)) ^ 2
}

#' Variance estimate by pooling
#' 
#' Use combining rules to estimate.
#' 
#' @param M Number of syntheses used in analysis stage (high for low Monte Carlo error)
#' @param bar.v "within" variance (average of variance point estimates)
#' @param b "between" variance (sample variance of point estimates)
#' @return Numeric value of variance estimate using pooling.
#' @keywords internal
#' 
var_by_pooling <- function(M, bar.v, b) {
  (1 + (1 / M)) * b - bar.v
}


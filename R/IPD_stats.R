
# create S3 class for each approach

#' @rdname strategy
#' 
#' @section Matching-adjusted indirect comparison (MAIC):
#' 
#' MAIC is a form of non-parametric likelihood reweighting method
#' which allows the propensity score logistic 
#' regression model to be estimated without IPD in the _AC_ population.
#' The mean outcomes \eqn{\mu_{t(AC)}} on treatment \eqn{t = A,B} in the _AC_
#' target population are estimated by taking a weighted average of the
#' outcomes \eqn{Y} of the \eqn{N} individuals in arm \eqn{t} of the _AB_ population
#' 
#' Used to compare marginal treatment effects where there are cross-trial
#' differences in effect modifiers and limited patient-level data.
#' 
#' \deqn{
#' \hat{Y}_{} = \frac{\sum_{i=1}^{N} Y_{it(AB)} w_{it}}{\sum _{i=1}^{N} w_{it}}
#' }
#' where the weight \eqn{w_{it}} assigned to the \eqn{i}-th individual receiving treatment
#' \eqn{t} is equal to the odds of being enrolled in the _AC_ trial vs the _AB_ trial.
#'
#'
#' The default formula is
#' \deqn{
#'  y = X_3 + X_4 + \beta_{t}X_1 + \beta_{t}X_2
#' }
#' 
#' @param R The number of resamples used for the non-parametric bootstrap
#' @template args-ald
#'
#' @return `maic` class object
#' @export
#'
strategy_maic <- function(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"),
                          R = 1000,
                          ald) {
  
  if (class(formula) != "formula")
    stop("formula argument must be of formula class.")
  
  default_args <- formals()
  args <- as.list(match.call())[-1]
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "maic", args))
}

#' @rdname strategy
#' 
#' @section Simulated treatment comparison (STC):
#' Outcome regression-based method which targets a conditional treatment effect.
#' STC is a modification of the covariate adjustment method.
#' An outcome model is fitted using IPD in the _AB_ trial
#' 
#' \deqn{
#' g(\mu_{t(AB)}(X)) = \beta_0 + \beta_1^T X + (\beta_B + \beta_2^T X^{EM}) I(t=B)
#' }
#' where \eqn{\beta_0} is an intercept term, \eqn{\beta_1} is a vector of coefficients for
#' prognostic variables, \eqn{\beta_B} is the relative effect of treatment _B_ compared
#' to _A_ at \eqn{X=0}, \eqn{\beta_2} is a vector of coefficients for effect
#' modifiers \eqn{X^{EM}} subvector of the full covariate vector \eqn{X}), and
#' \eqn{\mu_{t(AB)}(X)} is the expected outcome of an individual assigned
#' treatment \eqn{t} with covariate values \eqn{X} which is transformed onto a
#' chosen linear predictor scale with link function \eqn{g(\cdot)}.
#' 
#' The default formula is
#' \deqn{
#'  y = X_3 + X_4 + \beta_t(X_1 - \bar{X_1}) + \beta_t(X_2 - \bar{X2})
#' }
#' 
#' @return `stc` class object
#' @export
# 
strategy_stc <- function(formula =
                           as.formula("y ~ X3 + X4 + trt*I(X1 - mean(X1)) + trt*I(X2 - mean(X2))")) {
  
  if (class(formula) != "formula")
    stop("formula argument must be of formula class.")
  
  default_args <- formals()
  args <- as.list(match.call())[-1]
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "stc", args))
}

#' @rdname strategy
#' 
#' @section G-computation maximum likelihood:
#'
#' G-computation marginalizes the conditional estimates by separating the regression modelling
#' from the estimation of the marginal treatment effect for _A_ versus _C_.
#' First, a regression model of the observed outcome \eqn{y} on the covariates \eqn{x} and treatment \eqn{z} is fitted to the _AC_ IPD:
#' 
#' \deqn{
#' g(\mu_n) = \beta_0 + \boldsymbol{x}_n \boldsymbol{\beta_1} + (\beta_z + \boldsymbol{x_n^{EM}} \boldsymbol{\beta_2}) \mbox{I}(z_n=1)
#' }
#' In the context of G-computation, this regression model is often called the “Q-model.”
#' Having fitted the Q-model, the regression coefficients are treated as nuisance parameters.
#' The parameters are applied to the simulated covariates \eqn{x*} to predict hypothetical outcomes
#' for each subject under both possible treatments. Namely, a pair of predicted outcomes,
#' also called potential outcomes, under _A_ and under _C_, is generated for each subject.
#'  
#' By plugging treatment _C_ into the regression fit for every simulated observation,
#' we predict the marginal outcome mean in the hypothetical scenario in which all units are under treatment _C_:
#'
#' \deqn{
#' \hat{\mu}_0 = \int_{x^*} g^{-1} (\hat{\beta}_0 + x^* \hat{\beta}_1 ) p(x^*) dx^*
#' }
#' To estimate the marginal or population-average treatment effect for A versus C in the linear predictor scale,
#' one back-transforms to this scale the average predictions, taken over all subjects on the natural outcome scale,
#' and calculates the difference between the average linear predictions:
#'  
#' \deqn{
#' \hat{\Delta}^{(2)}_{10} = g(\hat{\mu}_1) - g(\hat{\mu}_0)
#' }
#' 
#' The default formula is
#' \deqn{
#'  y = X_3 + X_4 + \beta_{t}X_1 + \beta_{t}X_2
#' }
#'
#' @param R The number of resamples used for the non-parametric bootstrap
#' 
#' @return `gcomp_ml` class object
#' @export
#'
strategy_gcomp_ml <- function(formula =
                                as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"),
                              R = 1000) {
  
  if (class(formula) != "formula")
    stop("formula argument must be of formula class.")
  
  default_args <- formals()
  args <- as.list(match.call())[-1]
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "gcomp_ml", args))
}

#' @rdname strategy
#' 
#' @section G-computation Bayesian:
#' The difference between Bayesian G-computation and its maximum-likelihood
#' counterpart is in the estimated distribution of the predicted outcomes. The
#' Bayesian approach also marginalizes, integrates or standardizes over the
#' joint posterior distribution of the conditional nuisance parameters of the
#' outcome regression, as well as the joint covariate distribution.
#'
#' Draw a vector of size \eqn{N*} of predicted outcomes \eqn{y*z} under each set
#' intervention \eqn{z* \in \{0, 1\}} from its posterior predictive distribution
#' under the specific treatment. This is defined as \eqn{p(y*_{z*} |
#' \mathcal{D}_{AC}) = \int_{\beta} p(y*_{z*} | \beta) p(\beta | \mathcal{D}_{AC}) d\beta}
#' where \eqn{p(\beta | \mathcal{D}_{AC})} is the
#' posterior distribution of the outcome regression coefficients \eqn{\beta},
#' which encode the predictor-outcome relationships observed in the _AC_ trial IPD.
#' 
#' This is given by:
#'
#' \deqn{
#' p(y*_{z*} | \mathcal{D}_{AC}) = \int_{x*} p(y* | z*, x*, \mathcal{D}_{AC}) p(x* | \mathcal{D}_{AC}) dx*
#' }
#' 
#' \deqn{
#' = \int_{x*} \int_{\beta} p(y* | z*, x*, \beta) p(x* | \beta) p(\beta | \mathcal{D}_{AC}) d\beta dx*
#' }
#' In practice, the integrals above can be approximated numerically, using full Bayesian
#' estimation via Markov chain Monte Carlo (MCMC) sampling.
#' 
#' The default formula is
#' \deqn{
#'  y = X_3 + X_4 + \beta_{t}X_1 + \beta_{t}X_2
#' }
#'
#' @return `gcomp_stan` class object
#' @export
#'
strategy_gcomp_stan <- function(formula =
                                  as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")) {
  
  if (class(formula) != "formula")
    stop("formula argument must be of formula class.")
  
  default_args <- formals()
  args <- as.list(match.call())[-1]
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "gcomp_stan", args))
}

#' @name strategy
#' @title New strategy objects
#' 
#' @description
#' Create a type of strategy class for each modelling approach.
#'
#' @param strategy Class name from `strategy_maic`, `strategy_stc`, `strategy_gcomp_ml`, `strategy_gcomp_stan`
#' @param formula Linear regression `formula` object
#' @param ... Additional arguments
#'
#' @export
#'
new_strategy <- function(strategy, ...) {
  structure(list(...), class = c(strategy, "strategy", "list"))
}


#' @title Calculate the difference between treatments using all evidence
#' 
#' @description
#' This is the main, top-level wrapper for `hat_Delta_stats()`.
#' Methods taken from
#' \insertCite{RemiroAzocar2022}{mimR}.
#' 
#' @param AC.IPD Individual-level patient data. Suppose between studies _A_ and _C_.
#' @param BC.ALD Aggregate-level data. Suppose between studies _B_ and _C_. 
#' @param strategy Computation strategy function. These can be
#'    `strategy_maic()`, `strategy_stc()`, `strategy_gcomp_ml()` and  `strategy_gcomp_stan()`
#' @param CI Confidence interval; between 0,1
#' @param ... Additional arguments
#' @return List of length 3 of statistics as a `mimR` class object.
#'   Containing statistics between each pair of treatments.
#'   These are the mean contrasts, variances and confidence intervals,
#'   respectively.
#' @importFrom Rdpack reprompt
#' 
#' @references
#' \insertRef{RemiroAzocar2022}{mimR}
#' 
#' @export
#' @examples
#' 
#' data(AC_IPD)  # AC patient-level data
#' data(BC_ALD)  # BC aggregate-level data
#' 
#' # matching-adjusted indirect comparison
#' hat_Delta_stats_maic <- hat_Delta_stats(AC_IPD, BC_ALD, strategy = strategy_maic())
#' 
#' # simulated treatment comparison
#' hat_Delta_stats_stc <- hat_Delta_stats(AC_IPD, BC_ALD, strategy = strategy_stc())
#' 
#' # G-computation with maximum likelihood
#' # hat_Delta_stats_gcomp_ml <- hat_Delta_stats(AC_IPD, BC_ALD, strategy = strategy_gcomp_ml())
#' 
#' # G-computation with Bayesian inference
#' hat_Delta_stats_gcomp_stan <- hat_Delta_stats(AC_IPD, BC_ALD, strategy = strategy_gcomp_stan())
#' 
hat_Delta_stats <- function(AC.IPD, BC.ALD, strategy, CI = 0.95, ...) {

  if (CI <= 0 || CI >= 1) stop("CI argument must be between 0 and 1.")
  
  ##TODO: as method instead?
  if (!inherits(strategy, "strategy"))
    stop("strategy argument must be of a class strategy.")
  
  AC_hat_Delta_stats <- IPD_stats(strategy, ipd = AC.IPD, ald = BC.ALD, ...) 
  BC_hat_Delta_stats <- ALD_stats(ald = BC.ALD) 
  
  upper <- 0.5 + CI/2
  ci_range <- c(1-upper, upper)
  
  contrasts <- list(
    AB = AC_hat_Delta_stats$mean - BC_hat_Delta_stats$mean,
    AC = AC_hat_Delta_stats$mean,
    BC = BC_hat_Delta_stats$mean)
  
  contrast_variances <- list(
    AB = AC_hat_Delta_stats$var + BC_hat_Delta_stats$var,
    AC = AC_hat_Delta_stats$var,
    BC = BC_hat_Delta_stats$var)
  
  contrast_ci <- list(
    AB = contrasts$AB + qnorm(ci_range)*as.vector(sqrt(contrast_variances$AB)),
    AC = contrasts$AC + qnorm(ci_range)*as.vector(sqrt(contrast_variances$AC)),
    BC = contrasts$BC + qnorm(ci_range)*as.vector(sqrt(contrast_variances$BC)))
  
  stats <- list(contrasts = contrasts,
                variances = contrast_variances,
                CI = contrast_ci)
  
  structure(stats,
            CI = CI,
            class = c("mimR", class(stats)))
}


#' @name IPD_stats
#' @title Calculate individual-level patient data statistics
#' 
#' @description
#' Separate methods for each approach
#' MAIC, STC, G-computation via MLE or Bayesian inference.
#' 
#' @param strategy A list corresponding to different approaches
#' @template args-ipd
#' @template args-ald
#' @param ... Additional arguments
#' 
#' @return Mean and variance values
#' @export
#' 
IPD_stats <- function(strategy, ipd, ald, ...)
  UseMethod("IPD_stats", strategy)


#' @rdname IPD_stats
#' @export
IPD_stats.default <- function(...) {
  strategy_classes <- sub("IPD_stats\\.", "", methods(IPD_stats)[-1])
  avail_strategies <- paste0("strategy_", strategy_classes, "()", collapse = ", ")
  stop(paste0("strategy not available. Select from ", avail_strategies))
}


#' @rdname IPD_stats
#' @section Matching-adjusted indirect comparison statistics:
#' Marginal _A_ vs _C_ treatment effect estimates
#' using bootstrapping sampling.
#'
#' @export
#' 
IPD_stats.maic <- function(strategy,
                           ipd, ald) {
  # maic.boot(data = data,
  #           indices = 1:nrow(data),
  #           formula = strategy$formula,
  #           dat_ALD = strategy$dat_ALD)
  
  maic_boot <- boot::boot(data = ipd,
                          statistic = maic.boot,
                          R = strategy$R,
                          formula = strategy$formula,
                          ald = strategy$dat_ALD)
  
  list(mean = mean(maic_boot$t),
       var = var(maic_boot$t))
}


#' @rdname IPD_stats
#' @section Simulated treatment comparison statistics: 
#' IPD from the _AC_ trial are used to fit a regression model describing the
#' observed outcomes \eqn{y} in terms of the relevant baseline characteristics `x` and
#' the treatment variable \eqn{z}.
#' 
#' @export
#' 
IPD_stats.stc <- function(strategy,
                          ipd, ald) {
  fit <- glm(strategy$formula,
             data = ipd,
             family = binomial)
  
  treat_nm <- get_treatment_name(strategy$formula)
  
  # fitted treatment coefficient is relative A vs C conditional effect
  list(mean = coef(fit)[treat_nm],
       var = vcov(fit)[treat_nm, treat_nm])
}


#' @rdname IPD_stats
#' @section G-computation maximum likelihood statistics:
#' Compute a non-parametric bootstrap with \eqn{R=1000} resamples.
#'
#' @export
#'
IPD_stats.gcomp_ml <- function(strategy,
                               ipd, ald) {

  AC_maic_boot <- boot::boot(data = ipd,
                             statistic = gcomp_ml.boot,
                             R = strategy$R,
                             formula = strategy$formula)
  
  list(mean = mean(AC_maic_boot$t),
       var = var(AC_maic_boot$t))
}


#' @rdname IPD_stats
#' @section G-computation Bayesian statistics:
#' Using Stan, compute marginal log-odds ratio for _A_ vs _C_ for each MCMC sample
#' by transforming from probability to linear predictor scale
#'
#' @export
#'
IPD_stats.gcomp_stan <- function(strategy,
                                 ipd, ald) {
  
  ppv <- gcomp_stan(formula = strategy$formula,
                    ipd = ipd, ald = ald)
  
  hat.delta.AC <-
    qlogis(rowMeans(ppv$y.star.A)) - qlogis(rowMeans(ppv$y.star.C))
  
  list(mean = mean(hat.delta.AC),
       var = var(hat.delta.AC))
} 


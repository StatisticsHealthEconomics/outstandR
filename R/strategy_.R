
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
#' outcomes \eqn{Y} of the \eqn{N} individuals in arm \eqn{t} of the _AB_ population.
#' 
#' Used to compare marginal treatment effects where there are cross-trial
#' differences in effect modifiers and limited patient-level data.
#' 
#' \deqn{
#' \hat{Y}_{} = \frac{\sum_{i=1}^{N} Y_{it(AB)} w_{it}}{\sum_{i=1}^{N} w_{it}}
#' }
#' where the weight \eqn{w_{it}} assigned to the \eqn{i}-th individual receiving treatment
#' \eqn{t} is equal to the odds of being enrolled in the _AC_ trial vs the _AB_ trial.
#' 
#' @param R The number of resamples used for the non-parametric bootstrap
#' @return `maic` class object
#' 
#' @importFrom utils modifyList
#' @export
#'
strategy_maic <- function(formula = NULL,
                          family = gaussian(link = "identity"),
                          R = 1000L) {
  check_formula(formula)
  check_family(family)
   
  if (R <= 0 || R %% 1 != 0) {
    stop("R not positive whole number.")
  }
  
  default_args <- formals()
  args <- c(formula = formula, as.list(match.call())[-c(1,2)])
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
#' @return `stc` class object
#' @importFrom utils modifyList
#' @export
# 
strategy_stc <- function(formula = NULL,
                         family = gaussian(link = "identity")) {
  check_formula(formula)
  check_family(family)
  
  default_args <- formals()
  args <- c(formula = formula, as.list(match.call())[-c(1,2)])
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
#' To estimate the marginal or population-average treatment effect for _A_ versus _C_ in the linear predictor scale,
#' one back-transforms to this scale the average predictions, taken over all subjects on the natural outcome scale,
#' and calculates the difference between the average linear predictions:
#'  
#' \deqn{
#' \hat{\Delta}^{(2)}_{10} = g(\hat{\mu}_1) - g(\hat{\mu}_0)
#' }
#'
#' @param rho A named square matrix of covariate correlations; default NA.
#' @param R The number of resamples used for the non-parametric bootstrap
#' @param N Synthetic sample size for g-computation
#' 
#' @return `gcomp_ml` class object
#' @importFrom utils modifyList
#' @export
#'
strategy_gcomp_ml <- function(formula = NULL,
                              family = gaussian(link = "identity"),
                              rho = NA,
                              R = 1000L,
                              N = 1000L) {
  check_formula(formula)
  check_family(family)
  
  if (R <= 0 || R %% 1 != 0) {
    stop("R not positive whole number.")
  }
  if (N <= 0 || N %% 1 != 0) {
    stop("N not positive whole number.")
  }
  
  default_args <- formals()
  args <- c(formula = formula, as.list(match.call())[-c(1,2)])
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
#' p(y*_{z*} \mid \mathcal{D}_{AC}) = \int_{x*} p(y* \mid z*, x*, \mathcal{D}_{AC}) p(x* \mid \mathcal{D}_{AC}) dx*
#' }
#' 
#' \deqn{
#' = \int_{x*} \int_{\beta} p(y* \mid z*, x*, \beta) p(x* \mid \beta) p(\beta \mid \mathcal{D}_{AC}) d\beta dx*
#' }
#' In practice, the integrals above can be approximated numerically, using full Bayesian
#' estimation via Markov chain Monte Carlo (MCMC) sampling.
#' @return `gcomp_stan` class object
#' @importFrom utils modifyList
#' @export
#'
strategy_gcomp_stan <- function(formula = NULL,
                                family = gaussian(link = "identity")) {
  check_formula(formula)
  check_family(family)
  
  default_args <- formals()
  args <- c(formula = formula, as.list(match.call())[-c(1,2)])
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "gcomp_stan", args))
}

#' @rdname strategy
#' 
#' @section Multiple imputation marginalization (MIM):
#' @return `mim` class object
#' @importFrom utils modifyList
#' @export
# 
strategy_mim <- function(formula = NULL,
                         family = gaussian(link = "identity")) {
  check_formula(formula)
  check_family(family)
  
  default_args <- formals()
  args <- c(formula = formula, as.list(match.call())[-c(1,2)])
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "mim", args))
}

#' @name strategy
#' @title New strategy objects
#' 
#' @description
#' Create a type of strategy class for each modelling approach.
#'
#' @param strategy Class name from `strategy_maic`, `strategy_stc`, `strategy_gcomp_ml`, `strategy_gcomp_stan`
#' @param formula Linear regression `formula` object
#' @param family `family` object from the `stats` library
#' @param ... Additional arguments
#'
#' @export
#'
new_strategy <- function(strategy, ...) {
  structure(list(...), class = c(strategy, "strategy", "list"))
}

#
is_family <- function(obj) inherits(obj, "family")

#
check_family <- function(obj) {
  if (!is_family(obj)) {
    stop("family must be a family object")
  }
}

# create S3 class for each approach ---

#' @rdname strategy
#' 
#' @section Matching-adjusted indirect comparison (MAIC):
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
#' @param trt_var Treatment variable name
#' @param n_boot The number of resamples used for the non-parametric bootstrap
#' @param moments Integer. The number of moments of the covariates to balance. 
#'   Setting \code{moments = 2} includes both the original variables and their 
#'   squared terms, which effectively balances both the means and the variances.
#'   Default to 1.
#' @param int Logical. If \code{TRUE}, includes two-way interactions between all 
#'   covariates in the balancing model to effectively balance their covariances.
#'    Default \code{FALSE}
#' @param verbatim Logical. Output to console during running.
#' @return `maic` class object
#' 
#' @importFrom utils modifyList
#' @export
strategy_maic <- function(formula = NULL,
                          family = gaussian(link = "identity"),
                          trt_var = NULL,
                          n_boot = 1000L,
                          moments = 1,
                          int = FALSE,
                          verbatim = TRUE) {
  
  # parse formula depending on whether list or legacy formula
  if (!is.list(formula)) {
    if (verbatim) {
      message("Note: Using legacy 'formula' argument.")
      message(paste("--> Analysis Model:", deparse(formula)))
    }
    
    # guess trt_var safely because legacy formula exists
    trt_var <- get_treatment_name(formula, trt_var)
    
    # internally convert to the balance model (stripping 'y' and 'trt')
    rhs_vars <- all.vars(delete.response(terms(formula)))
    balance_vars <- setdiff(rhs_vars, trt_var)
    
    balance_model <- as.formula(paste("~", paste(balance_vars, collapse = " + ")))
    outcome_model <- formula
    
    if (verbatim) {
      message(paste("--> Inferred Balance Model: ~", paste(balance_vars, collapse = " + ")))
      message("    (Balancing on means of these covariates by default)")
    }
    
  } else {
    outcome_model <- formula$outcome_model
    balance_model <- formula$balance_model
    
    # Handle case where outcome_model is completely omitted from list
    if (is.null(outcome_model)) {
      if (is.null(trt_var)) {
        trt_var <- "trt"
        
        if (verbatim) {
          message("Outcome model and trt_var not provided. Defaulting trt_var to '", trt_var, "'.")
        }
      }
      
      outcome_model <- as.formula(paste("y ~", trt_var))
      
      if (verbatim) {
        message("Outcome model missing. Defaulting to: ", deparse(outcome_model))
      }      
    } else {
      trt_var <- get_treatment_name(outcome_model, trt_var)
    }
    
    # Handle missing balance_model with a delayed function
    if (is.null(balance_model)) {
      balance_model <- function(ald) {
        # Extract all unique covariates from the ALD, dropping NAs and outcome terms
        covars <- unique(ald$variable)
        covars <- covars[!is.na(covars) & !(covars %in% c("y", "N"))]
        
        if (length(covars) == 0) return(NULL) # Failsafe if ALD has no covariates
        
        as.formula(paste("~", paste(covars, collapse = " + ")))
      }
      if (verbatim) {
        message("Balance model not provided. Default to a linear sum of all covariates found in the ALD.")
      }
    }
  }
  # Ensure MAIC uses an unadjusted outcome model to prevent estimating conditional effects.
  resp_var <- all.vars(outcome_model)[1]
  out_vars <- all.vars(outcome_model)
  
  if (length(setdiff(out_vars, c(resp_var, trt_var))) > 0) {
    if (verbatim) {
      warning("Covariates detected in the MAIC outcome model. To ensure the estimation of compatible marginal treatment effects, MAIC requires an unadjusted outcome model. The outcome model is being automatically overridden to '", resp_var, " ~ ", trt_var, "'.", call. = FALSE)
    }
    outcome_model <- as.formula(paste(resp_var, "~", trt_var))
  }
  
  check_formula(outcome_model, trt_var)
  
  if (!is.null(balance_model) && !is.function(balance_model)) {
    check_balance_formula(balance_model, trt_var) 
  }
  
  check_family(family)
  
  if (n_boot <= 0 || n_boot %% 1 != 0) {
    stop("n_boot not positive whole number.")
  }
  
  args <- list(balance_model = balance_model,
               outcome_model = outcome_model,
               family = family,
               trt_var = trt_var,
               n_boot = n_boot,
               moments = moments,
               int = int)
  
  do.call(new_strategy, c(strategy = "maic", args))
}

#' Simulated treatment comparison (STC)
#' 
#' `r lifecycle::badge("deprecated")`
#' 
#' `strategy_stc()` was deprecated in outstandR version 1.X.X. 
#' We recommend using G-computation (`strategy_gcomp_ml()`) as a more robust 
#' alternative for this type of analysis.
#' 
#' Outcome regression-based method which targets a conditional treatment effect.
#' STC is a modification of the covariate adjustment method.
#' An outcome model is fitted using IPD in the _AB_ trial. For example,
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
#' @param formula Model formula.
#' @param family Model family.
#' @param trt_var Treatment variable name.
#'
#' @return `stc` class object
#' @importFrom utils modifyList
#' @export
# 
strategy_stc <- function(formula = NULL,
                         family = gaussian(link = "identity"),
                         trt_var = NULL) {
  lifecycle::deprecate_warn(
    when = "1.X.X",                           # version it is deprecated in
    what = "outstandR::strategy_stc()",       # function being deprecated
    with = "outstandR::strategy_gcomp_ml()"   # suggested alternative (optional)
  )
  
  # back-compatibility
  if (!is.list(formula)) {
    outcome_model <- formula
  } else {
    outcome_model <- formula$outcome_model
  }
  
  check_formula(outcome_model, trt_var)
  check_family(family)
  
  args <- list(outcome_model = outcome_model,
               family = family,
               trt_var = get_treatment_name(outcome_model, trt_var))
  
  do.call(new_strategy, c(strategy = "stc", args))
}

#' @rdname strategy
#' 
#' @section G-computation maximum likelihood:
#' G-computation marginalizes the conditional estimates by separating the regression modelling
#' from the estimation of the marginal treatment effect for _A_ versus _C_.
#' For example, a regression model of the observed outcome \eqn{y} on the covariates \eqn{x} and
#' treatment \eqn{z} is fitted to the _AC_ IPD:
#' 
#' \deqn{
#' g(\mu_n) = \beta_0 + \boldsymbol{x}_n \boldsymbol{\beta_1} + (\beta_z + \boldsymbol{x_n^{EM}} \boldsymbol{\beta_2}) \mbox{I}(z_n=1)
#' }
#' In the context of G-computation, this regression model is called the “Q-model".
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
#' @param trt_var Treatment variable name; string
#' @param rho A named square matrix of covariate correlations; default NA
#' @param marginal_distns Marginal distributions names; vector default NA
#' @param marginal_params Marginal distributions parameters; list of lists, default NA
#' @param n_boot The number of resamples used for the non-parametric bootstrap; integer
#' @param N Synthetic sample size for g-computation; integer; integer
#' 
#' @return `gcomp_ml` class object
#' @importFrom utils modifyList
#' @seealso [strategy_gcomp_bayes()]
#' @export
#'
strategy_gcomp_ml <- function(formula = NULL,
                              family = gaussian(link = "identity"),
                              trt_var = NULL,
                              rho = NA,
                              marginal_distns = NA,
                              marginal_params = NA,
                              n_boot = 1000L,
                              N = 1000L) {
  # back-compatibility
  if (!is.list(formula)) {
    outcome_model <- formula
  } else {
    outcome_model <- formula$outcome_model
  }
  
  check_formula(outcome_model, trt_var)
  check_family(family)
  check_distns(outcome_model, marginal_distns, marginal_params)
  check_rho(rho)
  
  if (n_boot <= 0 || n_boot %% 1 != 0) {
    stop("n_boot not positive whole number.")
  }
  if (N <= 0 || N %% 1 != 0) {
    stop("N not positive whole number.")
  }
  
  args <- list(outcome_model = outcome_model,
               family = family,
               rho = rho,
               trt_var = get_treatment_name(outcome_model, trt_var),
               marginal_distns = marginal_distns,
               marginal_params = marginal_params,
               n_boot = n_boot,
               N = N)
  
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
#' 
#' @param trt_var Treatment variable name
#' @param rho A named square matrix of covariate correlations; default NA
#' @param marginal_distns Marginal distributions names; vector default NA.
#'    Available distributions are given in stats::Distributions. See [copula::Mvdc()] for details
#' @param marginal_params Marginal distributions parameters; list of lists, default NA. See [copula::Mvdc()] for details
#' @param N Synthetic sample size for g-computation
#' 
#' @return `gcomp_bayes` class object
#' @importFrom utils modifyList
#' @seealso [strategy_gcomp_ml()] [copula::Mvdc()]
#' @export
#'
strategy_gcomp_bayes <- function(formula = NULL,
                                 family = gaussian(link = "identity"),
                                 trt_var = NULL,
                                 rho = NA,
                                 marginal_distns = NA,
                                 marginal_params = NA,
                                 N = 1000L) {
  # back-compatibility
  if (!is.list(formula)) {
    outcome_model <- formula
  } else {
    outcome_model <- formula$outcome_model
  }
  
  check_formula(outcome_model, trt_var)
  check_family(family)
  check_distns(outcome_model, marginal_distns, marginal_params)
  check_rho(rho)
  
  if (N <= 0 || N %% 1 != 0) {
    stop("N not positive whole number.")
  }
  
  args <- list(outcome_model = outcome_model,
               family = family,
               trt_var = get_treatment_name(outcome_model, trt_var),
               rho = rho,
               marginal_distns = marginal_distns,
               marginal_params = marginal_params,
               N = N)
  
  do.call(new_strategy, c(strategy = "gcomp_bayes", args))
}

#' @rdname strategy
#'
#' @param trt_var Treatment variable name; string
#' 
#' @section Multiple imputation marginalization (MIM):
#' MIM targets a marginal treatment effect by using parametric G-computation
#' within a multiple imputation framework. This approach views the covariate
#' adjustment regression as a nuisance model and separates its estimation from
#' the evaluation of the marginal treatment effect of interest. It is
#' particularly useful for ensuring compatibility in indirect comparisons
#' when adjusting for effect modifiers.
#' 
#' @return `mim` class object
#' @importFrom utils modifyList
#' @export
# 
strategy_mim <- function(formula = NULL,
                         family = gaussian(link = "identity"),
                         trt_var= NULL,
                         rho = NA,
                         marginal_distns = NA,
                         marginal_params = NA,
                         N = 1000L) {
  # back-compatibility
  if (!is.list(formula)) {
    outcome_model <- formula
  } else {
    outcome_model <- formula$outcome_model
  }
  
  check_formula(outcome_model, trt_var)
  check_family(family)
  check_distns(outcome_model, marginal_distns, marginal_params)
  check_rho(rho)
  
  if (N <= 0 || N %% 1 != 0) {
    stop("N not positive whole number.")
  }
  
  args <- list(outcome_model = outcome_model,
               family = family,
               trt_var = get_treatment_name(outcome_model, trt_var),
               rho = rho,
               marginal_distns = marginal_distns,
               marginal_params = marginal_params,
               N = N)
  
  do.call(new_strategy, c(strategy = "mim", args))
}

#' @name strategy
#' @title New strategy objects
#' 
#' @description
#' Create a type of strategy class for each modelling approach.
#'
#' @note While current implementations focus on binary, continuous, and count outcomes, 
#' support for survival data (using the \code{survival} package) is under active 
#' development and scheduled for a future version.
#' 
#' @param strategy Class name from `strategy_maic`, `strategy_stc`, `strategy_gcomp_ml`, `strategy_gcomp_bayes`, `strategy_mim`
#' @eval reg_args(include_formula = TRUE, include_family = TRUE)
#' @param ... Additional arguments
#' @returns Strategy list object
#'
#' @export
#'
new_strategy <- function(strategy, ...) {
  structure(list(...), class = c(strategy, "strategy", "list"))
}

#' @keywords internal
is_family <- function(obj) inherits(obj, "family")

#' @keywords internal
check_family <- function(obj) {
  if (!is_family(obj)) {
    stop("family must be a family object", call. = FALSE)
  }
}

#' @keywords internal
check_rho <- function(mat) {
  if (is.na(mat)) return()
  
  rn <- rownames(mat)
  cn <- colnames(mat)
  
  # names actually exist
  if (is.null(rn) || is.null(cn)) {
    stop("Validation Failed: Matrix must have both row and column names.", call. = FALSE)
  }
  
  # names are identical
  if (!identical(rn, cn)) {
    stop("Validation Failed: Row names and column names do not match.", call. = FALSE)
  }
}

#' @keywords internal
check_distns <- function(formula,
                         marginal_distns,
                         marginal_params) {
  
  if (length(marginal_distns) == 1 && is.na(marginal_distns) &&
      length(marginal_params) == 1 && is.na(marginal_params)) {
    return()
  }
  
  covariate_names <- get_covariate_names(formula)
  n_covariates <- length(covariate_names) - 1  # remove treatment
  
  if (length(marginal_distns) > n_covariates) {
    stop("Number of marginal distributions cannot be larger than
           the number of covariates in the formula.", call. = FALSE)
  }
  
  if (!all(is.na(marginal_params))) {
    if (length(marginal_params) > n_covariates) {
      stop("Number of marginal parameter lists cannot be larger than
           the number of covariates in the formula.", call. = FALSE)
    }
  }
}


##TODO:
## generic construction 
## could be useful if number of method gets big
#
# strategy_template <- function(formula = NULL,
#                                 family = gaussian(link = "identity"),
#                                 rho = NA,
#                                 N = 1000L) {
#   check_formula(formula)
#   check_family(family)
#   
#   if (N <= 0 || N %% 1 != 0) {
#     stop("N not positive whole number.")
#   }
#   
#   force(family)
#   force(formula) 
#   
#   default_args <- formals()
#   args <- c(formula = formula, as.list(match.call())[-c(1,2)])
#   args <- modifyList(default_args, args)
#   do.call(new_strategy, c(strategy = "gcomp_bayes", args))
# }



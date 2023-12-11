
# create class for each approach

#' @rdname strategy
#' 
#' @section Matching-adjusted indirect comparison (MAIC): 
#' Used to compare marginal
#' treatment effects where there are cross-trial differences in effect modifiers
#' and limited patient-level data.
#'
#' The default formula is
#' \deqn{
#'  y = X3 + X4 + trt*X1 + trt*X2
#' }
#' 
#' @param formula Linear regression formula object 
#' @param R 
#' @param ald Aggregate-level data 
#'
#' @return `maic` class object
#' @export
#'
strategy_maic <- function(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"),
                          R = 1000,
                          dat_ALD = BC.ALD) {
  default_args <- formals()
  args <- as.list(match.call())[-1]
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "maic", args))
}

#' @rdname strategy
#' 
#' @section Simulated treatment comparison (STC): 
#' Outcome regression-based method which targets a conditional treatment effect.
#' 
#' The default formula is
#' \deqn{
#'  y = X3 + X4 + trt*(X1 - mean(X1)) + trt*(X2 - mean(X2))
#' }
#' 
#' @param formula Linear regression formula object
#'
#' @return `stc` class object
#' @export
# 
strategy_stc <- function(formula =
                           as.formula("y ~ X3 + X4 +
                                   trt*I(X1 - mean(X1)) +
                                   trt*I(X2 - mean(X2))")) {
  default_args <- formals()
  args <- as.list(match.call())[-1]
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "stc", args))
}

#' @rdname strategy
#' 
#' @section G-computation maximum likelihood:
#'
#' The default formula is
#' \deqn{
#'  y = X3 + X4 + trt*X1 + trt*X2
#' }
#'
#' @param formula Linear regression formula object
#' @param R 
#' 
#' @return `gcomp_ml` class object
#' @export
#'
strategy_gcomp_ml <- function(formula =
                                as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"),
                              R = 1000) {
  default_args <- formals()
  args <- as.list(match.call())[-1]
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "gcomp_ml", args))
}

#' @rdname strategy
#' 
#' @section G-computation Bayesian:
#'
#' The default formula is
#' \deqn{
#'  y = X3 + X4 + trt*X1 + trt*X2
#' }
#' 
#' @param formula Linear regression formula object
#'
#' @return `gcomp_stan` class object
#' @export
#'
strategy_gcomp_stan <- function(formula =
                                  as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")) {
  default_args <- formals()
  args <- as.list(match.call())[-1]
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "gcomp_stan", args))
}

#' @name strategy
#' @title New strategy objects
#' 
#' @description
#' Create class for each approach
#'
#' @param strategy Class name
#' @param ... Additional arguments
#'
#' @export
#'
new_strategy <- function(strategy, ...) {
  structure(list(...), class = strategy)
}


#' @title Calculate the difference between treatments using all evidence
#' 
#' @description
#' This is the main wrapper for `hat_Delta_stats()`.
#' 
#' \insertCite{RemiroAzocar2022}{mimR}
#' 
#' @param AC.IPD Individual-level patient data. Suppose between studies A and C.
#' @param BC.ALD Aggregate-level data. Suppose between studies B and C. 
#' @param strategy Computation strategy function. These can be
#'    `strategy_maic()`, `strategy_stc()`, `strategy_gcomp_ml()` and  `strategy_gcomp_stan()`
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
#' 
hat_Delta_stats <- function(AC.IPD, BC.ALD, strategy, ...) {

  AC_hat_Delta_stats <- IPD_stats(strategy, ipd = AC.IPD, ald = BC.ALD, ...) 
  BC_hat_Delta_stats <- ALD_stats(data = BC.ALD) 
  
  ci_range <- c(0.025, 0.975)
  
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
  
  structure(stats, class = c("mimR", class(stats)))
}


#' @name IPD_stats
#' @title Individual level data statistics
#' 
#' @param strategy A list corresponding to different approaches
#' @template args-ipd
#' @template args-ald
#' @param ... Additional arguments
#' 
#' @return Mean and variance
#' @export
#' 
IPD_stats <- function(strategy, ipd, ald, ...)
  UseMethod("IPD_stats", strategy)


#' @rdname IPD_stats
#' 
IPD_stats.default <- function() {
  stop("strategy not available.")
}


#' @rdname IPD_stats
#' @title Matching-adjusted indirect comparison statistics
#' 
#' @description
#' Marginal A vs C treatment effect estimates
#' using bootstrapping
#'
#' @param strategy A list corresponding to different approaches
#' @template args-ipd
#' @template args-ald
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
                          dat_ALD = strategy$dat_ALD)   #TODO: swap ald
  
  list(mean = mean(maic_boot$t),
       var = var(maic_boot$t))
}


#' @rdname IPD_stats
#' @title Simulated treatment comparison statistics 
#' 
#' @description
#' IPD from the AC trial are used to fit a regression model describing the
#' observed outcomes `y` in terms of the relevant baseline characteristics `x` and
#' the treatment variable `z`.
#' 
#' @param strategy A list corresponding to different approaches
#' @template args-ipd
#' @template args-ald
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
#' @title G-computation maximum likelihood statistics
#'
#' @param strategy A list corresponding to different approaches
#' @template args-ipd
#' @template args-ald
#'
#' @export
#'
IPD_stats.gcomp_ml <- function(strategy,
                               ipd, ald) {

  # non-parametric bootstrap with 1000 resamples
  AC_maic_boot <- boot::boot(data = ipd,
                             statistic = gcomp_ml.boot,
                             R = strategy$R,
                             formula = strategy$formula)
  
  list(mean = mean(AC_maic_boot$t),
       var = var(AC_maic_boot$t))
}


#' @rdname IPD_stats
#' @title G-computation Bayesian statistics
#'
#' @param strategy A list corresponding to different approaches
#' @template args-ipd
#' @template args-ald
#'
#' @export
#'
IPD_stats.gcomp_stan <- function(strategy,
                                 ipd, ald) {
  
  ppv <- gcomp_stan(formula = strategy$formula,
                    ipd = ipd, ald = ald)
  
  # compute marginal log-odds ratio for A vs C for each MCMC sample
  # by transforming from probability to linear predictor scale
  hat.delta.AC <-
    qlogis(rowMeans(ppv$y.star.A)) - qlogis(rowMeans(ppv$y.star.C))
  
  list(mean = mean(hat.delta.AC),
       var = var(hat.delta.AC))
} 


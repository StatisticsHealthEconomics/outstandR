
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
                          ald = ald)
  
  list(mean = mean(maic_boot$t),
       var = var(maic_boot$t))
}


#' @rdname IPD_stats
#' @section Simulated treatment comparison statistics: 
#' IPD from the _AC_ trial are used to fit a regression model describing the
#' observed outcomes \eqn{y} in terms of the relevant baseline characteristics \eqn{x} and
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
#' by transforming from probability to linear predictor scale.
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


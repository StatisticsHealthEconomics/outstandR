
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
IPD_stats <- function(strategy, ipd, ald, scale, ...)
  UseMethod("IPD_stats", strategy)


#' @rdname IPD_stats
#' @importFrom utils methods
#' @export
IPD_stats.default <- function(...) {
  strategy_classes <- sub("IPD_stats\\.", "", methods(IPD_stats)[-1])
  avail_strategies <- paste0("strategy_", strategy_classes, "()", collapse = ", ")
  stop(paste0("strategy not available. Select from ", avail_strategies))
}


#' @rdname IPD_stats
#' 
#' @section Multiple imputation marginalisation:
#' Using Stan, compute marginal relative treatment effect for _A_ vs _C_ for each MCMC sample
#' by transforming from probability to linear predictor scale. Approximate by 
#' using imputation and combining estimates using Rubin's rules, in contrast to [IPD_stats.gcomp_stan()].
#' @import stats
#' @export
#'
IPD_stats.mim <- function(strategy,
                          ipd, ald,
                          scale, ...) {
  mis_res <-
    calc_mim(strategy,
             ipd, ald, ...)
  
  hat.delta.AC <-
    calculate_ate(mis_res$mean_A, mis_res$mean_C,
                  effect = scale)
  
  M <- mis_res$M
  
  # quantities originally defined by Rubin (1987) for multiple imputation
  coef_est <- mean(hat.delta.AC)   # average of treatment effect point estimates
  bar.v <- mean(mis_res$hats.v)    # "within" variance (average of variance point estimates)
  b <- var(hat.delta.AC)           # "between" variance (sample variance of point estimates)
  
  var_est <- var_by_pooling(M, bar.v, b)
  nu <- wald_type_interval(M, bar.v, b)
  
  ##TODO: how are these used?
  lci.Delta <- coef_est + qt(0.025, df = nu) * sqrt(var_est)
  uci.Delta <- coef_est + qt(0.975, df = nu) * sqrt(var_est)
  
  list(mean = coef_est,
       var = var_est)
} 

#' function operator
#'
IPD_stat_factory <- function(ipd_fun) {

  function(strategy, ipd, ald, scale, ...) {
    out <- ipd_fun(strategy, ipd, ald, ...)

    hat.delta.AC <- calculate_ate(out$mean_A, out$mean_C,
                                  effect = scale)

    coef_est <- mean(hat.delta.AC)
    var_est <- var(hat.delta.AC)

    list(mean = coef_est,
         var = var_est)
  }
}

#' @rdname IPD_stats
#' @section Simulated treatment comparison statistics:
#' IPD from the _AC_ trial are used to fit a regression model describing the
#' observed outcomes \eqn{y} in terms of the relevant baseline characteristics \eqn{x} and
#' the treatment variable \eqn{z}.
#' @export
#'
IPD_stats.stc <- IPD_stat_factory(outstandR:::calc_stc)

#' @rdname IPD_stats
#' @section Matching-adjusted indirect comparison statistics:
#' Marginal _A_ vs _C_ treatment effect estimates
#' using bootstrapping sampling.
#' @export
#'
IPD_stats.maic <- IPD_stat_factory(outstandR:::calc_maic)

#' @rdname IPD_stats
#' @section G-computation maximum likelihood statistics:
#' Compute a non-parametric bootstrap with default \eqn{R=1000} resamples.
#' @export
#'
IPD_stats.gcomp_ml <- IPD_stat_factory(outstandR:::calc_gcomp_ml)

#' @rdname IPD_stats
#' @section G-computation Bayesian statistics:
#' Using Stan, compute marginal log-odds ratio for _A_ vs _C_ for each MCMC sample
#' by transforming from probability to linear predictor scale.
#' @export
#'
IPD_stats.gcomp_stan <- IPD_stat_factory(outstandR:::calc_gcomp_stan)

# #' @export
#' IPD_stats.mim <- IPD_stat_factory(outstandR:::calc_mim)


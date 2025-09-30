
#' @name calc_IPD_stats
#' @title Calculate individual-level patient data statistics
#' 
#' @description
#' Computes mean and variance statistics for individual-level patient data using various approaches,
#' including Matching-Adjusted Indirect Comparison (MAIC), Simulated Treatment Comparison (STC),
#' and G-computation via Maximum Likelihood Estimation (MLE) or Bayesian inference.
#' 
#' @param strategy A list corresponding to different modelling approaches
#' @param analysis_params A list containing: 
#'   - `ald` Aggregate-level trial data
#'   - `ref_trt` Treatment labels reference (common; e.g. placebo)
#'   - `comp_trt` Treatment labels comparator
#'   - `scale` A scaling parameter for the calculation. From "log_odds", "risk_difference", "log_relative_risk".
#' @param ... Additional arguments
#' 
#' @return A list containing:
#' \describe{
#'   \item{mean}{Estimated mean treatment effect.}
#'   \item{var}{Estimated variance of the treatment effect.}
#' }
#' @examples
#' \dontrun{
#' strategy <- strategy_maic()
#' ipd <- data.frame(trt = sample(c("A", "C"), 100, replace = TRUE),
#'                   X1 = rnorm(100, 1, 1),
#'                   y = rnorm(100, 10, 2))
#' ald <- data.frame(trt = c(NA, "B", "C", "B", "C"),
#'                   variable = c("X1", "y", "y", NA, NA),
#'                   statistic = c("mean", "sum", "sum", "N", "N"),
#'                   value = c(0.5, 10, 12, 20, 25))
#' calc_IPD_stats(strategy, ipd, ald, scale = "log_odds")
#' }
#' @export
#' 
calc_IPD_stats <- function(strategy, analysis_params, ...)
  UseMethod("calc_IPD_stats", strategy)


#' @rdname calc_IPD_stats
#' @importFrom utils methods
#' @export
calc_IPD_stats.default <- function(...) {
  strategy_classes <- sub("calc_IPD_stats\\.", "", methods(calc_IPD_stats)[-1])
  avail_strategies <- paste0("strategy_", strategy_classes, "()", collapse = ", ")
  stop(paste0("strategy not available. Select from ", avail_strategies))
}


#' @param strategy Strategy
#' @param analysis_params Analysis parameters
#' @param ... Additional arguments
#'
#' @rdname calc_IPD_stats
#' 
#' @section Multiple imputation marginalisation:
#' Using Stan, compute marginal relative treatment effect for IPD
#' comparator "A" vs reference "C" arms for each MCMC sample
#' by transforming from probability to linear predictor scale. Approximate by 
#' using imputation and combining estimates using Rubin's rules,
#' in contrast to [calc_IPD_stats.gcomp_stan()].
#' 
#' @importFrom stats qt var
#' @export
#'
calc_IPD_stats.mim <- function(strategy,
                               analysis_params, ...) {
  
  ipd <- analysis_params$ipd
  ald <- analysis_params$ald
  scale <- analysis_params$scale
  ref_trt <- analysis_params$ref_trt
  comp_trt <- analysis_params$ipd_comp
  
  mis_res <-
    calc_mim(strategy,
             ipd, ald, 
             ref_trt, comp_trt, 
             ...)
  
  hat.delta.AC <-
    calculate_ate(mis_res$mean_comp, mis_res$mean_ref,
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
  
  list(
    contrasts = list(
      mean = coef_est,
      var = var_est),
    absolute = list(
      mean = NA,  #p_est,  ##TODO:
      var = NA)   #p_var)
  )
} 

#' Factory function for creating calc_IPD_stats methods
#'
#' Creates a method for computing IPD mean and variance statistics based on the supplied function.
#'
#' @param ipd_fun A function that computes mean and variance statistics for individual-level patient data.
#' @return A function that computes mean and variance statistics for a given strategy.
#' @keywords internal
#'
IPD_stat_factory <- function(ipd_fun) {
  
  function(strategy, analysis_params,
           var_method = "sample", ...) {
    
    ipd <- analysis_params$ipd
    ald <- analysis_params$ald
    scale <- analysis_params$scale
    
    out <- ipd_fun(strategy, analysis_params, ...)
    
    # relative treatment effect
    hat.delta.AC <- calculate_ate(out$mean_A, out$mean_C,
                                  effect = scale)
    
    coef_est <- mean(hat.delta.AC, na.rm = TRUE)
    
    if (var_method == "sandwich") {
      ##TODO:
      var_est <- estimate_var_sandwich(strategy, ipd, ...)
    } else if (var_method == "sample") {
      var_est <- var(hat.delta.AC, na.rm = TRUE)
    }
    
    p_est <- sapply(out, mean, na.rm = TRUE)
    p_var <- sapply(out, var, na.rm = TRUE)
    
    list(
      contrasts = list(
        mean = coef_est,
        var = var_est),
      absolute = list(
        mean = p_est,
        var = p_var)
    )
  }
}

#' @rdname calc_IPD_stats
#' @section Simulated treatment comparison statistics:
#' IPD for reference "C" and comparator "A" trial arms are used to fit a regression model describing the
#' observed outcomes \eqn{y} in terms of the relevant baseline characteristics \eqn{x} and
#' the treatment variable \eqn{z}.
#' @export
#'
calc_IPD_stats.stc <- IPD_stat_factory(outstandR:::calc_stc)

#' @rdname calc_IPD_stats
#' @section Matching-adjusted indirect comparison statistics:
#' Marginal IPD comparator treatment "A" vs reference treatment "C" treatment effect estimates
#' using bootstrapping sampling.
#' @export
#'
calc_IPD_stats.maic <- IPD_stat_factory(outstandR:::calc_maic)

#' @rdname calc_IPD_stats
#' @section G-computation maximum likelihood statistics:
#' Compute a non-parametric bootstrap with default \eqn{R=1000} resamples.
#' @export
#'
calc_IPD_stats.gcomp_ml <- IPD_stat_factory(outstandR:::calc_gcomp_ml)

#' @rdname calc_IPD_stats
#' @section G-computation Bayesian statistics:
#' Using Stan, compute marginal relative effects for IPD comparator "A" vs reference "C" treatment arms for each MCMC sample
#' by transforming from probability to linear predictor scale.
#' @export
#'
calc_IPD_stats.gcomp_stan <- IPD_stat_factory(outstandR:::calc_gcomp_stan)

# #' @export
#' calc_IPD_stats.mim <- IPD_stat_factory(outstandR:::calc_mim)


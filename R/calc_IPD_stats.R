
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
#'   - `ipd`: Individual-level patient data (data frame)
#'   - `ald`: Aggregate-level trial data (data frame)
#'   - `ref_trt`: Treatment label for the reference arm (common; e.g., "C")
#'   - `ipd_comp`: Treatment label for the comparator arm in the IPD (e.g., "A")
#'   - `scale`: Scaling parameter ("log_odds", "risk_difference", "log_relative_risk")
#' @param ... Additional arguments
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{contrasts}: A list with elements \code{mean} and \code{var}.
#'   \item \code{absolute}: A list with elements \code{mean} and \code{var}.
#' }
#' @examples
#' strategy <- strategy_maic(formula = as.formula(y~trt:X1), family = binomial())
# 
#' ipd <- data.frame(trt = sample(c("A", "C"), size = 100, replace = TRUE),
#'                   X1 = rnorm(100, 1, 1),
#'                   y = sample(c(1,0), size = 100, prob = c(0.7,0.3), replace = TRUE))
#' 
#' ald <- data.frame(trt = c(NA, "B", "C", "B", "C"),
#'                   variable = c("X1", "y", "y", NA, NA),
#'                   statistic = c("mean", "sum", "sum", "N", "N"),
#'                   value = c(0.5, 10, 12, 20, 25))
#' 
#' calc_IPD_stats(strategy,
#'   analysis_params = list(ipd = ipd, ald = ald, scale = "log_odds"))
#'   
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
  stop(paste0("strategy not available. Select from ", avail_strategies), call. = FALSE)
}


#' @param strategy An object of class `strategy` created by functions such as 
#'   [strategy_maic()], [strategy_stc()], or [strategy_mim()]. 
#'   Contains modelling details like the formula and family.
#' @param analysis_params Analysis parameters
#' @param ... Additional arguments
#'
#' @rdname calc_IPD_stats
#' 
#' @section Multiple imputation marginalisation:
#' Using Stan, compute marginal relative treatment effect for IPD
#' comparator "A" vs reference "C" arms for each MCMC sample
#' by transforming from probability to linear predictor scale. Approximate by 
#' using imputation and combining estimates using Rubin's rules.
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
    calculate_ate(mis_res$means$A, mis_res$means$C,
                  effect = scale)
  
  # posterior predictive sample size
  M <- mis_res$model$M
  
  # quantities originally defined by Rubin (1987) for multiple imputation
  coef_est <- mean(hat.delta.AC)   # average of treatment effect point estimates
  bar.v <- mean(mis_res$model$hats.v)    # "within" variance (average of variance point estimates)
  b <- var(hat.delta.AC)           # "between" variance (sample variance of point estimates)
  
  var_est <- var_by_pooling(M, bar.v, b)
  nu <- wald_type_interval(M, bar.v, b)

  p_est <- sapply(mis_res$means, mean, na.rm = TRUE)
  p_var <- sapply(mis_res$means, var, na.rm = TRUE)
  
  list(
    contrasts = list(
      mean = coef_est,
      var  = var_est),
    absolute = list(
      mean = p_est,
      var  = p_var),
    model = c(method_name = "MIM",
              nu = nu, 
              mis_res$model)
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
  
  # capture name of the function passed in
  method_raw <- deparse(substitute(ipd_fun))
  method_name <- method_raw |> 
    gsub(pattern = ".*[:]", replacement = "") |> 
    gsub(pattern = "calc_", replacement = "") |> 
    toupper()
  
  function(strategy, analysis_params,
           var_method = "sample", ...) {
    
    ipd <- analysis_params$ipd
    ald <- analysis_params$ald
    scale <- analysis_params$scale
    
    var_method <- analysis_params$var_method
    if (is.null(var_method)) var_method <- "sample"
    
    out <- ipd_fun(strategy, analysis_params, ...)
    
    mean_comp <- out$means$A
    mean_ref <- out$means$C
    
    # relative treatment effect
    hat.delta.AC <- calculate_ate(mean_comp, mean_ref,
                                  effect = scale)
    
    coef_est <- mean(hat.delta.AC, na.rm = TRUE)
    
    if (var_method == "sandwich") {
      var_est <- estimate_var_sandwich(strategy, analysis_params, ...)
    } else if (var_method == "sample") {
      var_est <- var(hat.delta.AC, na.rm = TRUE)
    } else {
      stop("Variance method not known.", call. = FALSE)
    }
    
    p_est <- sapply(out$means, mean, na.rm = TRUE)
    p_var <- sapply(out$means, var, na.rm = TRUE)
    
    list(
      contrasts = list(
        mean = coef_est,
        var  = var_est),
      absolute = list(
        mean = p_est,
        var  = p_var),
      method_name = method_name,
      model = out$model
    )
  }
}

#' @rdname calc_IPD_stats
#' @param var_method A string specifying the variance estimation method,
#'   either "sample" (default) or "sandwich".
#' @section Simulated treatment comparison statistics:
#' IPD for reference "C" and comparator "A" trial arms are used to fit a regression model describing the
#' observed outcomes \eqn{y} in terms of the relevant baseline characteristics \eqn{x} and
#' the treatment variable \eqn{z}.
#' @export
#'
calc_IPD_stats.stc <- IPD_stat_factory(outstandR:::calc_stc)

#' @rdname calc_IPD_stats
#' @param var_method A string specifying the variance estimation method,
#'   either "sample" (default) or "sandwich".
#' @section Matching-adjusted indirect comparison statistics:
#' Marginal IPD comparator treatment "A" vs reference treatment "C" treatment effect estimates
#' using bootstrapping sampling.
#' @export
#'
calc_IPD_stats.maic <- IPD_stat_factory(outstandR:::calc_maic)

#' @rdname calc_IPD_stats
#' @param var_method A string specifying the variance estimation method,
#'   either "sample" (default) or "sandwich".
#' @section G-computation maximum likelihood statistics:
#' Compute a non-parametric bootstrap with default \eqn{R=1000} resamples.
#' @export
#'
calc_IPD_stats.gcomp_ml <- IPD_stat_factory(outstandR:::calc_gcomp_ml)

#' @rdname calc_IPD_stats
#' @param var_method A string specifying the variance estimation method,
#'   either "sample" (default) or "sandwich".
#' @section G-computation Bayesian statistics:
#' Using Stan, compute marginal relative effects for IPD comparator "A" vs reference "C" treatment arms for each MCMC sample
#' by transforming from probability to linear predictor scale.
#' @export
#'
calc_IPD_stats.gcomp_bayes <- IPD_stat_factory(outstandR:::calc_gcomp_bayes)

# #' @export
#' calc_IPD_stats.mim <- IPD_stat_factory(outstandR:::calc_mim)


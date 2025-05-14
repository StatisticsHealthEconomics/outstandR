
#' @title Calculate the difference between treatments using all evidence
#' 
#' @description
#' This is the main, top-level wrapper for `{outstandR}`.
#' Methods taken from
#' \insertCite{RemiroAzocar2022}{outstandR}.
#' 
#' @param ipd_trial Individual-level patient data. Suppose between studies _A_ and _C_.
#' @param ald_trial Aggregate-level data. Suppose between studies _B_ and _C_. 
#' @param strategy Computation strategy function. These can be
#'    `strategy_maic()`, `strategy_stc()`, `strategy_gcomp_ml()` and  `strategy_gcomp_stan()`
#' @param CI Confidence interval; between 0,1
#' @param scale Relative treatment effect scale. If `NULL`, the scale is automatically determined from the model.
#' @param ... Additional arguments
#' @return List of length 3 of statistics as a `outstandR` class object.
#'   Containing statistics between each pair of treatments.
#'   These are the mean, variances and confidence intervals,
#'   for contrasts and absolute values.
#' @importFrom Rdpack reprompt
#' 
#' @references
#' \insertRef{RemiroAzocar2022}{outstandR}
#' 
#' @export
#' @examples
#' data(AC_IPD)  # AC patient-level data
#' data(BC_ALD)  # BC aggregate-level data
#' 
#' lin_form <- as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")
#' 
#' # matching-adjusted indirect comparison
#' outstandR_maic <- outstandR(AC_IPD, BC_ALD, strategy = strategy_maic(formula = lin_form))
#' 
#' # simulated treatment comparison
#' outstandR_stc <- outstandR(AC_IPD, BC_ALD, strategy = strategy_stc(lin_form))
#' 
#' # G-computation with maximum likelihood
#' # outstandR_gcomp_ml <- outstandR(AC_IPD, BC_ALD, strategy = strategy_gcomp_ml(lin_form))
#' 
#' # G-computation with Bayesian inference
#' outstandR_gcomp_stan <- outstandR(AC_IPD, BC_ALD, strategy = strategy_gcomp_stan(lin_form))
#' 
#' # Multiple imputation marginalization
#' outstandR_gcomp_stan <- outstandR(AC_IPD, BC_ALD, strategy = strategy_mim(lin_form))
#' 
outstandR <- function(ipd_trial, ald_trial, ref_trt = "C",
                      strategy, CI = 0.95, scale = NULL, ...) {
  
  validate_outstandr(ipd_trial, ald_trial, strategy, CI, scale)
  
  ipd <- prep_ipd(strategy$formula, ipd_trial)
  ald <- prep_ald(strategy$formula, ald_trial, trt_var = strategy$trt_var)

  ipd_trts <- get_ipd_trts(ipd, ref_trt)
  ald_trts <- get_ald_trts(ald, ref_trt)
  
  if (is.null(scale)) scale <- get_treatment_effect(strategy$family$link)
  
  ipd_stats <- calc_IPD_stats(strategy, ipd = ipd, ald = ald, scale, ...) 
  ald_stats <- calc_ALD_stats(strategy, ald = ald, treatments = ald_trts, scale = scale) 
  
  stats <- result_stats(ipd_stats, ald_stats, CI)
  
  structure(stats,
            CI = CI,
            scale = scale,
            model = strategy$family$family,
            class = c("outstandR", class(stats)))
}

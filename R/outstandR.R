
#' @title Calculate the difference between treatments using all evidence
#' 
#' @description
#' This is the main, top-level wrapper for `{outstandR}`.
#' Methods taken from
#' \insertCite{RemiroAzocar2022}{outstandR}.
#' 
#' @param AC.IPD Individual-level patient data. Suppose between studies _A_ and _C_.
#' @param BC.ALD Aggregate-level data. Suppose between studies _B_ and _C_. 
#' @param strategy Computation strategy function. These can be
#'    `strategy_maic()`, `strategy_stc()`, `strategy_gcomp_ml()` and  `strategy_gcomp_stan()`
#' @param CI Confidence interval; between 0,1
#' @param scale Relative treatment effect scale. If `NULL`, the scale is automatically determined from the model.
#' @param ... Additional arguments
#' @return List of length 3 of statistics as a `outstandR` class object.
#'   Containing statistics between each pair of treatments.
#'   These are the mean contrasts, variances and confidence intervals,
#'   respectively.
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
outstandR <- function(AC.IPD, BC.ALD, strategy, CI = 0.95, scale = NULL, ...) {
  
  if (CI <= 0 || CI >= 1) stop("CI argument must be between 0 and 1.")
  ##TODO: as method instead?
  if (!inherits(strategy, "strategy"))
    stop("strategy argument must be a class strategy.")
  
  ipd <- prep_ipd(strategy$formula, AC.IPD)
  ald <- prep_ald(strategy$formula, BC.ALD)

  if (is.null(scale)) scale <- get_treatment_effect(strategy$family$link)
  
  AC_stats <- IPD_stats(strategy, ipd = ipd, ald = ald, scale, ...) 
  BC_stats <- ALD_stats(strategy, ald = ald, scale = scale) 
  
  stats <- contrast_stats(AC_stats, BC_stats, CI)
  
  structure(stats,
            CI = CI,
            class = c("outstandR", class(stats)))
}

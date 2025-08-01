
#' @title Calculate the difference between treatments using all evidence
#' 
#' @description
#' This is the main, top-level wrapper for `{outstandR}`.
#' Methods taken from
#' \insertCite{RemiroAzocar2022}{outstandR}.
#' 
#' @param ipd_trial Individual-level patient data. For example, suppose between studies _A_ and _C_.
#'   In a long format and must contain a treatment column and outcome column consistent with the formula object.
#'   The labels in the treatment are used internally so there must be a common treatment with the aggregate-level data trial.
#' @param ald_trial Aggregate-level data. For example, suppose between studies _B_ and _C_. The column names are
#'  - `variable`: Covariate name. In the case of treatment arm sample size this is `NA`
#'  - `statistic`: Summary statistic name from "mean", standard deviation "sd", probability "prop", or "sum"
#'  - `value`: Numerical value of summary statistic
#'  - `trt`: Treatment label. Because we assume a common covariate distribution between treatment arms this is `NA`
#' @param strategy Computation strategy function. These can be
#'    `strategy_maic()`, `strategy_stc()`, `strategy_gcomp_ml()` and `strategy_gcomp_stan()`.
#' @param ref_trt Reference / common / anchoring treatment name.
#' @param CI Confidence interval; between 0,1
#' @param scale Relative treatment effect scale. If `NULL`, the scale is automatically determined from the model.
#'   Choose from "log-odds", "log_relative_risk", "risk_difference", "delta_z", "mean_difference", "rate_difference" depending on the data type.
#' @param ... Additional arguments. Currently, can pass named arguments to `rstanarm::stan_glm()` via `strategy_gcomp_stan()`.
#' 
#' @return List of length 3 of statistics as a `outstandR` class object.
#'   Containing statistics between each pair of treatments.
#'   These are the mean, variances and confidence intervals,
#'   for contrasts and absolute values.
#' @importFrom Rdpack reprompt
#' @seealso [strategy_maic()] [strategy_stc()] [strategy_gcomp_ml()] [strategy_gcomp_stan()]
#' 
#' @references
#' \insertRef{RemiroAzocar2022}{outstandR}
#' 
#' @export
#' @examples
#' data(AC_IPD)  # AC patient-level data
#' data(BC_ALD)  # BC aggregate-level data
#' 
#' # linear formula
#' lin_form <- as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")
#' 
#' # matching-adjusted indirect comparison
#' outstandR_maic <- outstandR(AC_IPD, BC_ALD,
#'                             strategy = strategy_maic(formula = lin_form))
#' 
#' # simulated treatment comparison
#' outstandR_stc <- outstandR(AC_IPD, BC_ALD,
#'                            strategy = strategy_stc(lin_form))
#' 
#' # G-computation with maximum likelihood
#' outstandR_gcomp_ml <- outstandR(AC_IPD, BC_ALD,
#'                                 strategy = strategy_gcomp_ml(lin_form))
#' 
#' # G-computation with Bayesian inference
#' outstandR_gcomp_stan <- outstandR(AC_IPD, BC_ALD,
#'                                   strategy = strategy_gcomp_stan(lin_form))
#' 
#' # Multiple imputation marginalization
#' outstandR_mim <- outstandR(AC_IPD, BC_ALD,
#'                            strategy = strategy_mim(lin_form))
#' 
outstandR <- function(ipd_trial, ald_trial, strategy,
                      ref_trt = NA,
                      CI = 0.95, scale = NULL, ...) {
  
  validate_outstandr(ipd_trial, ald_trial, strategy, CI, scale)

  trt_var <- strategy$trt_var
  
  ipd <- prep_ipd(strategy$formula, ipd_trial)
  ald <- prep_ald(strategy$formula, ald_trial, trt_var = trt_var)

  ref_trt <- get_ref_trt(ref_trt, trt_var, ipd_trial, ald_trial)
  
  # treatment names for each study
  ipd_comp <- get_comparator(ipd, ref_trt, trt_var)
  ald_comp <- get_comparator(ald, ref_trt, trt_var)
  
  if (is.null(scale)) scale <- get_treatment_effect(strategy$family$link)
  
  analysis_params <- list(
    ipd = ipd, 
    ald = ald,
    scale = scale,
    trt_var = trt_var,
    ref_trt = ref_trt,
    ipd_comp = ipd_comp,
    ald_comp = ald_comp
  )
  
  ipd_stats <- calc_IPD_stats(strategy, analysis_params, ...) 
  ald_stats <- calc_ALD_stats(strategy, analysis_params) 
  
  stats <- result_stats(ipd_stats, ald_stats, CI)
  
  structure(
    stats,
    CI = CI,
    ref_trt = ref_trt,
    scale = scale,
    model = strategy$family$family,
    class = c("outstandR", class(stats)))
}

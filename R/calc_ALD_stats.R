
#' Aggregate-level data mean and variance statistics
#'
#' Computes the mean and variance of marginal treatment effects for aggregate-level trial data.
#'
#' @param strategy A list containing the strategy details, including the family distribution.
#' @param analysis_params A list containing: 
#'   - `ald` Aggregate-level trial data
#'   - `ref_trt` Treatment labels reference (common; e.g. placebo)
#'   - `comp_trt` Treatment labels comparator
#'   - `scale` A scaling parameter for the calculation. From "log_odds", "risk_difference", "log_relative_risk".
#'
#' @return A list containing:
#' \describe{
#'   \item{`mean`}{The marginal treatment effect mean.}
#'   \item{`var`}{The marginal treatment effect variance.}
#' }
#' @seealso [marginal_treatment_effect()], [marginal_variance()]
#' @export
#' @examples
#' \dontrun{
#' strategy <- list(family = list(family = "binomial"))  # basic version
#' ald <- data.frame(trt = c("B","C","B","C"),
#'                   variable = c(NA, NA, "y", "y"),
#'                   statistic = c("N", "N", "sum", "sum"),
#'                   value = c(100, 100, 50, 60) 
#'
#' calc_ALD_stats(strategy = strategy,
#'                list(ald = ald,
#'                     ref_trt = "C",
#'                     comp_trt = "B",
#'                     scale = "log"))
#' }
#'
calc_ALD_stats <- function(strategy, analysis_params) {
  
  ald <- analysis_params$ald
  ref_trt <- analysis_params$ref_trt
  comp_trt <- analysis_params$ald_comp
  scale <- analysis_params$scale
  
  family <- strategy$family$family
  
  ald_cc <- continuity_correction(ald)
  
  mean_eff <- marginal_treatment_effect(ald_cc, ref_trt, comp_trt, scale, family)
  var_eff <- marginal_variance(ald_cc, ref_trt, comp_trt, scale, family)
  
  list(mean = mean_eff,
       var = var_eff)
}


#' Marginal effect variance using the delta method
#' 
#' Computes the total variance of marginal treatment effects using the delta method.
#'
#' @param ald Aggregate-level data
#' @param ref_trt Treatment labels reference (common; e.g. placebo)
#' @param comp_trt Treatment labels comparator
#' @param scale A scaling parameter for the calculation.
#' @param family A character string specifying the family distribution (e.g., "binomial").
#' 
#' @return The total variance of marginal treatment effects.
#' @examples
#' \dontrun{
#' ald <- data.frame(trt = c("B","C","B","C"),
#'                   variable = c(NA, NA, "y", "y"),
#'                   statistic = c("N", "N", "sum", "sum"),
#'                   value = c(100, 100, 50, 60)
#'                   
#' marginal_variance(ald, ref_trt = "C", comp_trt = "B",
#'                   scale = "log", family = "binomial")
#' }
#' @export
#' 
marginal_variance <- function(ald,
                              ref_trt = NA,
                              comp_trt = NA,
                              scale,
                              family) {
  
  v_ref <- calculate_trial_variance(ald, ref_trt, scale, family)
  v_comp <- calculate_trial_variance(ald, comp_trt, scale, family)
  
  v_comp + v_ref
}


#' Marginal treatment effect from reported event counts
#' 
#' Computes the relative treatment effect from aggregate-level data using event counts.
#'
#' @param ald Aggregate-level data
#' @param ref_trt Treatment labels reference (common; e.g. placebo)
#' @param comp_trt Treatment labels comparator
#' @param scale A scaling parameter for the calculation.
#' @param family A character string specifying the family distribution (e.g., "binomial").
#' 
#' @return The relative treatment effect.
#' @examples
#' \dontrun{
#' ald <- data.frame(trt = c("B","C","B","C"),
#'                   variable = c(NA, NA, "y", "y"),
#'                   statistic = c("N", "N", "sum", "sum"),
#'                   value = c(100, 100, 50, 60)
#'                   
#' marginal_treatment_effect(ald, ref_trt = "C", comp_trt = "B",
#'                           scale = "log", family = "binomial")
#' }
#' @export
#' 
marginal_treatment_effect <- function(ald,
                                      ref_trt = NA,
                                      comp_trt = NA,
                                      scale,
                                      family) {
  
  m_ref <- calculate_trial_mean(ald, ref_trt, scale, family)
  m_comp <- calculate_trial_mean(ald, comp_trt, scale, family)
  
  m_comp - m_ref
}


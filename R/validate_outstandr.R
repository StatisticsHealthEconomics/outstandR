
#' Input data validator
#' 
#' @param ipd_trial Individual patient data
#' @param ald_trial Aggregate level data
#' @param strategy Strategy object
#' @param CI Confidence interval
#' @param scale Outcome scale
#' @returns No return value, called for side effects
#'
#' @keywords internal
#' 
validate_outstandr <- function(ipd_trial, ald_trial,
                               strategy,
                               CI, scale) {
  
  if (CI <= 0 || CI >= 1) {
    stop("CI argument must be between 0 and 1.")
  }
  
  ##TODO: link this to actual functions
  available_scales <-
    c("log_odds", "log_relative_risk", "risk_difference", "mean_difference")
  
  if (!is.null(scale) && !any(scale %in% available_scales)) {
    stop("scale not in available list.")
  }
  
  if (!inherits(strategy, "strategy")) {
    stop("strategy argument must be a class strategy.")
  }
 
  return()  
}


#' Input data validator
#' @keywords internal
#' 
validate_outstandr <- function(ipd_trial, ald_trial,
                               strategy,
                               CI, scale) {
  
  if (CI <= 0 || CI >= 1) {
    stop("CI argument must be between 0 and 1.")
  }
  
  if (!is.null(scale) && !any(scale %in% c("log_odds", "log_relative_risk", "risk_difference",
                                           "mean_difference"))) {
    stop("scale not in available list.")
  }
  
  if (!inherits(strategy, "strategy")) {
    stop("strategy argument must be a class strategy.")
  }
 
  return()  
}

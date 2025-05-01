
#' Aggregate-level data mean and variance statistics
#'
#' Computes the mean and variance of marginal treatment effects for aggregate-level trial data.
#'
#' @param strategy A list containing the strategy details, including the family distribution.
#' @param ald Aggregate-level trial data
#' @param treatments Treatment labels list; default `B`, `C` (common; e.g. placebo)
#' @param scale A scaling parameter for the calculation.
#'
#' @return A list containing:
#' \describe{
#'   \item{mean}{The marginal treatment effect mean.}
#'   \item{var}{The marginal treatment effect variance.}
#' }
#' @seealso [marginal_treatment_effect()], [marginal_variance()]
#' @export
#' @examples
#' \dontrun{
#' strategy <- list(family = list(family = "binomial"))  # basic version
#' ald <- data.frame(trial = 1:5, n_B = c(10, 20, 15, 30, 25), n_C = c(12, 18, 20, 25, 22))
#' ALD_stats(strategy, ald, treatments = list("B", "C"), scale = "log")
#' }
#'
ALD_stats <- function(strategy,
                      ald,
                      treatments = list("B", "C"),
                      scale) {
  family <- strategy$family$family
  
  ald_cc <- continuity_correction(ald, treatments)
  
  mean_eff <- marginal_treatment_effect(ald_cc, treatments, scale, family)
  var_eff <- marginal_variance(ald_cc, treatments, scale, family)
  
  list(mean = mean_eff,
       var = var_eff)
}


#' Marginal effect variance using the delta method
#' 
#' Computes the total variance of marginal treatment effects using the delta method.
#' For binomial data, calculates:
#' \deqn{\frac{1}{n_C} + \frac{1}{n_{\bar{C}}} + \frac{1}{n_B} + \frac{1}{n_{\bar{B}}}}.
#'
#' @param ald Aggregate-level data
#' @param treatments A list of treatment labels; default _B_ vs _C_
#' @param scale A scaling parameter for the calculation.
#' @param family A character string specifying the family distribution (e.g., "binomial").
#' 
#' @return The total variance of marginal treatment effects.
#' @examples
#' \dontrun{
#' ald <- data.frame(trial = 1:5, n_B = c(10, 20, 15, 30, 25), n_C = c(12, 18, 20, 25, 22))
#' marginal_variance(ald, treatments = list("B", "C"), scale = "log", family = "binomial")
#' }
#' @export
#' 
marginal_variance <- function(ald,
                              treatments = list("B", "C"),
                              scale,
                              family) {
  v1 <- calculate_trial_variance(ald, treatments[[1]], scale, family)
  v2 <- calculate_trial_variance(ald, treatments[[2]], scale, family)
  
  v1 + v2
}


#' Marginal treatment effect from reported event counts
#' 
#' Computes the relative treatment effect from aggregate-level data using event counts.
#' For binomial data, calculates:
#' \deqn{
#' \log\left( \frac{n_B/(N_B-n_B)}{n_C/(N_B-n_{B})} \right) = \log(n_B n_{\bar{C}}) - \log(n_C n_{\bar{B}})
#' }
#' where \eqn{\bar{C}} is the compliment of \eqn{C}
#' so e.g. \eqn{n_{\bar{C}} = N_C - n_c}.
#'
#' @param ald Aggregate-level data
#' @param treatments A list of treatment labels. Last variable is reference; default `B`, `C` (common; e.g. placebo)
#' @param scale A scaling parameter for the calculation.
#' @param family A character string specifying the family distribution (e.g., "binomial").
#' 
#' @return The relative treatment effect.
#' @examples
#' \dontrun{
#' ald <- data.frame(trial = 1:5, n_B = c(10, 20, 15, 30, 25), n_C = c(12, 18, 20, 25, 22))
#' marginal_treatment_effect(ald, treatments = list("B", "C"), scale = "log", family = "binomial")
#' }
#' @export
#' 
marginal_treatment_effect <- function(ald,
                                      treatments = list("B", "C"),
                                      scale,
                                      family) {
  m1 <- calculate_trial_mean(ald, treatments[[1]], scale, family)
  m2 <- calculate_trial_mean(ald, treatments[[2]], scale, family)
  
  m1 - m2
}


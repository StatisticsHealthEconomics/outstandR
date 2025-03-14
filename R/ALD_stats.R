
#' Aggregate-level data mean and variance statistics
#'
#' @param ald Aggregate-level trial data
#' @param treatments Treatment labels list; default `B`, `C`
#'
#' @return List of marginal treatment effect mean and variance
#' @seealso [marginal_treatment_effect()], [marginal_variance()]
#' @export
#'
ALD_stats <- function(strategy, ald, treatments = list("B", "C"), scale) {
  list(mean = marginal_treatment_effect(ald, treatments, scale),
       var = marginal_variance(ald, treatments, scale))
}


#' Marginal effect variance using the delta method
#' 
#' For binomial data calculate
#' \deqn{\frac{1}{n_C} + \frac{1}{n_{\bar{C}}} + \frac{1}{n_B} + \frac{1}{n_{\bar{B}}}}.
#'
#' @param ald Aggregate-level data
#' @param treatments Treatment labels list; default _B_ vs _C_
#' @return Total variance
#' @export
#' 
marginal_variance <- function(ald, treatments = list("B", "C"), scale) {
  trial_vars <- purrr::map_dbl(treatments, ~calculate_trial_variance(ald, .x, scale))
  sum(trial_vars)
}


#' Marginal treatment effect from reported event counts
#' 
#' For binomial data calculate
#' \deqn{
#' \log\left( \frac{n_B/(N_B-n_B)}{n_C/(N_B-n_{B})} \right) = \log(n_B n_{\bar{C}}) - \log(n_C n_{\bar{B}})
#' }
#' where \eqn{\bar{C}} is the compliment of \eqn{C}
#' so e.g. \eqn{n_{\bar{C}} = N_C - n_c}.
#'
#' @param ald Aggregate-level data
#' @param treatments Treatment labels list. Last variable is reference; default _B_ vs _C_
#' @return Relative treatment effect
#' @export
#' 
marginal_treatment_effect <- function(ald, treatments = list("B", "C"), scale) {
  trial_means <- purrr::map_dbl(treatments, ~calculate_trial_mean(ald, .x, scale))
  trial_means[1] - trial_means[2]
}


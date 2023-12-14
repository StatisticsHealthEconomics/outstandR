
#' Marginal effect variance using the delta method
#' 
#' Calculate
#' \deqn{\frac{1}{n_C} + \frac{1}{n_{\bar{C}}} + \frac{1}{n_B} + \frac{1}{n_{\bar{B}}}}.
#'
#' @param ald Aggregate-level data
#' @param treatments Treatment labels; default _B_ vs _C_
#' @return Sum of variances
#' @export
#' 
marginal_variance <- function(ald, treatments = list("B", "C")) {
  trial_vars <- purrr::map_dbl(treatments, ~trial_variance(ald, .x))
  sum(trial_vars)
}


#' Marginal treatment effect from reported event counts
#' 
#' Calculate
#' \deqn{
#' \log\left( \frac{n_B/(N_B-n_B)}{n_C/(N_B-n_{B})} \right) = \log(n_B n_{\bar{C}}) - log(n_C n_{\bar{B}})
#' }
#' where \eqn{\bar{C}} is the compliment of \eqn{C}
#' so e.g. \eqn{n_{\bar{C}} = N_C - n_c}.
#'
#' @param ald Aggregate-level data
#' @param treatments Treatment labels; default _B_ vs _C_
#' @return Trial effect difference
#' @export
#' 
marginal_treatment_effect <- function(ald, treatments = list("B", "C")) {
  trial_effect <- purrr::map_dbl(treatments, ~trial_treatment_effect(ald, .x))
  trial_effect[2] - trial_effect[1]
}


#' Trial variance with aggregate data
#'
#' Calculate
#' \deqn{1/(\sum y_k) + 1/(N_k - \sum y_k)}.
#' 
#' @param ald Aggregate-level data
#' @param tid Treatment label
#'
#' @return Value
#' @export
#'
trial_variance <- function(ald, tid) {
  var_string <- glue::glue("1/ald$y.{tid}.sum + 1/(ald$N.{tid} - ald$y.{tid}.sum)")
  eval(parse(text = var_string))
}


#' Trial treatment effect with aggregate data
#' 
#' Calculate
#' \deqn{\log(\sum y_k (N_k - \sum y_k))}.
#' 
#' @param ald Aggregate-level data
#' @param tid Treatment label
#'
#' @return Value
#' @export
#'
trial_treatment_effect <- function(ald, tid) {
  var_string <- glue::glue("log(ald$y.{tid}.sum*(ald$N.{tid} - ald$y.{tid}.sum))")
  eval(parse(text = var_string))
}

  
#' Aggregate-level data mean and variance statistics
#'
#' @param ald Aggregate-level trial data
#' @param treatments Treatment labels; default `B`, `C`
#'
#' @return List of marginal treatment effect mean and variance
#' @seealso [marginal_treatment_effect()], [marginal_variance()]
#' @export
#'
ALD_stats <- function(ald, treatments = list("B", "C")) {
  list(mean = marginal_treatment_effect(ald, treatments),
       var = marginal_variance(ald, treatments))
}


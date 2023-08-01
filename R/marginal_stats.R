
#' marginal effect variance using the delta method
#'
#' \eqn{\frac{1}{n_C} + \frac{1}{n_{\bar{C}}} + \frac{1}{n_B} + \frac{1}{n_{\bar{B}}}}
#'
#' @param x 
#' @param trials 
#' @export
#' 
marginal_variance <- function(x, trials = list("B", "C")) {
  
  trial_vars <- purrr::map_dbl(trials, ~trial_variance(x, .x))
  sum(trial_vars)
}


#' B vs C marginal treatment effect from reported event counts
#'
#' \leqn{\log(n_B n_{\bar{C}})  -log(n_C n_{\bar{B}})}
#'
#' @param x 
#' @param trials
#' @export
#' 
marginal_treatment_effect <- function(x, trials = list("B", "C")) {
  trial_effect <- purrr::map_dbl(trials, ~trial_treatment_effect(x, .x))
  trial_effect[2] - trial_effect[1]
}


#' trial_variance
#'
#' @param x 
#' @param k 
#'
#' @return
#' @export
#'
trial_variance <- function(x, k) {
  var_string <- glue::glue("1/x$y.{k}.sum + 1/(x$N.{k} - x$y.{k}.sum)")
  eval(parse(text = var_string))
}


#' trial_treatment_effect
#'
#' @param x 
#' @param k 
#'
#' @return
#' @export
#'
trial_treatment_effect <- function(x, k) {
  var_string <- glue::glue("log(x$y.{k}.sum*(x$N.{k} - x$y.{k}.sum))")
  eval(parse(text = var_string))
}

  
#' ALD_stats
#'
#' @param data 
#' @param trials 
#'
#' @return
#' @export
#'
ALD_stats <- function(data, trials = list("B", "C")) {
  list(mean = marginal_treatment_effect(data, trials),
       var = marginal_variance(data, trials))
}


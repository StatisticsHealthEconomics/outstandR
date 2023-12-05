
#' marginal effect variance using the delta method
#'
#' \deqn{\frac{1}{n_C} + \frac{1}{n_{\bar{C}}} + \frac{1}{n_B} + \frac{1}{n_{\bar{B}}}}
#'
#' @param x x
#' @param trials Trial labels
#' @export
#' 
marginal_variance <- function(x, trials = list("B", "C")) {
  
  trial_vars <- purrr::map_dbl(trials, ~trial_variance(x, .x))
  sum(trial_vars)
}


#' B vs C marginal treatment effect from reported event counts
#'
#' \deqn{\log(n_B n_{\bar{C}})  -log(n_C n_{\bar{B}})}
#'
#' @param x x
#' @param trials Trial labels
#' @return Trial effect difference
#' @export
#' 
marginal_treatment_effect <- function(x, trials = list("B", "C")) {
  trial_effect <- purrr::map_dbl(trials, ~trial_treatment_effect(x, .x))
  trial_effect[2] - trial_effect[1]
}


#' Trial variance
#'
#' @param x x
#' @param k k
#'
#' @return Value
#' @export
#'
trial_variance <- function(x, k) {
  var_string <- glue::glue("1/x$y.{k}.sum + 1/(x$N.{k} - x$y.{k}.sum)")
  eval(parse(text = var_string))
}


#' Trial treatment effect
#'
#' @param x x
#' @param k k
#'
#' @return Value
#' @export
#'
trial_treatment_effect <- function(x, k) {
  var_string <- glue::glue("log(x$y.{k}.sum*(x$N.{k} - x$y.{k}.sum))")
  eval(parse(text = var_string))
}

  
#' Aggregate-level data statistics
#'
#' @param data Data
#' @param trials Trial labels; default `B`, `C`
#'
#' @return List of marginal treatment effect and variance
#' @export
#'
ALD_stats <- function(data, trials = list("B", "C")) {
  list(mean = marginal_treatment_effect(data, trials),
       var = marginal_variance(data, trials))
}


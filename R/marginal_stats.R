
#' marginal effect variance using the delta method
#'
#' \eqn{\frac{1}{n_C} + \frac{1}{n_{\bar{C}}} + \frac{1}{n_B} + \frac{1}{n_{\bar{B}}}}
#'
marginal_variance <- function(x)
  1/x$y.C.sum + 1/(x$N.C - x$y.C.sum) + 1/x$y.B.sum + 1/(x$N.B - x$y.B.sum)


#' B vs C marginal treatment effect from reported event counts
#'
#' \leqn{\log(n_B n_{\bar{C}})  -log(n_C n_{\bar{B}})}
#'
marginal_treatment_effect <- function(x)
  log(x$y.B.sum*(x$N.C - x$y.C.sum)) - log((x$y.C.sum*(x$N.B - x$y.B.sum)))


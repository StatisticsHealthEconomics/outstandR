
#' Aggregate-level data mean and variance statistics
#'
#' @param ald Aggregate-level trial data
#' @param treatments Treatment labels list; default `B`, `C`
#'
#' @return List of marginal treatment effect mean and variance
#' @seealso [marginal_treatment_effect()], [marginal_variance()]
#' @export
#'
ALD_stats <- function(strategy, ald, treatments = list("B", "C")) {
  list(mean = marginal_treatment_effect(ald, treatments, link = strategy$family$link),
       var = marginal_variance(ald, treatments, link = strategy$family$link))
}


#' Marginal effect variance using the delta method
#' 
#' Calculate
#' \deqn{\frac{1}{n_C} + \frac{1}{n_{\bar{C}}} + \frac{1}{n_B} + \frac{1}{n_{\bar{B}}}}.
#'
#' @param ald Aggregate-level data
#' @param treatments Treatment labels list; default _B_ vs _C_
#' @return Sum of variances
#' @export
#' 
marginal_variance <- function(ald, treatments = list("B", "C"), link) {
  trial_vars <- purrr::map_dbl(treatments, ~trial_variance(ald, .x, link))
  sum(trial_vars)
}


#' Marginal treatment effect from reported event counts
#' 
#' Calculate
#' \deqn{
#' \log\left( \frac{n_B/(N_B-n_B)}{n_C/(N_B-n_{B})} \right) = \log(n_B n_{\bar{C}}) - \log(n_C n_{\bar{B}})
#' }
#' where \eqn{\bar{C}} is the compliment of \eqn{C}
#' so e.g. \eqn{n_{\bar{C}} = N_C - n_c}.
#'
#' @param ald Aggregate-level data
#' @param treatments Treatment labels list; default _B_ vs _C_
#' @return Trial effect difference
#' @export
#' 
marginal_treatment_effect <- function(ald, treatments = list("B", "C"), link) {
  trial_effect <- purrr::map_dbl(treatments, ~trial_treatment_effect(ald, .x, link))
  trial_effect[2] - trial_effect[1]
}


#' Trial variance of the log-odds (logit) estimate with aggregate data
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
trial_variance <- function(ald, tid, link = "logit") {
  
  y <- ald[[paste0("y.", tid, ".sum")]]
  N <- ald[[paste0("N.", tid)]]
  
  link_transform_var(y, N, link)
}


#' Trial treatment effect with aggregate data
#' 
#' Calculate
#' \deqn{\log(\sum y_k (N_k - \sum y_k))}.
#' 
#' @param ald Aggregate-level data
#' @param tid Treatment label
#' @param link Link function; default "logit"
#'
#' @return Value
#' @export
#'
trial_treatment_effect <- function(ald, tid, link = "logit") {
  ##TODO: should this be instead i.e. log odds? it was * before
  # var_string <- glue::glue("log(ald$y.{tid}.sum / (ald$N.{tid} - ald$y.{tid}.sum))")
  
  # estimated probability
  p_hat <- ald[[paste0("y.", tid, ".sum")]] / ald[[paste0("N.", tid)]]
  
  link_transform(p_hat, link)
}


#' mean
#'
link_transform <- function(p, link) {
  if (link == "logit") {
    # log-OR
    return(qlogis(p))  # log(p / (1 - p))
  } else if (link == "log") {
    # log-Relative Risk (log-RR)
    return(log(p))
  }
}

#' variance
#'
link_transform_var <- function(y, N, link) {
  if (link == "logit") {
    # log-OR
    return(1/y + 1/(N - y))
  } else if (link == "log") {
    # log-RR
    return(1/y)
  }
}

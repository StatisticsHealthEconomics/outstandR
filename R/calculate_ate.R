#

#' Calculate average treatment effect
#'
#' @param ppv model prediction samples
#' @param family family object of the model
#'
#' @returns ATE
#' @export
#'
calculate_ate <- function(mean_A, mean_C, effect) {

  if (effect == "log_odds") {
    ate <- qlogis(mean_A) - qlogis(mean_C)
  } else if (effect == "risk_difference") {
    ate <- mean_A - mean_C
  } else if (effect == "delta_z") {
    ate <- qnorm(mean_A) - qnorm(mean_C)
  } else if (effect == "log_relative_risk_rare_events") {
    ate <- log(-log(1 - mean_A)) - log(-log(1 - mean_C))
  } else if (effect == "log_relative_risk") {  # Poisson log link
    ate <- log(mean_A) - log(mean_C)
  } else {
    stop("Unsupported link function.")
  }
  
  ate
}

#' @export
calculate_trial_variance <- function(ald, tid, effect, family) {

  if (family == "binomial") {
    return(
      calculate_trial_variance_binary(ald, tid, effect))
  } else if (family == "gaussian") {
    return(
      calculate_trial_variance_continuous(ald, tid, effect))
  } else {
    stop("family not recognised.")
  } 
}

#' @export
calculate_trial_variance_binary <- function(ald, tid, effect) {
  
  y <- ald[[paste0("y.", tid, ".sum")]]
  N <- ald[[paste0("N.", tid)]]
  
  if (effect == "log_odds") {
    ##TODO: double check these
    # res <- y/N + (N-y)/N
    res <- 1/y + 1/(N-y)
  } else if (effect == "log_relative_risk") {
    # using delta method
    res <- 1/y - 1/N
  } else if (effect == "risk_difference") {
    res <- y * (1 - y/N) / N
  } else if (effect == "delta_z") {
    res <- 1/y + 1/(N - y)
  } else if (effect == "log_relative_risk_rare_events") {
    res <- 1/y - 1/N
  } else {
    stop("Unsupported link function.")
  }
  
  res
}

#' @export
calculate_trial_variance_continuous <- function(ald, tid, effect) {
  
  ybar <- ald[[paste0("y.", tid, ".bar")]]
  ysd <- ald[[paste0("y.", tid, ".sd")]]
  N <- ald[[paste0("N.", tid)]]
  
  if (effect == "log_odds") {
    res <- pi^2/3 * (1/N)
  } else if (effect == "log_relative_risk") {
    message("log mean used\n")
    res <- log(ybar)
  } else if (effect == "risk_difference") {
    res <- (ysd^2)/N
  } else if (effect == "delta_z") {
    ##TODO:
    stop("Unsupported link function.")
  } else if (effect == "log_relative_risk_rare_events") {
    ##TODO:
    stop("Unsupported link function.")
  } else {
    stop("Unsupported link function.")
  }
  
  res
}

#' @export
calculate_trial_mean <- function(ald, tid, effect, family) {
  
  if (family == "binomial") {
    return(
      calculate_trial_mean_binary(ald, tid, effect))
  } else if (family == "gaussian") {
    return(
      calculate_trial_mean_continuous(ald, tid, effect))
  } else {
    stop("family not recognised.")
  } 
}

#' @export
calculate_trial_mean_binary <- function(ald, tid, effect) {
  
  y <- ald[[paste0("y.", tid, ".sum")]]
  N <- ald[[paste0("N.", tid)]]
  p <- y/N
  
  if (effect == "log_odds") {
    res <- qlogis(p)
  } else if (effect == "risk_difference") {
    res <- p
  } else if (effect == "delta_z") {
    res <- qnorm(p)
  } else if (effect == "log_relative_risk_rare_events") {
    res <- log(-log(1 - p))
  } else if (effect == "log_relative_risk") {
    res <- log(p)
  } else {
    stop("Unsupported link function.")
  }
  
  res
}

#' @export
calculate_trial_mean_continuous <- function(ald, tid, effect) {
  
  ybar <- ald[[paste0("y.", tid, ".bar")]]
  ysd <- ald[[paste0("y.", tid, ".sd")]]
  N <- ald[[paste0("N.", tid)]]
  
  if (effect == "log_odds") {
    message("log mean used\n")
    res <- log(ybar)
  } else if (effect == "risk_difference") {
    res <- ybar
  } else if (effect == "delta_z") {
    res <- ybar/ysd
  } else if (effect == "log_relative_risk_rare_events") {
    ##TODO:
    stop("Unsupported link function.")
  } else if (effect == "log_relative_risk") {
    message("log mean used\n")
    res <- log(ybar)
  } else {
    stop("Unsupported link function.")
  }
  
  res
}


#' Get treatment effect scale corresponding to a link function
#'
#' @export
get_treatment_effect <- function(link) {

  if (link == "logit") {
    rte <- "log_odds"
  } else if (link == "identity") {
    rte <- "risk_difference"
  } else if (link == "probit") {
    rte <- "delta_z"
  } else if (link == "cloglog") {  # binomial
    rte <- "log_relative_risk_rare_events"
  } else if (link == "log") {  # Poisson log link
    rte <- "log_relative_risk"
  } else {
    stop("Unsupported link function. Choose from
         'logit', 'identity', 'probit', 'cloglog', or 'log'.")
  }
  
  rte
}

calc_log_odds_ratio <- function(mean_A, mean_C) {
  qlogis(mean_A) - qlogis(mean_C)
}

calc_risk_difference <- function(mean_A, mean_C) {
  mean_A - mean_C
}

calc_delta_z <- function(mean_A, mean_C) {
  qnorm(mean_A) - qnorm(mean_C)
}

calc_log_relative_risk_rare_events <- function(mean_A, mean_C) {
  log(-log(1 - mean_A)) - log(-log(1 - mean_C))
}

calc_log_relative_risk <- function(mean_A, mean_C) {
  log(mean_A) - log(mean_C)
}

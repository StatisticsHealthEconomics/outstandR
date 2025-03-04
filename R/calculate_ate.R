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

#
calculate_trial_variance <- function(ald, tid, effect) {
  
  y <- ald[[paste0("y.", tid, ".sum")]]
  N <- ald[[paste0("N.", tid)]]
  
  if (effect == "log_odds") {
    res <- y/N + (N-y)/N
  } else if (effect == "log_relative_risk") {
    res <- 1/(N-y) - 1/N
  } else if (effect == "risk_difference") {
    ##TODO:
  } else if (effect == "delta_z") {
    ##TODO:
  } else if (effect == "log_relative_risk_rare_events") {
    ##TODO:
  } else {
    stop("Unsupported link function.")
  }
  
  res
}

#
calculate_trial_mean <- function(ald, tid, effect) {
  
  y <- ald[[paste0("y.", tid, ".sum")]]
  N <- ald[[paste0("N.", tid)]]
  p <- y/N
  
  if (effect == "log_odds") {
    # res <- log(p/(1-p)
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


#' Get treatment effect scale corresponding to a link function
#'
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

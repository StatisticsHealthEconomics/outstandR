#

#' Calculate Average Treatment Effect
#'
#' Computes the average treatment effect (ATE) based on the specified effect scale.
#'
#' @param mean_A,mean_C Mean of the outcome for the treatment and control
#' @param effect A character string specifying the effect scale. Options are:
#'   \describe{
#'     \item{"log_odds"}{Log-odds difference.}
#'     \item{"risk_difference"}{Risk difference.}
#'     \item{"delta_z"}{Probit scale difference (z-scores).}
#'     \item{"log_relative_risk_rare_events"}{Log relative risk for rare events.}
#'     \item{"log_relative_risk"}{Log relative risk.}
#'   }
#'
#' @return The computed average treatment effect on the specified scale.
#' @examples
#' \dontrun{
#' calculate_ate(mean_A = 0.7, mean_C = 0.5, effect = "log_odds")
#' calculate_ate(mean_A = 0.7, mean_C = 0.5, effect = "risk_difference")
#' }
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

#' Calculate trial variance
#'
#' Computes the variance of treatment effects for a trial based on the specified family distribution.
#'
#' @param ald Aggregate-level data.
#' @param tid Treatment identifier used to extract relevant columns from `ald`.
#' @param effect A character string specifying the effect scale (e.g., "log_odds", "risk_difference").
#' @param family A character string specifying the model family (e.g., "binomial", "gaussian").
#'
#' @return The computed variance of treatment effects.
#' @examples
#' \dontrun{
#' ald <- data.frame(y.B.sum = c(10), N.B = c(100))
#' calculate_trial_variance(ald, tid = "B", effect = "log_odds", family = "binomial")
#' }
#' @export
calculate_trial_variance <- function(ald, tid, effect, family) {
  
  if (family == "binomial") {
    return(
      calculate_trial_variance_binary(ald, tid, effect))
  } else if (family == "gaussian") {
    return(
      calculate_trial_variance_continuous(ald, tid, effect))
  }
  
  stop("family not recognised.")
}

#' @export
calculate_trial_variance_binary <- function(ald, tid, effect) {
  
  y <- ald[[paste0("y.", tid, ".sum")]]
  N <- ald[[paste0("N.", tid)]]
  
  effect_functions <- list(
    "log_odds" = function() 1/y + 1/(N-y),
    "log_relative_risk" = function() 1/y - 1/N,
    "risk_difference" = function() y * (1 - y/N) / N,
    "delta_z" = function() 1/y + 1/(N - y),
    "log_relative_risk_rare_events" = function() 1/y - 1/N
  )
  
  if (!effect %in% names(effect_functions)) {
    stop(paste0("Unsupported effect function. Choose from ",
                names(effect_functions)))
  }
  
  effect_functions[[effect]]()
}

#' @export
calculate_trial_variance_continuous <- function(ald, tid, effect) {
  
  ybar <- ald[[paste0("y.", tid, ".bar")]]
  ysd <- ald[[paste0("y.", tid, ".sd")]]
  N <- ald[[paste0("N.", tid)]]
  
  effect_functions <- list(
    "log_odds" = function() pi^2/3 * (1/N),
    "log_relative_risk" = function() {
      message("log mean used\n")
      log(ybar)
    },
    "risk_difference" = function() (ysd^2)/N
  )
  
  if (!effect %in% names(effect_functions)) {
    stop(paste0("Unsupported effect function. Choose from ",
                names(effect_functions)))
  }
  
  effect_functions[[effect]]()
}

#' @export
calculate_trial_mean <- function(ald, tid, effect, family) {
  
  if (family == "binomial") {
    return(
      calculate_trial_mean_binary(ald, tid, effect))
  } else if (family == "gaussian") {
    return(
      calculate_trial_mean_continuous(ald, tid, effect))
  }
  
  stop("family not recognised.")
}

#' @export
calculate_trial_mean_binary <- function(ald, tid, effect) {
  
  y <- ald[[paste0("y.", tid, ".sum")]]
  N <- ald[[paste0("N.", tid)]]
  p <- y/N
  
  effect_fns <- list(
    log_odds = function() qlogis(p),
    risk_difference = function() p,
    delta_z = function() qnorm(p),
    log_relative_risk_rare_events = function() log(-log(1 - p)),
    log_relative_risk = function() log(p)
  )
  
  if (!effect %in% names(effect_fns)) {
    stop(paste0("Unsupported link function. Choose from ",
                names(effect_fns)))
  }
  
  effect_fns[[effect]]()
}

#' @export
calculate_trial_mean_continuous <- function(ald, tid, effect) {
  
  ybar <- ald[[paste0("y.", tid, ".bar")]]
  ysd <- ald[[paste0("y.", tid, ".sd")]]
  N <- ald[[paste0("N.", tid)]]
  
  effect_fns <- list(
    log_odds = function() {
      message("log mean used\n")
      log(ybar)
    },
    risk_difference = function() ybar,
    delta_z = function() ybar / ysd,
    log_relative_risk = function() {
      message("log mean used\n")
      log(ybar)
    }
    # log_relative_risk_rare_events intentionally unsupported
  )
  
  if (!effect %in% names(effect_fns)) {
    stop(paste0("Unsupported link function. Choose from ",
         names(effect_fns)))
  }
  
  effect_fns[[effect]]()
}


#' Get treatment effect scale corresponding to a link function
#'
#' Maps a given link function to its corresponding treatment effect scale.
#'
#' @param link A character string specifying the link function. Options are:
#'   \describe{
#'     \item{"logit"}{Log-odds scale.}
#'     \item{"identity"}{Risk difference.}
#'     \item{"probit"}{Probit scale.}
#'     \item{"cloglog"}{Log relative risk for rare events.}
#'     \item{"log"}{Log relative risk.}
#'   }
#'
#' @return A character string representing the treatment effect scale.
#' @examples
#' \dontrun{
#' get_treatment_effect(link = "logit")
#' get_treatment_effect(link = "identity")
#' }
#' @export
get_treatment_effect <- function(link) {
  
  link_map <- list(
    logit = "log_odds",
    identity = "risk_difference",
    probit = "delta_z",
    cloglog = "log_relative_risk_rare_events",
    log = "log_relative_risk"
  )
  
  if (!link %in% names(link_map)) {
    stop(paste0("Unsupported link function. Choose from ",
         names(link_map)))
  }
  
  link_map[[link]]
}

# individual effects

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

#' @keywords internal
continuity_correction <- function(ald,
                                  treatments = list("B", "C"),
                                  correction = 0.5) {
  for (t in treatments) {
    y_name <- paste0("y.", t, ".sum")
    N_name <- paste0("N.", t)
    
    y <- ald[[y_name]]
    N <- ald[[N_name]]
    failures <- N - y
    
    # apply correction to both successes and failures
    ald[[y_name]] <- y + correction
    ald[[N_name]] <- N + 2 * correction  # since both y and failures are adjusted
  }
  
  ald
}

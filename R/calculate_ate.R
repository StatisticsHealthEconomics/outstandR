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
#' calculate_ate(mean_A = 0.7, mean_ref = 0.5, effect = "log_odds")
#' calculate_ate(mean_A = 0.7, mean_ref = 0.5, effect = "risk_difference")
#' }
#' @export
#'
calculate_ate <- function(mean_A, mean_ref, effect) {
  
  if (effect == "log_odds") {
    ate <- qlogis(mean_A) - qlogis(mean_ref)
  } else if (effect == "risk_difference") {
    ate <- mean_A - mean_ref
  } else if (effect == "delta_z") {
    ate <- qnorm(mean_A) - qnorm(mean_ref)
  } else if (effect == "log_relative_risk_rare_events") {
    ate <- log(-log(1 - mean_A)) - log(-log(1 - mean_ref))
  } else if (effect == "log_relative_risk") {  # Poisson log link
    ate <- log(mean_A) - log(mean_ref)
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
#' ald <- data.frame(trt = c("B","C","B","C"),
#'                   variable = c(NA, NA, "y", "y"),
#'                   statistic = c("N", "N", "sum", "sum"),
#'                   value = c(100, 100, 50, 60)
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
  
  y <- dplyr::filter(
    ald,
    variable == "y",
    trt == tid,
    statistic == "sum")$value
  
  N <- dplyr::filter(
    ald,
    trt == tid,
    statistic == "N")$value
  
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
  
  ybar <- dplyr::filter(
    ald,
    variable == "y",
    trt == tid,
    statistic == "mean")$value
  
  ysd <- dplyr::filter(
    ald,
    variable == "y",
    trt == tid,
    statistic == "sd")$value
  
  N <- dplyr::filter(
    ald,
    trt == tid,
    statistic == "N")$value
  
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
  
  y <- dplyr::filter(
    ald,
    variable == "y",
    trt == tid,
    statistic == "sum")$value
  
  N <- dplyr::filter(
    ald,
    trt == tid,
    statistic == "N")$value
  
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
  
  ybar <- dplyr::filter(
    ald,
    variable == "y",
    trt == tid,
    statistic == "mean")$value
  
  ysd <- dplyr::filter(
    ald,
    variable == "y",
    trt == tid,
    statistic == "sd")$value
  
  N <- dplyr::filter(
    ald,
    trt == tid,
    statistic == "N")$value
  
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

calc_log_odds_ratio <- function(mean_A, mean_ref) {
  qlogis(mean_A) - qlogis(mean_ref)
}

calc_risk_difference <- function(mean_A, mean_ref) {
  mean_A - mean_ref
}

calc_delta_z <- function(mean_A, mean_ref) {
  qnorm(mean_A) - qnorm(mean_ref)
}

calc_log_relative_risk_rare_events <- function(mean_A, mean_ref) {
  log(-log(1 - mean_A)) - log(-log(1 - mean_ref))
}

calc_log_relative_risk <- function(mean_A, mean_ref) {
  log(mean_A) - log(mean_ref)
}

#' @keywords internal
continuity_correction <- function(ald,
                                  treatments = list("B", "C"),
                                  correction = 0.5) {
  # check if correction is needed in any group
  needs_correction <- 
    ald |> 
    dplyr::filter((variable == "y" & statistic == "sum") | statistic == "N") |>
    dplyr::group_by(trt, variable) |>
    spread(statistic, value) |>  # Spread sd and N into separate columns
    dplyr::mutate(need_contcorr = case_when(
      sum == 0 ~ TRUE,
      sum == N ~ TRUE,
      TRUE ~ FALSE
    )) |> pull() |> any()
  
  if (!needs_correction) {
    return(ald)
  }
  
  message(sprintf(
    "Applying continuity correction: %d", correction
  ))
  
  ald_corrected <- ald %>%
    dplyr::mutate(
      value = case_when(
        statistic == "sum" & variable == "y" ~ value + correction,
        statistic == "N" ~ value + 2 * correction,
        TRUE ~ value
      )
    )
  
  ald_corrected
}

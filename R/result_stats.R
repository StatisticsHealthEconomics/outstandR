
#' Calculate and arrange result statistics
#'
#' Combining output from aggregate level data studies BC and 
#' adjusted individual level data studies AC into a single object.
#'
#' @param ipd_stats,ald_stats Trial data for IPD and ALD
#' @param CI Confidence interval level, i.e. 1-alpha; default 0.95
#' @param scale Relative treatment effect scale; default "mean_difference"
#'
#' @returns List of ITC output:
#' * `contrasts`: A list for relative effects containing:
#'     * `means`
#'     * `variances`
#'     * `CI`
#' * `absolute`: A list for absolute effects containing:
#'     * `means`
#'     * `variances`
#'     * `CI`
#'     
#' @keywords internal
#'
result_stats <- function(ipd_stats,
                         ald_stats,
                         CI = 0.95,
                         scale = "mean_difference") {

  AC_contrasts <- ipd_stats$contrasts
  AC_absolute <- ipd_stats$absolute
  
  # contrasts ---
  
  contrasts <- list(
    AB = AC_contrasts$mean - ald_stats$mean,
    AC = AC_contrasts$mean,
    BC = ald_stats$mean)
  
  contrast_variances <- list(
    AB = AC_contrasts$var + ald_stats$var,
    AC = AC_contrasts$var,
    BC = ald_stats$var)
  
  contrast_ci <- list(
    AB = calc_ci(mean_val = contrasts$AB, sd_val = sqrt(contrast_variances$AB), level = CI),
    AC = calc_ci(mean_val = contrasts$AC, sd_val = sqrt(contrast_variances$AC), level = CI),
    BC = calc_ci(mean_val = contrasts$BC, sd_val = sqrt(contrast_variances$BC), level = CI)
  )
  
  ##TODO: MIM CI
  # lci.mim <- coef_est + qt(0.025, df = model$nu) * sqrt(var_est)
  # uci.mim <- coef_est + qt(0.975, df = model$nu) * sqrt(var_est)
  
  # absolute values ---
  
  abs_B <- calc_absolute_B(
    mean_C   = unname(AC_absolute$mean[2]),
    var_C    = unname(AC_absolute$var[2]),
    delta_BC = ald_stats$mean,
    var_BC   = ald_stats$var,
    scale    = scale
  )
  
  absolute <- list(
    A = unname(AC_absolute$mean[1]),
    B = abs_B$mean,
    C = unname(AC_absolute$mean[2])
  )
  
  absolute_var <- list(
    A = unname(AC_absolute$var[1]),
    B = abs_B$var,
    C = unname(AC_absolute$var[2])
  )
  
  absolute_ci <- list(
    A = calc_ci(mean_val = absolute$A, sd_val = sqrt(absolute_var$A), level = CI),
    B = calc_ci(mean_val = absolute$B, sd_val = sqrt(absolute_var$B), level = CI),
    C = calc_ci(mean_val = absolute$C, sd_val = sqrt(absolute_var$C), level = CI)
  )
  
  list(
    contrasts = list(
      means = contrasts,
      variances = contrast_variances,
      CI = contrast_ci),
    absolute = list(
      means = absolute,
      variances = absolute_var,
      CI = absolute_ci
    ))
}

#' Calculate absolute mean and variance for ALD comparator arm B on the natural outcome scale
#' @keywords internal
calc_absolute_B <- function(mean_C, var_C, delta_BC, var_BC, scale = NULL) {
  if (is.null(scale) || scale %in% c("mean_difference", "risk_difference", "rate_difference")) {
    return(list(
      mean = unname(mean_C + delta_BC),
      var  = unname(var_C + var_BC)
    ))
  }
  
  if (scale == "log_odds") {
    eta_C <- qlogis(mean_C)
    var_eta_C <- var_C / ((mean_C * (1 - mean_C))^2)
    eta_B <- eta_C + delta_BC
    var_eta_B <- var_eta_C + var_BC
    mean_B <- plogis(eta_B)
    var_B <- ((mean_B * (1 - mean_B))^2) * var_eta_B
    return(list(mean = unname(mean_B), var = unname(var_B)))
  }
  
  if (scale == "log_relative_risk") {
    eta_C <- log(mean_C)
    var_eta_C <- var_C / (mean_C^2)
    eta_B <- eta_C + delta_BC
    var_eta_B <- var_eta_C + var_BC
    mean_B <- exp(eta_B)
    var_B <- (mean_B^2) * var_eta_B
    return(list(mean = unname(mean_B), var = unname(var_B)))
  }
  
  if (scale == "delta_z") {
    eta_C <- qnorm(mean_C)
    var_eta_C <- var_C / (dnorm(eta_C)^2)
    eta_B <- eta_C + delta_BC
    var_eta_B <- var_eta_C + var_BC
    mean_B <- pnorm(eta_B)
    var_B <- (dnorm(eta_B)^2) * var_eta_B
    return(list(mean = unname(mean_B), var = unname(var_B)))
  }
  
  if (scale == "log_relative_risk_rare_events") {
    eta_C <- log(-log(1 - mean_C))
    deriv_C <- 1 / ((1 - mean_C) * (-log(1 - mean_C)))
    var_eta_C <- (deriv_C^2) * var_C
    eta_B <- eta_C + delta_BC
    var_eta_B <- var_eta_C + var_BC
    mean_B <- 1 - exp(-exp(eta_B))
    deriv_B <- 1 / ((1 - mean_B) * (-log(1 - mean_B)))
    var_B <- var_eta_B / (deriv_B^2)
    return(list(mean = unname(mean_B), var = unname(var_B)))
  }
  
  list(
    mean = unname(mean_C + delta_BC),
    var  = unname(var_C + var_BC)
  )
}

#' @keywords internal
calc_ci <- function(mean_val, sd_val, level = 0.95) {
  upper <- 0.5 + level/2
  ci_range <- c(1 - upper, upper)
  z_vals <- qnorm(ci_range)
  
  mean_val + z_vals*as.vector(sd_val)
}

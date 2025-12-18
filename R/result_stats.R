
#' Calculate and arrange result statistics
#'
#' Combining output from aggregate level data studies BC and 
#' adjusted individual level data studies AC into a single object.
#'
#' @param ipd_stats,ald_stats 
#' @param CI Confidence interval level, i.e. 1-alpha; default 0.95
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
                         CI = 0.95) {

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
  
  absolute <- list(
    A = AC_absolute$mean["mean_A"],
    # B = AB_absolute$mean["mean_B"],
    C = AC_absolute$mean["mean_C"]
  )
  
  absolute_var <- list(
    A = AC_absolute$var["mean_A"],
    # B = AB_absolute$var["mean_B"],
    C = AC_absolute$var["mean_C"]
  )
  
  ##TODO:
  absolute_ci <- list(
    AB = NA,
    AC = NA,
    BC = NA
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

#' @keywords internal
calc_ci <- function(mean_val, sd_val, level = 0.95) {
  upper <- 0.5 + level/2
  ci_range <- c(1 - upper, upper)
  z_vals <- qnorm(ci_range)
  
  mean_val + z_vals*as.vector(sd_val)
}

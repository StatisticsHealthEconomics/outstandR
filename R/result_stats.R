
#' Result Statistics
#'
#' Combining output from aggregate level data studies BC and 
#' adjusted individual level data studies AC into a single object.
#'
#' @param AC_stats,BC_stats 
#' @param CI Confidence interval 1-alpha
#'
#' @returns List
#' @keywords internal
#'
result_stats <- function(AC_stats,
                         BC_stats,
                         CI = 0.95) {
  upper <- 0.5 + CI/2
  ci_range <- c(1-upper, upper)
  z_vals <- qnorm(ci_range)
  
  AC_contrasts <- AC_stats$contrasts
  AC_absolute <- AC_stats$absolute
  
  # contrasts
  
  contrasts <- list(
    AB = AC_contrasts$mean - BC_stats$mean,
    AC = AC_contrasts$mean,
    BC = BC_stats$mean)
  
  contrast_variances <- list(
    AB = AC_contrasts$var + BC_stats$var,
    AC = AC_contrasts$var,
    BC = BC_stats$var)
  
  contrast_ci <- list(
    AB = contrasts$AB + z_vals*as.vector(sqrt(contrast_variances$AB)),
    AC = contrasts$AC + z_vals*as.vector(sqrt(contrast_variances$AC)),
    BC = contrasts$BC + z_vals*as.vector(sqrt(contrast_variances$BC)))
  
  # absolute values
  
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
  
  list(
    contrasts = list(
      means = contrasts,
      variances = contrast_variances,
      CI = contrast_ci),
    absolute = list(   ##TODO:
      means = absolute,
      variances = absolute_var
      # CI = contrast_ci
    ))
}

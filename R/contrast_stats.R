
#
contrast_stats <- function(AC_stats,
                           BC_stats,
                           CI = 0.95) {
  upper <- 0.5 + CI/2
  ci_range <- c(1-upper, upper)
  
  contrasts <- list(
    AB = AC_stats$mean - BC_stats$mean,
    AC = AC_stats$mean,
    BC = BC_stats$mean)
  
  contrast_variances <- list(
    AB = AC_stats$var + BC_stats$var,
    AC = AC_stats$var,
    BC = BC_stats$var)
  
  contrast_ci <- list(
    AB = contrasts$AB + qnorm(ci_range)*as.vector(sqrt(contrast_variances$AB)),
    AC = contrasts$AC + qnorm(ci_range)*as.vector(sqrt(contrast_variances$AC)),
    BC = contrasts$BC + qnorm(ci_range)*as.vector(sqrt(contrast_variances$BC)))
  
  list(contrasts = contrasts,
       variances = contrast_variances,
       CI = contrast_ci)
}

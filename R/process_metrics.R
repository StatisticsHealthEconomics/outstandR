
#' Compute performance metrics for a given method
#' 
#' @param means Mean
#' @param variances Variances
#' @param truth Truth
#' 
#' @keywords internal
#' 
process_metrics <- function(means, variances, truth) {
  # remove NAs (an issue for MAIC in Scenario 7, separation issues)
  NAs <- is.na(means)
  means <- means[!NAs]  
  variances <- variances[!NAs]
  number.NAs <- sum(NAs) # number that do not converge
  replicates <- length(means)
  bias.metric <- bias(means, truth)
  bias.metric.mcse <- bias.mcse(means)
  # # bias is equal to estimated marginal treatment effect (truth is zero)
  # ate <- mean(means)
  # ate.mcse <- mcse.estimate(means)
  mae.metric <- mae(means, truth)
  abs.err <- abs(means - truth)
  mae.mcse <- mcse.estimate(abs.err)
  mse.metric <- mse(means, truth) 
  mse.metric.mcse <- mse.mcse(means, truth) 
  # construct Wald-type interval estimates using normal distribution
  lci <- means + qnorm(0.025)*sqrt(variances)
  uci <- means + qnorm(0.975)*sqrt(variances)
  lci.mean <- mean(lci)
  lci.mcse <- mcse.estimate(lci)
  uci.mean <- mean(uci)
  uci.mcse <- mcse.estimate(uci)
  cov <- coverage(lci, uci, truth)
  cov.mcse <- coverage.mcse(cov, replicates)
  empse.metric <- empse(means)
  empse.metric.mcse <- empse.mcse(empse.metric, replicates)
  vr <- var.ratio(means, sqrt(variances))
  vr.mcse <- var.ratio.mcse(avg.se=mean(sqrt(variances)), 
                            emp.se=empse.metric,
                            var.avg.se=mcse.estimate(sqrt(variances))^2,
                            var.emp.se=empse.metric.mcse^2)
  tibble::lst(
    bias.metric,
    bias.metric.mcse,
    lci.mean,
    lci.mcse,
    uci.mean,
    uci.mcse,
    vr,
    vr.mcse,
    cov,
    cov.mcse,
    empse.metric,
    empse.metric.mcse,
    mse.metric,
    mse.metric.mcse,
    mae.metric,
    mae.mcse,
    number.NAs)
} 


#' @name IPD_stats
#' @title Calculate individual-level patient data statistics
#' 
#' @description
#' Separate methods for each approach
#' MAIC, STC, G-computation via MLE or Bayesian inference.
#' 
#' @param strategy A list corresponding to different approaches
#' @template args-ipd
#' @template args-ald
#' @param ... Additional arguments
#' 
#' @return Mean and variance values
#' @export
#' 
IPD_stats <- function(strategy, ipd, ald, scale, ...)
  UseMethod("IPD_stats", strategy)


#' @rdname IPD_stats
#' @importFrom utils methods
#' @export
IPD_stats.default <- function(...) {
  strategy_classes <- sub("IPD_stats\\.", "", methods(IPD_stats)[-1])
  avail_strategies <- paste0("strategy_", strategy_classes, "()", collapse = ", ")
  stop(paste0("strategy not available. Select from ", avail_strategies))
}


#' @rdname IPD_stats
#' @section Matching-adjusted indirect comparison statistics:
#' Marginal _A_ vs _C_ treatment effect estimates
#' using bootstrapping sampling.
#' @importFrom boot boot
#' @export
#' 
IPD_stats.maic <- function(strategy,
                           ipd, ald,
                           scale) {
  
  maic_boot <- boot::boot(data = ipd,
                          statistic = maic.boot,
                          R = strategy$R,
                          formula = strategy$formula,
                          family = strategy$family,
                          ald = ald)

  hat.delta.AC <- calculate_ate(maic_boot$mean_A, maic_boot$mean_C,
                                effect = scale)
  
  treat_nm <- get_treatment_name(strategy$formula)
  coef_est <- mean(hat.delta.AC)
  var_est <- var(hat.delta.AC)

  list(mean = coef_est,
       var = var_est)
}


#' @rdname IPD_stats
#' @section Simulated treatment comparison statistics: 
#' IPD from the _AC_ trial are used to fit a regression model describing the
#' observed outcomes \eqn{y} in terms of the relevant baseline characteristics \eqn{x} and
#' the treatment variable \eqn{z}.
#' @importFrom stats glm
#' @export
#' 
IPD_stats.stc <- function(strategy,
                          ipd, ald,
                          scale) {

  # centre covariates
  term.labels <- attr(terms(strategy$formula), "term.labels")
  centre_vars <- gsub("trt:", "", term.labels[grepl(":", term.labels)])

  ipd[, centre_vars] <- scale(ipd[, centre_vars], scale = FALSE)
  
  fit <- glm(formula = strategy$formula,
             family = strategy$family,
             data = ipd)
  
  # extract model coefficients
  coef_fit <- coef(fit)
  
  mean_A <- family$linkinv(coef_fit[1])                # probability for control group
  mean_C <- family$linkinv(coef_fit[1] + coef_fit[2])  # probability for treatment group
  
  hat.delta.AC <- calculate_ate(mean_A, mean_C,
                                effect = scale)
  
  treat_nm <- get_treatment_name(strategy$formula)
  coef_est <- mean(hat.delta.AC)
  var_est <- var(hat.delta.AC)
  
  # coef_est <- coef(fit)[treat_nm]
  # var_est <- vcov(fit)[treat_nm, treat_nm]
  
  list(mean = coef_est,
       var = var_est)
}


#' @rdname IPD_stats
#' @section G-computation maximum likelihood statistics:
#' Compute a non-parametric bootstrap with \eqn{R=1000} resamples.
#' @importFrom boot boot
#' @export
#'
IPD_stats.gcomp_ml <- function(strategy,
                               ipd, ald,
                               scale) {

  gcomp_boot <- boot::boot(data = ipd,
                           statistic = gcomp_ml.boot,
                           R = strategy$R,
                           formula = strategy$formula,
                           family = strategy$family,
                           ald = ald)
  
  hat.delta.AC <- calculate_ate(gcomp_boot$A_mean, gcomp_boot$C_mean,
                                effect = scale)
  
  treat_nm <- get_treatment_name(strategy$formula)
  coef_est <- mean(hat.delta.AC)
  var_est <- var(hat.delta.AC)
  
  list(mean = coef_est,
       var = var_est)
}


#' @rdname IPD_stats
#' @section G-computation Bayesian statistics:
#' Using Stan, compute marginal log-odds ratio for _A_ vs _C_ for each MCMC sample
#' by transforming from probability to linear predictor scale.
#' @importFrom stats qlogis
#' @export
#'
IPD_stats.gcomp_stan <- function(strategy,
                                 ipd, ald,
                                 scale) {
  
  ppv <- gcomp_stan(formula = strategy$formula,
                    family = strategy$family,
                    ipd = ipd, ald = ald)

  # posterior means for each treatment group
  mean_A <- rowMeans(ppv$y.star.A)
  mean_C <- rowMeans(ppv$y.star.C)
  
  hat.delta.AC <- calculate_ate(mean_A, mean_C, effect = scale)
  
  treat_nm <- get_treatment_name(strategy$formula)
  coef_est <- mean(hat.delta.AC)
  var_est <- var(hat.delta.AC)
  
  list(mean = coef_est,
       var = var_est)
} 


#' @rdname IPD_stats
#' @section Multiple imputation marginalisation:
#' Using Stan, compute marginal log-odds ratio for _A_ vs _C_ for each MCMC sample
#' by transforming from probability to linear predictor scale. Approximate by 
#' using imputation and combining estimates using Rubin's rules, in contrast to [IPD_stats.gcomp_stan()].
#' @import stats
#' @export
#'
IPD_stats.mim <- function(strategy,
                          ipd, ald,
                          scale) {
  
  mis_res <- mim(formula = strategy$formula,
                 family = strategy$family,
                 ipd, ald)
  
  M <- mis_res$M
  
  # quantities originally defined by Rubin (1987) for multiple imputation
  bar.delta <- mean(mis_res$hats.delta)  # average of treatment effect point estimates
  bar.v <- mean(mis_res$hats.v)          # "within" variance (average of variance point estimates)
  b <- var(mis_res$hats.delta)           # "between" variance (sample variance of point estimates)
  
  # pooling: average of point estimates is marginal log odds ratio
  coef_est <- bar.delta
  
  var_est <- var_by_pooling(M, bar.v, b)
  nu <- wald_type_interval(M, bar.v, b)
  
  ##TODO: how are these used?
  lci.Delta <- coef_est + qt(0.025, df = nu) * sqrt(var_est)
  uci.Delta <- coef_est + qt(0.975, df = nu) * sqrt(var_est)
  
  treat_nm <- get_treatment_name(strategy$formula)
  
  list(mean = coef_est,
       var = var_est)
} 


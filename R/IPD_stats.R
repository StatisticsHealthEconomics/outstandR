
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
IPD_stats <- function(strategy, ipd, ald, ...)
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
                           ipd, ald) {
  # maic.boot(data = data,
  #           indices = 1:nrow(data),
  #           formula = strategy$formula,
  #           dat_ALD = strategy$dat_ALD)
  
  maic_boot <- boot::boot(data = ipd,
                          statistic = maic.boot,
                          R = strategy$R,
                          formula = strategy$formula,
                          ald = ald)

  treat_nm <- get_treatment_name(strategy$formula)
  coef_est <- mean(maic_boot$t)
  var_est <- var(maic_boot$t)
  
  # compute baseline probability in control group (P0)
  newdat <- ipd[ipd[[treat_nm]] == 0, ]
  P0 <- mean(predict(fit, newdata = newdat, type = "response"))
  
  converted_effect <- convert_effect(coef_est, from_scale, to_scale, P0)
  
  list(mean = converted_effect,
       var = var_est)  ##TODO: variance conversion
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
                          ipd, ald) {

  # centre covariates
  term.labels <- attr(terms(strategy$formula), "term.labels")
  centre_vars <- gsub("trt:", "", term.labels[grepl(":", term.labels)])

  ipd[, centre_vars] <- scale(ipd[, centre_vars], scale = FALSE)
  
  fit <- glm(strategy$formula,
             data = ipd,
             family = strategy$family)
  
  treat_nm <- get_treatment_name(strategy$formula)
  coef_est <- coef(fit)[treat_nm]
  var_est <- vcov(fit)[treat_nm, treat_nm]
  
  # compute baseline probability in control group (P0)
  newdat <- ipd[ipd[[treat_nm]] == 0, ]
  P0 <- mean(predict(fit, newdata = newdat, type = "response"))
  
  converted_effect <- convert_effect(coef_est, from_scale, to_scale, P0)
  
  # fitted treatment coefficient is relative A vs C conditional effect
  list(mean = converted_effect,
       var = var_est)  ##TODO: variance conversion
}


#' @rdname IPD_stats
#' @section G-computation maximum likelihood statistics:
#' Compute a non-parametric bootstrap with \eqn{R=1000} resamples.
#' @importFrom boot boot
#' @export
#'
IPD_stats.gcomp_ml <- function(strategy,
                               ipd, ald) {

  AC_maic_boot <- boot::boot(data = ipd,
                             statistic = gcomp_ml.boot,
                             R = strategy$R,
                             formula = strategy$formula,
                             family = strategy$family,
                             ald = ald)
  
  treat_nm <- get_treatment_name(strategy$formula)
  coef_est <- mean(AC_maic_boot$t)
  var_est <- var(AC_maic_boot$t)
  
  # compute baseline probability in control group (P0)
  newdat <- ipd[ipd[[treat_nm]] == 0, ]
  P0 <- mean(predict(fit, newdata = newdat, type = "response"))
  
  converted_effect <- convert_effect(coef_est, from_scale, to_scale, P0)
  
  list(mean = converted_effect,
       var = var_est)  ##TODO: variance conversion
}


#' @rdname IPD_stats
#' @section G-computation Bayesian statistics:
#' Using Stan, compute marginal log-odds ratio for _A_ vs _C_ for each MCMC sample
#' by transforming from probability to linear predictor scale.
#' @importFrom stats qlogis
#' @export
#'
IPD_stats.gcomp_stan <- function(strategy,
                                 ipd, ald) {
  
  ppv <- gcomp_stan(formula = strategy$formula,
                    family = strategy$family,
                    ipd = ipd, ald = ald)

  # posterior means for each treatment group
  mean_A <- rowMeans(ppv$y.star.A)
  mean_C <- rowMeans(ppv$y.star.C)
  
  hat.delta.AC <- calculate_ate(mean_A, mean_C, family = strategy$family)
  
  treat_nm <- get_treatment_name(strategy$formula)
  coef_est <- mean(hat.delta.AC)
  var_est <- var(hat.delta.AC)
  
  # compute baseline probability in control group (P0)
  newdat <- ipd[ipd[[treat_nm]] == 0, ]
  P0 <- mean(predict(fit, newdata = newdat, type = "response"))
  
  converted_effect <- convert_effect(coef_est, from_scale, to_scale, P0)
  
  list(mean = converted_effect,
       var = var_est)  ##TODO: variance conversion
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
                          ipd, ald) {
  
  mis_res <- mim(formula = strategy$formula,
                 family = strategy$family,
                 ipd, ald)
  
  M <- mis_res$M
  
  # quantities originally defined by Rubin (1987) for multiple imputation
  bar.delta <- mean(mis_res$hats.delta)  # average of treatment effect point estimates
  bar.v <- mean(mis_res$hats.v)          # "within" variance (average of variance point estimates)
  b <- var(mis_res$hats.delta)           # "between" variance (sample variance of point estimates)
  
  # pooling: average of point estimates is marginal log odds ratio
  hat.Delta <- bar.delta
  
  hat.var.Delta <- var_by_pooling(M, bar.v, b)
  nu <- wald_type_interval(M, bar.v, b)
  
  lci.Delta <- hat.Delta + qt(0.025, df = nu) * sqrt(hat.var.Delta)
  uci.Delta <- hat.Delta + qt(0.975, df = nu) * sqrt(hat.var.Delta)
  
  treat_nm <- get_treatment_name(strategy$formula)
  coef_est <- hat.Delta
  var_est <- hat.var.Delta
  
  # compute baseline probability in control group (P0)
  newdat <- ipd[ipd[[treat_nm]] == 0, ]
  P0 <- mean(predict(fit, newdata = newdat, type = "response"))
  
  converted_effect <- convert_effect(coef_est, from_scale, to_scale, P0)
  
  list(mean = converted_effect,
       var = var_est)  ##TODO: variance conversion
} 


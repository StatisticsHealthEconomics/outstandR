
# create class for each approach

strategy_maic <- function() {
  new_strategy("maic")
}

strategy_stcc <- function() {
  new_strategy("stc")
}

strategy_gcomp_ml <- function() {
  new_strategy("gcomp_ml")
}

strategy_gcomp_stan <- function() {
  new_strategy("gcomp_stan")
}

strategy_gcomp_stan <- function(strategy) {
  structure(strategy, class = strategy)
}


#' main wrapper
#' 
hat_Delta_stats <- function(AC.IPD, BC.ALD, strategy, ...) {
  
  AC_hat_Delta_stats <- IPD_stats(strategy, data = AC.IPD, ...) 
  
  hat.mean.Delta.BC <- marginal_treatment_effect(BC.ALD)
  hat.var.Delta.BC <- marginal_variance(BC.ALD)
  
  stats <- c(hat.mean.Delta.AB = AC_hat_Delta_stats$mean - hat.mean.Delta.BC,
             hat.var.Delta.AB = AC_hat_Delta_stats$var + hat.var.Delta.BC)
  
  structure(stats, class = c("mimR", class(stats)))
}


#' individual level data statistics
#'
IPD_stats <- function(strategy, data, ...)
  UseMethod("IPD_stats", strategy)

#
IPD_stats.default <- function() {
  stop("strategy not available.")
}


#' marginal A vs C treatment effect estimates
#' using bootstrapping
#'
IPD_stats.maic <- function(strategy,
                           formula = as.formula("y ~ trt"),
                           data = AC.IPD,
                           R = 1000) {
  maic_boot <- boot::boot(data = data,
                          statistic = maic.boot,
                          R = R,
                          formula = formula)
  
  list(mean =  mean(maic_boot$t),
       var = var(maic_boot$t))
}


#'
IPD_stats.stc <- function(strategy,
                          formula =
                            as.formula("y ~ X3 + X4 +
                                   trt*I(X1-BC.ALD$mean.X1) +
                                   trt*I(X2-BC.ALD$mean.X2)"),
                          data = AC.IPD) {
  fit_stc <- glm(formula, data = data,
                 family = binomial)
  
  # fitted treatment coefficient is relative A vs C conditional effect
  list(mean = coef(fit_stc)["trt"],
       var = vcov(fit_stc)["trt", "trt"])
}


#
IPD_stats.gcomp_ml <- function(strategy,
                               formula =
                                 as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"),
                               data = AC.IPD,
                               R = 1000) {
  
  # non-parametric bootstrap with 1000 resamples
  AC_maic_boot <- boot::boot(data = data,
                             statistic = gcomp_ml.boot,
                             R = R,
                             formula = formula)
  
  list(mean = mean(AC_maic_boot$t),
       var = var(AC_maic_boot$t))
}


#
IPD_stats.gcomp_stan <- function(strategy,
                                 formula =
                                   as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"),
                                 data) {
  
  ppv <- gcomp_stan(formula = formula,
                    dat = data)
  
  # compute marginal log-odds ratio for A vs C for each MCMC sample
  # by transforming from probability to linear predictor scale
  hat.delta.AC <-
    qlogis(rowMeans(ppv$y.star.A)) - qlogis(rowMeans(ppv$y.star.C))
  
  list(mean = mean(hat.delta.AC),
       var = var(hat.delta.AC))
} 



# create class for each approach

#' strategy_maic
#'
#' @param formula 
#' @param R 
#' @param dat_ALD 
#'
#' @return
#' @export
#'
strategy_maic <- function(formula = as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"),
                          R = 1000,
                          dat_ALD = BC.ALD) {
  default_args <- formals()
  args <- as.list(match.call())[-1]
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "maic", args))
}

#' strategy_stc
#'
#' @param formula 
#'
#' @return
#' @export
#'
strategy_stc <- function(formula =
                           as.formula("y ~ X3 + X4 +
                                   trt*I(X1 - mean(X1)) +
                                   trt*I(X2 - mean(X2))")) {
  default_args <- formals()
  args <- as.list(match.call())[-1]
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "stc", args))
}


#' strategy_gcomp_ml
#'
#' @param formula 
#' @param R 
#'
#' @return
#' @export
#'
strategy_gcomp_ml <- function(formula =
                                as.formula("y ~ X3 + X4 + trt*X1 + trt*X2"),
                              R = 1000) {
  default_args <- formals()
  args <- as.list(match.call())[-1]
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "gcomp_ml", args))
}


#' strategy_gcomp_stan
#'
#' @param formula 
#'
#' @return
#' @export
#'
strategy_gcomp_stan <- function(formula =
                                  as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")) {
  default_args <- formals()
  args <- as.list(match.call())[-1]
  args <- modifyList(default_args, args)
  do.call(new_strategy, c(strategy = "gcomp_stan", args))
}


#' new_strategy
#'
#' @param strategy 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
new_strategy <- function(strategy, ...) {
  structure(list(...), class = strategy)
}


#' Main wrapper for hat_Delta_stats
#' 
#' @param AC.IPD 
#' @param BC.ALD 
#' @param strategy 
#' @param ... Additional arguments
#' @export
#' 
hat_Delta_stats <- function(AC.IPD, BC.ALD, strategy, ...) {

  AC_hat_Delta_stats <- IPD_stats(strategy, data = AC.IPD, ...) 
  BC_hat_Delta_stats <- ALD_stats(data = BC.ALD) 
  
  ci_range <- c(0.025, 0.975)
  
  contrasts <- list(
    AB = AC_hat_Delta_stats$mean - BC_hat_Delta_stats$mean,
    AC = AC_hat_Delta_stats$mean,
    BC = BC_hat_Delta_stats$mean)
  
  contrast_variances <- list(
    AB = AC_hat_Delta_stats$var + BC_hat_Delta_stats$var,
    AC = AC_hat_Delta_stats$var,
    BC = BC_hat_Delta_stats$var)
  
  contrast_ci <- list(
    AB = contrasts$AB + qnorm(ci_range)*sqrt(contrast_variances$AB),
    AC = contrasts$AC + qnorm(ci_range)*sqrt(contrast_variances$AC),
    BC = contrasts$BC + qnorm(ci_range)*sqrt(contrast_variances$BC))
  
  stats <- list(contrasts,
                contrast_variances,
                contrast_ci)
  
  structure(stats, class = c("mimR", class(stats)))
}


#' @name IPD_stats
#' Individual level data statistics
#' @export
#' 
IPD_stats <- function(strategy, data, ...)
  UseMethod("IPD_stats", strategy)



#' @rdname IPD_stats
#' 
IPD_stats.default <- function() {
  stop("strategy not available.")
}


#' @rdname IPD_stats
#' IPD_stats.maic
#' 
#' marginal A vs C treatment effect estimates
#' using bootstrapping
#'
#' @param strategy 
#' @param data 
#'
#' @export
#' 
IPD_stats.maic <- function(strategy,
                           data = AC.IPD) {
  # maic.boot(data = data,
  #           indices = 1:nrow(data),
  #           formula = strategy$formula,
  #           dat_ALD = strategy$dat_ALD)
  
  maic_boot <- boot::boot(data = data,
                          statistic = maic.boot,
                          R = strategy$R,
                          formula = strategy$formula,
                          dat_ALD = strategy$dat_ALD)
  
  list(mean =  mean(maic_boot$t),
       var = var(maic_boot$t))
}


#' @rdname IPD_stats
#' IPD_stats.stc
#' 
#' @param strategy 
#' @param data 
#' @export
#' 
IPD_stats.stc <- function(strategy,
                          data = AC.IPD) {
  fit <- glm(strategy$formula,
             data = data,
             family = binomial)
  
  treat_nm <- get_treatment_name(strategy$formula)
  
  # fitted treatment coefficient is relative A vs C conditional effect
  list(mean = coef(fit)[treat_nm],
       var = vcov(fit)[treat_nm, treat_nm])
}


#' @rdname IPD_stats
#' IPD_stats.gcomp_ml
#'
#' @param strategy 
#' @param data 
#'
#' @return
#' @export
#'
IPD_stats.gcomp_ml <- function(strategy,
                               data = AC.IPD) {

  # non-parametric bootstrap with 1000 resamples
  AC_maic_boot <- boot::boot(data = data,
                             statistic = gcomp_ml.boot,
                             R = strategy$R,
                             formula = strategy$formula)
  
  list(mean = mean(AC_maic_boot$t),
       var = var(AC_maic_boot$t))
}


#' @rdname IPD_stats
#' IPD_stats.gcomp_stan
#'
#' @param strategy 
#' @param data 
#'
#' @return
#' @export
#'
IPD_stats.gcomp_stan <- function(strategy,
                                 data) {
  
  ppv <- gcomp_stan(formula = strategy$formula,
                    dat = data)
  
  # compute marginal log-odds ratio for A vs C for each MCMC sample
  # by transforming from probability to linear predictor scale
  hat.delta.AC <-
    qlogis(rowMeans(ppv$y.star.A)) - qlogis(rowMeans(ppv$y.star.C))
  
  list(mean = mean(hat.delta.AC),
       var = var(hat.delta.AC))
} 


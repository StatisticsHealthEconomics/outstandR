
#' @title Calculate the difference between treatments using all evidence
#' 
#' @description
#' This is the main, top-level wrapper for `{mimR}`.
#' Methods taken from
#' \insertCite{RemiroAzocar2022}{mimR}.
#' 
#' @param AC.IPD Individual-level patient data. Suppose between studies _A_ and _C_.
#' @param BC.ALD Aggregate-level data. Suppose between studies _B_ and _C_. 
#' @param strategy Computation strategy function. These can be
#'    `strategy_maic()`, `strategy_stc()`, `strategy_gcomp_ml()` and  `strategy_gcomp_stan()`
#' @param CI Confidence interval; between 0,1
#' @param ... Additional arguments
#' @return List of length 3 of statistics as a `mimR` class object.
#'   Containing statistics between each pair of treatments.
#'   These are the mean contrasts, variances and confidence intervals,
#'   respectively.
#' @importFrom Rdpack reprompt
#' 
#' @references
#' \insertRef{RemiroAzocar2022}{mimR}
#' 
#' @export
#' @examples
#' data(AC_IPD)  # AC patient-level data
#' data(BC_ALD)  # BC aggregate-level data
#' 
#' lin_form <- as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")
#' 
#' # matching-adjusted indirect comparison
#' mimR_maic <- mimR(AC_IPD, BC_ALD, strategy = strategy_maic(lin_form))
#' 
#' # simulated treatment comparison
#' mimR_stc <- mimR(AC_IPD, BC_ALD, strategy = strategy_stc(lin_form))
#' 
#' # G-computation with maximum likelihood
#' # mimR_gcomp_ml <- mimR(AC_IPD, BC_ALD, strategy = strategy_gcomp_ml(lin_form))
#' 
#' # G-computation with Bayesian inference
#' mimR_gcomp_stan <- mimR(AC_IPD, BC_ALD, strategy = strategy_gcomp_stan(lin_form))
#' 
mimR <- function(AC.IPD, BC.ALD, strategy, CI = 0.95, ...) {
  
  if (CI <= 0 || CI >= 1) stop("CI argument must be between 0 and 1.")
  ##TODO: as method instead?
  if (!inherits(strategy, "strategy"))
    stop("strategy argument must be a class strategy.")
  
  # select data according to formula
  ipd <- model.frame(strategy$formula, data = AC.IPD)
  
  term.labels <- attr(terms(strategy$formula), "term.labels")
  mean_names <- paste0("mean.", term.labels)
  sd_names <- paste0("sd.", term.labels)
  term_names <- c(mean_names, sd_names)
  
  # remove treatment labels
  term_names <- sort(term_names[!grepl(pattern = "trt", term_names)])
  
  # replace outcome variable name
  response_var <- all.vars(strategy$formula)[1]
  response_names <- gsub(pattern = "y", replacement = response_var,
                         x = c("y.B.sum", "y.B.bar", "N.B", "y.C.sum", "y.C.bar", "N.C")) 
  
  keep_names <- c(term_names, response_names)
  
  ald <- BC.ALD[keep_names]
  
  AC_mimR <- IPD_stats(strategy, ipd = ipd, ald = ald, ...) 
  BC_mimR <- ALD_stats(ald = ald) 
  
  upper <- 0.5 + CI/2
  ci_range <- c(1-upper, upper)
  
  contrasts <- list(
    AB = AC_mimR$mean - BC_mimR$mean,
    AC = AC_mimR$mean,
    BC = BC_mimR$mean)
  
  contrast_variances <- list(
    AB = AC_mimR$var + BC_mimR$var,
    AC = AC_mimR$var,
    BC = BC_mimR$var)
  
  contrast_ci <- list(
    AB = contrasts$AB + qnorm(ci_range)*as.vector(sqrt(contrast_variances$AB)),
    AC = contrasts$AC + qnorm(ci_range)*as.vector(sqrt(contrast_variances$AC)),
    BC = contrasts$BC + qnorm(ci_range)*as.vector(sqrt(contrast_variances$BC)))
  
  stats <- list(contrasts = contrasts,
                variances = contrast_variances,
                CI = contrast_ci)
  
  structure(stats,
            CI = CI,
            class = c("mimR", class(stats)))
}

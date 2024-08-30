
#' @title Calculate the difference between treatments using all evidence
#' 
#' @description
#' This is the main, top-level wrapper for `{outstandR}`.
#' Methods taken from
#' \insertCite{RemiroAzocar2022}{outstandR}.
#' 
#' @param AC.IPD Individual-level patient data. Suppose between studies _A_ and _C_.
#' @param BC.ALD Aggregate-level data. Suppose between studies _B_ and _C_. 
#' @param strategy Computation strategy function. These can be
#'    `strategy_maic()`, `strategy_stc()`, `strategy_gcomp_ml()` and  `strategy_gcomp_stan()`
#' @param CI Confidence interval; between 0,1
#' @param ... Additional arguments
#' @return List of length 3 of statistics as a `outstandR` class object.
#'   Containing statistics between each pair of treatments.
#'   These are the mean contrasts, variances and confidence intervals,
#'   respectively.
#' @importFrom Rdpack reprompt
#' 
#' @references
#' \insertRef{RemiroAzocar2022}{outstandR}
#' 
#' @export
#' @examples
#' data(AC_IPD)  # AC patient-level data
#' data(BC_ALD)  # BC aggregate-level data
#' 
#' lin_form <- as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")
#' 
#' # matching-adjusted indirect comparison
#' outstandR_maic <- outstandR(AC_IPD, BC_ALD, strategy = strategy_maic(formula = lin_form))
#' 
#' # simulated treatment comparison
#' outstandR_stc <- outstandR(AC_IPD, BC_ALD, strategy = strategy_stc(lin_form))
#' 
#' # G-computation with maximum likelihood
#' # outstandR_gcomp_ml <- outstandR(AC_IPD, BC_ALD, strategy = strategy_gcomp_ml(lin_form))
#' 
#' # G-computation with Bayesian inference
#' outstandR_gcomp_stan <- outstandR(AC_IPD, BC_ALD, strategy = strategy_gcomp_stan(lin_form))
#' 
#' # Multiple imputation marginalization
#' outstandR_gcomp_stan <- outstandR(AC_IPD, BC_ALD, strategy = strategy_mim(lin_form))
#' 
outstandR <- function(AC.IPD, BC.ALD, strategy, CI = 0.95, ...) {
  
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

  AC_outstandR <- IPD_stats(strategy, ipd = ipd, ald = ald, ...) 
  BC_outstandR <- ALD_stats(ald = ald) 
  
  upper <- 0.5 + CI/2
  ci_range <- c(1-upper, upper)
  
  contrasts <- list(
    AB = AC_outstandR$mean - BC_outstandR$mean,
    AC = AC_outstandR$mean,
    BC = BC_outstandR$mean)
  
  contrast_variances <- list(
    AB = AC_outstandR$var + BC_outstandR$var,
    AC = AC_outstandR$var,
    BC = BC_outstandR$var)
  
  contrast_ci <- list(
    AB = contrasts$AB + qnorm(ci_range)*as.vector(sqrt(contrast_variances$AB)),
    AC = contrasts$AC + qnorm(ci_range)*as.vector(sqrt(contrast_variances$AC)),
    BC = contrasts$BC + qnorm(ci_range)*as.vector(sqrt(contrast_variances$BC)))
  
  stats <- list(contrasts = contrasts,
                variances = contrast_variances,
                CI = contrast_ci)
  
  structure(stats,
            CI = CI,
            class = c("outstandR", class(stats)))
}

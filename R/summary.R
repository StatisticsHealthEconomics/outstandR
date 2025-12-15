
#' Summary method for outstandR
#'  
#' @param object [outstandR()] output object
#' @param ... Additional arguments
#' @return List of class `summary.outstandR`
#' 
#' @export
summary.outstandR <- function(object, CI = NA, ...) {
  if(!inherits(object, "outstandR")) stop("Object supplied is not of class 'outstandR'", call. = FALSE)
  
  res <- object
  
  # lazy evaluation confidence interval
  if (!is.na(CI)) {
    res$contrast$CI <- list(
      AB = calc_ci(mean_val = res$contrasts$AB, sd_val = sqrt(res$contrast_variances$AB), level = CI),
      AC = calc_ci(mean_val = res$contrasts$AC, sd_val = sqrt(res$contrast_variances$AC), level = CI),
      BC = calc_ci(mean_val = res$contrasts$BC, sd_val = sqrt(res$contrast_variances$BC), level = CI)
    )
    ##TODO:
    res$absolute$CI <- list(
      AB = NA,
      AC = NA,
      BC = NA
    )
    
    attr(res, "CI") <- CI
  }
  
  class(res) <- "summary.outstandR"
  return(res)
}

#' @rdname summary.outstandR
#' @return Original argument, but mainly called for side effects
#' @export
#' @method print summary.outstandR
print.summary.outstandR <- function(x, digits = 3, ...) {
  
  # header Information
  cat("\nSummary\n")
  cat("---------------------\n")
  cat("Formula: ", deparse(x$formula), "\n") # deparse handles complex formulas better
  cat("Family:  ", x$family, paste0("(", x$link, ")"), "\n")
  cat("Contrast:", x$contrast, "\n")
  
  # main estimates
  cat("\nParameter Estimates:\n")
  cat(paste0("\nParameter Estimates (CI Level: ", ci_label, "):\n"))
  
  # model diagnostics
  cat("\nUnderlying Model Fit:\n")
  if(!is.null(x$aic)) {
    cat("AIC:", format(x$aic, digits = digits), "\n")
  } else {
    cat("AIC: Not available\n")
  }
  
  invisible(x)
}

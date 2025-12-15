
#' Summary method for outstandR
#'  
#' @param object [outstandR()] output object
#' @param ... Additional arguments
#' @return No return value, called for side effects
#' 
#' @export
summary.outstandR <- function(object, CI = NA, ...) {
  if(!inherits(object, "outstandR")) stop("Object supplied is not of class 'outstandR'", call. = FALSE)
  
  res <- object
  
  # lazy evaluation confidence interval
  if (!is.na(CI)) {
    ##TODO:
    # res$coefficients <- calc_ci(object$boot_dist, level = CI)
    attr(res, "ci_level") <- CI 
  }
  
  class(res) <- "summary.outstandR"
  return(res)
}

#' @rdname summary.outstandR
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

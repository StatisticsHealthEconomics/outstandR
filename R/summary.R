#' Summary method for outstandR
#' 
#' @param object [outstandR()] output object.
#' @param CI Confidence interval level.
#' @param ... Additional arguments.
#' @return List of class `summary.outstandR`
#' 
#' @export
summary.outstandR <- function(object, CI = NA, ...) {
  if(!inherits(object, "outstandR")) stop("Object supplied is not of class 'outstandR'", call. = FALSE)
  
  res <- object
  
  # lazy evaluation confidence interval
  if (!is.na(CI)) {
    # Check if 'contrasts' (plural) exists in object, user code referenced both singular/plural
    # Assuming standard object structure here:
    res$contrast$CI <- list(
      AB = calc_ci(mean_val = res$contrasts$AB, 
                   sd_val = sqrt(res$contrast_variances$AB), 
                   level = CI),
      AC = calc_ci(mean_val = res$contrasts$AC, 
                   sd_val = sqrt(res$contrast_variances$AC), 
                   level = CI),
      BC = calc_ci(mean_val = res$contrasts$BC, 
                   sd_val = sqrt(res$contrast_variances$BC), 
                   level = CI)
    )
    
    res$absolute$CI <- lapply(names(res$absolute), function(grp_name) {
      calc_ci(
        mean_val = res$absolute[[grp_name]], 
        sd_val   = sqrt(res$absolute_variances[[grp_name]]), 
        level    = CI)
    })
    
    names(res$absolute$CI) <- names(res$absolute)
  }
  
  attr(res, "CI") <- CI
  
  class(res) <- "summary.outstandR"
  return(res)
}

#' @rdname summary.outstandR
#' @param x An object used to select a method.
#' @param digits Minimal number of significant digits, see `print.default`.
#' @return Original argument, but mainly called for side effects
#' @export
#' @method print summary.outstandR
print.summary.outstandR <- function(x, digits = 3, ...) {
  
  # header Information
  cat("\nSummary\n")
  cat("---------------------\n")
  cat("Formula: ", deparse(x$formula), "\n") 
  cat("Family:  ", x$family, paste0("(", x$link, ")"), "\n")
  cat("Contrast:", x$contrasts, "\n")
  
  ci_val <- attr(x, "CI")
  
  # Format the label based on whether CI exists
  if (!is.null(ci_val) && !is.na(ci_val)) {
    ci_label <- paste0(ci_val * 100, "%")
  } else {
    ci_label <- "None"
  }
  
  # Main estimates (Consolidated duplicate headers)
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

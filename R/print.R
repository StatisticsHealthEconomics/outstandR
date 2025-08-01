
#' Print a Summary of a outstandR Object
#' 
#' This is a method for the function [print()] for objects of the
#' class "outstandR" created by a call to [outstandR()]
#' 
#' @param x Objects of the class "outstandR"
#' @param ... Additional arguments passed to other methods
#' @seealso [outstandR()]
#' @import pillar
#' @importFrom tibble tibble
#' @export
#' 
print.outstandR <- function(x, ...) {
  if (!requireNamespace("tibble", quietly = TRUE) || !requireNamespace("pillar", quietly = TRUE)) {
    stop("Please install the 'tibble' and 'pillar' packages for colored tibble output.")
  }
  
  ref_trt <- attr(x, "ref_trt")
  
  cat(pillar::style_bold("Object of class 'outstandR'"), "\n")
  cat("Model:", pillar::style_subtle(attr(x, "model")), "\n")
  cat("Scale:", pillar::style_subtle(attr(x, "scale")), "\n")
  cat("Common treatment:", pillar::style_subtle(ref_trt), "\n")
  # cat(pillar::style_subtle("Common treatment:"), attr(x, "reference"), "\n")  ##TODO:
  cat("Individual patient data study:", pillar::style_subtle("AC"), "\n")
  cat("Aggregate level data study:", pillar::style_subtle("BC"), "\n")
  cat("Confidence interval level:", pillar::style_subtle(attr(x, "CI")), "\n\n")
  
  # Function to color CI values
  color_ci <- function(value) {
    if (is.na(value)) pillar::style_subtle("NA")
    else if (value > 0) pillar::style_success(sprintf("%.3f", value))
    else if (value < 0) pillar::style_warning(sprintf("%.3f", value))
    else sprintf("%.3f", value)
  }
  
  contrasts <- x$contrasts
  absolute <- x$absolute
  
  con_tab <- tibble::tibble(
    Treatments = names(contrasts$means),
    Estimate = unlist(contrasts$means),
    Std.Error = unlist(contrasts$variances),
    lower.0.95 = sapply(contrasts$CI, \(x) x[1]),
    upper.0.95 = sapply(contrasts$CI, \(x) x[2])
  )
  
  abs_tab <- tibble::tibble(
    Treatments = names(absolute$means),
    Estimate = unlist(absolute$means),
    Std.Error = unlist(absolute$variances),
    lower.0.95 = c(NA,NA),  #sapply(absolute$CI, \(x) x[1]),  ##TODO:
    upper.0.95 = c(NA,NA)   #sapply(absolute$CI, \(x) x[2])
  )
  
  cat("Contrasts:\n\n")
  print(con_tab, n = Inf)

  cat("\nAbsolute:\n\n")
  print(abs_tab, n = Inf)
  
  invisible()
}


#' @export
#' 
print.outstandR <- function(x, ...) {
  if (!requireNamespace("tibble", quietly = TRUE) || !requireNamespace("pillar", quietly = TRUE)) {
    stop("Please install the 'tibble' and 'pillar' packages for colored tibble output.")
  }
  
  cat(pillar::style_bold("Object of class 'outstandR'"), "\n")
  cat("Scale:", pillar::style_subtle(attr(x, "scale")), "\n")
  cat("Common treatment:", pillar::style_subtle("C"), "\n")
  # cat(pillar::style_subtle("Common treatment:"), attr(x, "reference"), "\n")
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
  
  tib <- tibble::tibble(
    Treatments = names(x$contrasts),
    Estimate = unlist(x$contrasts),
    Std.Error = unlist(x$variances),
    lower.0.95 = sapply(x$CI, \(x) x[1]),
    upper.0.95 = sapply(x$CI, \(x) x[2])
  )
  
  print(tib, n = Inf)
  
  invisible(tib)
}

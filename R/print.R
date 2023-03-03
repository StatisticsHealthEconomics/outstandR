
#'
print.mimR <- function(x, newline = TRUE) {
  
  cat("\nContrasts:", x$contrasts, "\n")
  cat("\nVariances:", x$contrasts_variances, "\n")
  cat("\nConfidence intervals:", x$contrasts_ci, "\n")
  
  if (newline) {
    cat("\n")
  }
  
  invisible(x)
}


#' Summary method for outstandR
#'  
#' @param object [outstandR()] output object
#' @param ... Additional arguments
#' @return No return value, called for side effects
#' 
#' @export
summary.outstandR <- function(object, ...) {
  
  ans <- NULL
  class(ans) <- "summary.outstandR"
  ans
} 

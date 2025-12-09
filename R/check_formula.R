
#' Check formula
#' 
#' @param formula Formula
#' @param trt_var Treatment variable
#' @return No return value, called for side effects
#' @keywords internal
#' 
check_formula <- function(formula, trt_var = NULL) {
  
  if (!inherits(formula, "formula"))
    stop("`formula` must be of class 'formula'.", call. = FALSE)

  terms_labels <- attr(terms(formula), "term.labels")
  
  if (!is.null(trt_var)) {
    if (!(trt_var %in% terms_labels))
      stop(sprintf("Treatment term '%s' is missing in the formula", trt_var), call. = FALSE)
  }
  
  invisible()
}

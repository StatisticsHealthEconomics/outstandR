
#' @keywords internal
#' 
check_formula <- function(formula, treatment = "trt") {
  
  if (!inherits(formula, "formula"))
    stop("`formula` must be of class 'formula'.")

  terms_labels <- attr(terms(formula), "term.labels")
  
  if (!(treatment %in% terms_labels))
    stop(sprintf("Treatment term '%s' is missing in the formula", treatment))
  
  invisible()
}

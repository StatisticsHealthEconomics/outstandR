
#
check_formula <- function(formula) {
  
  if (!inherit(formula, "formula"))
    stop("formula argument must be of formula class.")
  
  if (!"trt" %in% attr(terms(formula), "term.labels"))
    stop("Treatment term, trt, is missing in the formula")
  
  invisible()
}

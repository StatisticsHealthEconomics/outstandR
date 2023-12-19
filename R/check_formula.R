
#
check_formula <- function(formula) {
  
  if (!"trt" %in% attr(terms(formula), "term.labels"))
    stop("Treatment term, trt, is missing in the formula")
  
  invisible()
}

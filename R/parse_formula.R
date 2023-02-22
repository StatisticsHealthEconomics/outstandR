
#
get_effect_modifiers <- function(formula) {
  
  formula <- as.formula(formula)# 
  term.labels <- attr(terms(formula), "term.labels")
  
  modifiers <- term.labels[grepl(":", term.labels)]
  modifiers <- gsub(".+:", "", modifiers)
  
  modifiers
}


# Extract RHS of a formula
#
# @param x A formula object.
# @return An expression.
rhs <- function(x, as_formula = FALSE) {
  len <- length(x)
  if (len == 3L) {
    out <- x[[3L]]
  } else {
    out <- x[[2L]]
  }
  out
}

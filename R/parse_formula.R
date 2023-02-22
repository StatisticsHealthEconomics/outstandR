
#
get_effect_modifiers <- function(formula) {
  
  formula <- as.formula(formula)
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


#
get_treatment_name <- function(formula) {
  formula <- as.formula(formula)
  term.labels <- attr(terms(formula), "term.labels")
  
  treat_nm <- term.labels[grepl(":", term.labels)]
  treat_nm <- gsub(":.+", "", treat_nm[1])
  treat_nm
}


#
get_mean_names <- function(formula, dat) {
  
}

#
get_sd_names <- function(formula, dat) {
  
}

#
get_covariate_names <- function(formula) {

}

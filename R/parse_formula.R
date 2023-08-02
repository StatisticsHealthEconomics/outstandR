
#' Get effect modifiers
#'
#' @param formula 
#'
#' @return effect modifiers names
#' @export
#'
get_effect_modifiers <- function(formula) {
  
  formula <- as.formula(formula)
  term.labels <- attr(terms(formula), "term.labels")
  
  modifiers <- term.labels[grepl(":", term.labels)]
  modifiers <- gsub(".+:", "", modifiers)
  
  modifiers
}

#' Get treatment name
#'
#' @param formula 
#'
#' @return treatment name
#' @export
#'
get_treatment_name <- function(formula) {
  formula <- as.formula(formula)
  term.labels <- attr(terms(formula), "term.labels")
  
  treat_nm <- term.labels[grepl(":", term.labels)]
  treat_nm <- gsub(":.+", "", treat_nm[1])
  treat_nm
}

#' Get mean names
#'
#' @param dat 
#' @param var_names 
#'
#' @return mean names
#' @export
#'
get_mean_names <- function(dat, var_names) {
  dat_names <- names(dat)
  # is_sd_name <- grepl(pattern = "\\.mean", dat_names)
  is_mean_name <- grepl(pattern = "mean\\.", dat_names)
  is_var_name <- grepl(pattern = paste(var_names, collapse = "|"), dat_names)
  
  dat_names[is_mean_name & is_var_name]
}

#' Get SD names
#'
#' @param dat 
#' @param var_names 
#'
#' @return SD names
#' @export
#'
get_sd_names <- function(dat, var_names) {
  dat_names <- names(dat)
  # is_sd_name <- grepl(pattern = "\\.sd", dat_names)
  is_sd_name <- grepl(pattern = "sd\\.", dat_names)
  is_var_name <- grepl(pattern = paste(var_names, collapse = "|"), dat_names)
  
  dat_names[is_sd_name & is_var_name]
}

#' Get covariate names
#'
#' @param formula 
#'
#' @return covariate names
#' @export
#'
get_covariate_names <- function(formula) {
  all.vars(formula)[-1]
}

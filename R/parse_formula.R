
#' Get effect modifiers
#'
#' @param formula Linear regression formula string
#'
#' @return Effect modifiers names
#' @keywords internal
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
#' @param formula Linear regression formula string
#'
#' @return Treatment name
#' @keywords internal
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
#' @template args-ald
#' @param var_names Variable names character vector
#'
#' @return Mean names vector
#' @keywords internal
#'
get_mean_names <- function(ald, var_names) {
  dat_names <- names(ald)
  # is_sd_name <- grepl(pattern = "\\.mean", dat_names)
  is_mean_name <- grepl(pattern = "mean\\.", dat_names)
  is_var_name <- grepl(pattern = paste(var_names, collapse = "|"), dat_names)
  
  dat_names[is_mean_name & is_var_name]
}

#' Get standard deviation names
#'
#' @template args-ald
#' @param var_names Variable names character vector
#'
#' @return Standard deviation names vector
#' @keywords internal
#'
get_sd_names <- function(ald, var_names) {
  dat_names <- names(ald)
  # is_sd_name <- grepl(pattern = "\\.sd", dat_names)
  is_sd_name <- grepl(pattern = "sd\\.", dat_names)
  is_var_name <- grepl(pattern = paste(var_names, collapse = "|"), dat_names)
  
  dat_names[is_sd_name & is_var_name]
}

#' Get covariate names
#'
#' @param formula Linear regression formula object
#'
#' @return covariate names vector
#' @keywords internal
#'
get_covariate_names <- function(formula) {
  
  if (!inherits(formula, "formula"))
    stop("formula argument must be of formula class.")
  
  all.vars(formula)[-1]
}

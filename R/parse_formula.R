
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
  
  if (is.na(treat_nm))
    warning("Treatment name missing from formula.")
  
  treat_nm
}

#' Get mean names
#'
#' @template args-ald
#' @param keep_nms Variable names character vector
#'
#' @return Mean names vector
#' @keywords internal
#'
get_mean_names <- function(ald, keep_nms) {
  
  dat_names <- names(ald)
  # is_sd_name <- grepl(pattern = "\\.mean", dat_names)
  is_mean_name <- grepl(pattern = "mean\\.", dat_names)
  is_var_name <- grepl(pattern = paste(keep_nms, collapse = "|"), dat_names)
  keep_mean_nm <- is_mean_name & is_var_name
  
  if (all(!keep_mean_nm))
    warning("No matching mean names found.")
  
  dat_names[keep_mean_nm]
}

#' Get standard deviation names
#'
#' @template args-ald
#' @param keep_nms Variable names character vector
#'
#' @return Standard deviation names vector
#' @keywords internal
#'
get_sd_names <- function(ald, keep_nms) {
  
  dat_names <- names(ald)
  # is_sd_name <- grepl(pattern = "\\.sd", dat_names)
  is_sd_name <- grepl(pattern = "sd\\.", dat_names)
  is_var_name <- grepl(pattern = paste(keep_nms, collapse = "|"), dat_names)
  keep_sd_nm <- is_sd_name & is_var_name
  
  if (all(!keep_sd_nm))
    warning("No matching sd names found.")
  
  dat_names[keep_sd_nm]
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

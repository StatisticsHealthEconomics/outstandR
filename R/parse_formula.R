
#' Guess treatment name
#' 
#' Does a variable appear more than once in interactions?
#' Otherwise pick first LHS interaction term.
#' So this doesnt work without an effect modifier.
#' 
#' @eval reg_args(include_formula = TRUE, include_family = FALSE)
#'
#' @return Treatment name
#' @keywords internal
#'
guess_treatment_name <- function(formula) {
  
  formula <- as.formula(formula)
  term.labels <- attr(terms(formula), "term.labels")
  
  interaction_terms <- term.labels[grepl(":", term.labels)]
  
  # Split each interaction term by ":" and flatten
  vars <- unlist(strsplit(interaction_terms, ":"))
  
  # Find the first duplicated variable (assumed to be treatment)
  treat_nm <- vars[duplicated(vars)][1]
  
  # no duplicate
  if (length(treat_nm) == 0) {
    if (length(interaction_terms) > 0) {
      # take the first term
      treat_nm <- strsplit(interaction_terms[1], ":")[[1]][1]
    } else {
      stop("Unable to guess the treatment name from formula.")
    }
  }
  
  treat_nm
}

#
get_treatment_name <- function(formula, trt_var = NULL) {
  if (!is.null(trt_var) && is.character(trt_var)) {
    return(trt_var)
  }
  
  guess_treatment_name(formula)
}

#' Get mean names
#'
#' @eval study_data_args(include_ipd = FALSE, include_ald = TRUE)
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
  
  mean_nms <- dat_names[keep_mean_nm]
   
  covariate_nms <- sub(".*mean\\.", "", mean_nms)
  
  setNames(mean_nms, covariate_nms)
}

#' Get standard deviation names
#'
#' @eval study_data_args(include_ipd = FALSE, include_ald = TRUE)
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
  
  sd_nms <- dat_names[keep_sd_nm]
  
  covariate_nms <- sub(".*sd\\.", "", sd_nms)
  
  setNames(sd_nms, covariate_nms)
}

#' Get covariate names
#'
#' @eval reg_args(include_formula = TRUE, include_family = FALSE)
#' 
#' @return covariate names vector
#' @keywords internal
#'
get_covariate_names <- function(formula) {
  
  if (!inherits(formula, "formula"))
    stop("formula argument must be of formula class.")
  
  all.vars(formula)[-1]
}

#' Get effect modifiers
#'
#' @eval reg_args(include_formula = TRUE, include_family = FALSE)
#' 
#' @return Effect modifiers names
#' @keywords internal
#'
get_eff_mod_names <- function(formula, trt_var) {
  
  formula <- as.formula(formula)
  
  term.labels <- attr(terms(formula), "term.labels")
  
  # effect modifier terms only
  eff_mod_terms <- term.labels[grepl(":", term.labels)]
  
  gsub(paste0("^", trt_var, ":"), "", eff_mod_terms)
}

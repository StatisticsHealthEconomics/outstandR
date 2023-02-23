
#
get_effect_modifiers <- function(formula) {
  
  formula <- as.formula(formula)
  term.labels <- attr(terms(formula), "term.labels")
  
  modifiers <- term.labels[grepl(":", term.labels)]
  modifiers <- gsub(".+:", "", modifiers)
  
  modifiers
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
get_mean_names <- function(dat, var_names) {
  dat_names <- names(dat)
  is_sd_name <- grepl(pattern = "\\.mean", dat_names)
  is_var_name <- grepl(pattern = var_names, dat_names)
  
  dat_names[is_sd_name & is_var_name]
}

#
get_sd_names <- function(dat, var_names) {
  dat_names <- names(dat)
  is_sd_name <- grepl(pattern = "\\.sd", dat_names)
  is_var_name <- grepl(pattern = var_names, dat_names)
  
  dat_names[is_sd_name & is_var_name]
}

#
get_covariate_names <- function(formula) {
  all.vars(formula)[-1]
}


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
get_mean_names <- function(formula, dat) {
  dat_names <- names(dat)
  cov_names <- get_covariate_names(formula)
  is_sd_name <- grepl(pattern = "\\.mean", dat_names)
  is_cov_name <- grepl(pattern = cov_names, dat_names)
  
  dat_names[is_sd_name & is_cov_name]
}

#
get_sd_names <- function(formula, dat) {
  dat_names <- names(dat)
  cov_names <- get_covariate_names(formula)
  is_sd_name <- grepl(pattern = "\\.sd", dat_names)
  is_cov_name <- grepl(pattern = cov_names, dat_names)
  
  dat_names[is_sd_name & is_cov_name]
}

#
get_covariate_names <- function(formula) {
  all.vars(formula)[-1]
}

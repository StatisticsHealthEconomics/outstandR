
#' Guess treatment name
#' 
#' Does a variable appear more than once in interactions?
#' If not then pick first LHS interaction term.
#' Finally, if there are no interactions then pick last main effect term.
#' 
#' @eval reg_args(include_formula = TRUE, include_family = FALSE)
#' @return Treatment name string.
#' 
#' @importFrom crayon yellow
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
  if (is.na(treat_nm) || length(treat_nm) == 0) {
    if (length(interaction_terms) > 0) {
      # take the first term
      treat_nm <- strsplit(interaction_terms[1], ":")[[1]][1]
    } else {
      # Identify main effect terms (i.e., no colon)
      main_effect_terms <- term.labels[!grepl(":", term.labels)]

      # Pick the last main effect term
      treat_nm <- tail(main_effect_terms, 1)
      
      if (length(treat_nm) == 0) {
        stop("Treatment term 'trt' is missing in the formula")
      }
      
      message("Treatment is guessed as: ", crayon::yellow(treat_nm))
    }
  }
  
  treat_nm
}

#' Get treatment name
#'
#' @param formula Formula
#' @param trt_var Treatment variable
#' @return Treatment name string.
#' @keywords internal
#' 
get_treatment_name <- function(formula, trt_var = NULL) {
  if (!is.null(trt_var) && is.character(trt_var)) {
    return(trt_var)
  }
  
  guess_treatment_name(formula)
}

#' Get covariate names
#'
#' @eval reg_args(include_formula = TRUE, include_family = FALSE)
#' 
#' @return Covariate names character vector
#' @keywords internal
#'
get_covariate_names <- function(formula) {
  
  if (!inherits(formula, "formula")) {
    stop("formula argument must be of formula class.", call. = FALSE)
  }
  
  all.vars(formula)[-1]
}

#' Get effect modifiers
#'
#' @param trt_var Treatment variable name; Default 'trt'.
#' @eval reg_args(include_formula = TRUE, include_family = FALSE)
#' 
#' @return Effect modifiers strings names
#' @keywords internal
#'
get_eff_mod_names <- function(formula, trt_var = "trt") {
  
  formula <- as.formula(formula)
  
  all_terms <- attr(terms(formula), "term.labels")
  
  effect_modifiers <- character()
  
  for (term in all_terms) {
    is_interaction <- grepl(":", term)
    
    if (is_interaction) {
      components <- strsplit(term, ":")[[1]]
      
      if (trt_var %in% components) {
        other_components <- components[components != trt_var]
        effect_modifiers <- c(effect_modifiers, other_components)
      }
    }
  }
  
  unique(effect_modifiers)
}

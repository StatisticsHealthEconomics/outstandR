# -------------------------------------------------------------------------
# Variance Method Validation
# -------------------------------------------------------------------------

#' Determine and validate variance method for a strategy
#' 
#' @param strategy Strategy object
#' @param user_method User supplied method (string) or NULL
#' @return String name of variance method to use
#' @keywords internal
get_var_method <- function(strategy, user_method = NULL) {
  
  allowed <- get_allowed_var_methods(strategy)
  
  if (length(allowed) == 0) {
    stop("No variance methods defined for this strategy.", call. = FALSE)
  }
  
  # Default to the first method in the allowed list
  if (is.null(user_method)) {
    return(allowed[1])
  }
  
  if (!user_method %in% allowed) {
    stop(paste("Variance method for", class(strategy)[1], "must be one of:", 
               paste(allowed, collapse = ", ")), call. = FALSE)
  }
  
  user_method
}

#' Retrieve list of allowed variance methods
#' 
#' @param strategy Strategy object
#' @return Character vector of allowed methods. First element is the default.
#' @keywords internal
get_allowed_var_methods <- function(strategy) UseMethod("get_allowed_var_methods", strategy)

#' @keywords internal
get_allowed_var_methods.default <- function(strategy) {
  c("sample", "sandwich")
}

#' @keywords internal
get_allowed_var_methods.maic <- function(strategy) {
  c("sample", "sandwich")
}

#' @keywords internal
get_allowed_var_methods.stc <- function(strategy) {
  c("sample", "sandwich")
}

#' @keywords internal
get_allowed_var_methods.gcomp_ml <- function(strategy) {
  c("sample", "sandwich")
}

#' @keywords internal
get_allowed_var_methods.gcomp_bayes <- function(strategy) {
  c("sample", "sandwich")
}

#' @keywords internal
get_allowed_var_methods.mim <- function(strategy) {
  c("rubin")
}

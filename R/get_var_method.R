# -------------------------------------------------------------------------
# Variance Method Validation (S3 Generics)
# -------------------------------------------------------------------------

#' Determine and validate variance method for a strategy
#' 
#' @param strategy Strategy object
#' @param user_method User supplied method (string) or NULL
#' @return String name of variance method
#' @keywords internal
get_var_method <- function(strategy, user_method = NULL) UseMethod("get_var_method", strategy)

#' @keywords internal
get_var_method.default <- function(strategy, user_method = NULL) {
  # Fallback allowed methods
  allowed <- c("sample", "sandwich")
  
  if (is.null(user_method)) {
    return("sample")
  }
  
  if (!user_method %in% allowed) {
    stop(paste("Variance method must be one of:", paste(allowed, collapse = ", ")), call. = FALSE)
  }
  
  user_method
}

#' @keywords internal
get_var_method.maic <- function(strategy, user_method = NULL) {
  allowed <- c("sample", "sandwich")
  
  if (is.null(user_method)) {
    return("sample")
  }
  
  if (!user_method %in% allowed) {
    stop(paste("MAIC variance method must be one of:", paste(allowed, collapse = ", ")), call. = FALSE)
  }
  
  user_method
}

#' @keywords internal
get_var_method.stc <- function(strategy, user_method = NULL) {
  allowed <- c("sample", "sandwich")
  
  if (is.null(user_method)) {
    return("sample")
  }
  
  if (!user_method %in% allowed) {
    stop(paste("STC variance method must be one of:", paste(allowed, collapse = ", ")), call. = FALSE)
  }
  
  user_method
}

#' @keywords internal
get_var_method.gcomp_ml <- function(strategy, user_method = NULL) {
  allowed <- c("sample", "sandwich")
  
  if (is.null(user_method)) {
    return("sample")
  }
  
  if (!user_method %in% allowed) {
    stop(paste("G-computation (ML) variance method must be one of:", paste(allowed, collapse = ", ")), call. = FALSE)
  }
  
  user_method
}

#' @keywords internal
get_var_method.gcomp_bayes <- function(strategy, user_method = NULL) {
  # Typically sample-based from posterior, but sandwich allows checking frequentist properties if implemented
  allowed <- c("sample", "sandwich")
  
  if (is.null(user_method)) {
    return("sample")
  }
  
  if (!user_method %in% allowed) {
    stop(paste("G-computation (Bayes) variance method must be one of:", paste(allowed, collapse = ", ")), call. = FALSE)
  }
  
  user_method
}

#' @keywords internal
get_var_method.mim <- function(strategy, user_method = NULL) {
  allowed <- c("rubin")
  
  if (is.null(user_method)) {
    return("rubin")
  }
  
  if (!user_method %in% allowed) {
    stop(paste("MIM variance method must be one of:", paste(allowed, collapse = ", ")), call. = FALSE)
  }
  
  user_method
}

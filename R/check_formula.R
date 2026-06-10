#' Check outcome formula
#' 
#' @param formula Formula for the outcome model
#' @param trt_var Treatment variable
#' @return No return value, called for side effects
#' @keywords internal
#' @export
check_formula <- function(formula, trt_var = NULL) {
  check_formula_base(
    formula, 
    trt_var = trt_var, 
    must_have_response = TRUE, 
    must_have_trt = TRUE
  )
}

#' Check balance formula
#' 
#' @param formula Formula for the balance model
#' @param trt_var Treatment variable (used to warn if mistakenly included)
#' @return No return value, called for side effects
#' @keywords internal
#' @export
check_balance_formula <- function(formula, trt_var = NULL) {
  check_formula_base(
    formula, 
    trt_var = trt_var, 
    must_have_response = FALSE, 
    must_have_trt = FALSE
  )
}

#' Base formula validator
#' @keywords internal
check_formula_base <- function(formula, 
                               trt_var = NULL, 
                               must_have_response = TRUE, 
                               must_have_trt = TRUE) {
  
  if (!inherits(formula, "formula")) {
    stop("`formula` must be of class 'formula'.", call. = FALSE)
  }
  
  # R interprets a formula as having a response if attr(terms, "response") == 1
  has_response <- attr(terms(formula), "response") == 1
  terms_labels <- attr(terms(formula), "term.labels")
  
  # --- 1. Response Variable Checks ---
  if (must_have_response && !has_response) {
    stop("Outcome formula must have a response variable (e.g., y ~ ...).", call. = FALSE)
  }
  
  if (!must_have_response && has_response) {
    stop("Balance formula should be one-sided (e.g., ~ X1 + X2) and not contain a response variable.", call. = FALSE)
  }
  
  # --- 2. Treatment Variable Checks ---
  if (!is.null(trt_var)) {
    has_trt <- trt_var %in% terms_labels
    
    if (must_have_trt && !has_trt) {
      stop(sprintf("Treatment term '%s' is missing in the outcome formula.", trt_var), call. = FALSE)
    }
    
    if (!must_have_trt && has_trt) {
      warning(sprintf("Treatment term '%s' was found in the balance formula. Balance models should typically only contain covariates.", trt_var), call. = FALSE)
    }
  }
  
  # --- 3. Survival Object Checks ---
  if (has_response) {
    response <- formula[[2]]
    # Check if evaluated response inherits "Surv" OR if the formula explicitly calls Surv()
    is_surv <- inherits(response, "Surv") || 
      (is.call(response) && as.character(response[[1]]) == "Surv")
    
    if (is_surv) {
      stop(paste(
        "Survival data (Surv objects) are not yet supported in outstandR v2.0.0.",
        "This feature is scheduled for a future version."
      ), call. = FALSE)
    }
  }
  
  invisible()
}
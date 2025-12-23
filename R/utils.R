# Utility functions for outstandR package

#' Create counterfactual datasets for treatment comparison
#'
#' Creates two copies of a dataset with treatment variable set to different values.
#'
#' @param x_star Base dataset to create counterfactuals from
#' @param trt_var Name of the treatment variable
#' @param comp_trt Comparator treatment value
#' @param ref_trt Reference treatment value
#' @return A list with two elements:
#'   \describe{
#'     \item{comp}{Dataset with treatment set to comp_trt}
#'     \item{ref}{Dataset with treatment set to ref_trt}
#'   }
#' @keywords internal
#'
create_counterfactual_datasets <- function(x_star, trt_var, comp_trt, ref_trt) {
  data_comp <- data_ref <- x_star
  
  # intervene on treatment while keeping set covariates fixed
  data_comp[[trt_var]] <- comp_trt
  data_ref[[trt_var]] <- ref_trt
  
  list(comp = data_comp, ref = data_ref)
}


#' Extract treatment coefficient name from model coefficients
#'
#' Safely extracts the treatment coefficient name from model coefficients,
#' handling potential factor level appends.
#'
#' @param coef_names Names of model coefficients
#' @param trt_var Name of the treatment variable
#' @return The name of the treatment coefficient
#' @keywords internal
#'
extract_treatment_coef_name <- function(coef_names, trt_var) {
  grep(pattern = paste0("^", trt_var, "[^:]*$"),
       coef_names, value = TRUE)
}

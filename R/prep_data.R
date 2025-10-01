# prepare data functions

#' Prepare Individual Patient Data
#'
#' @param form Formula
#' @param data Individual patient data
#' 
prep_ipd <- function(form, data) {
  # select data according to formula
  model.frame(form, data = data)
}

#' Prepare Aggregate Level Data
#' 
#' @param form Formula
#' @param data A dataframe of aggregate level data
#' @param trt_var Treatment variable name
#' 
prep_ald <- function(form, data, trt_var = "trt") {
  
  all_term_labels <- attr(terms(form), "term.labels")
  
  # pull out original of transformed variable
  # remove duplicates (since X and I(X^2) will both appear)
  term.labels <- unique(gsub("I\\(([^)]+)\\^2\\)", "\\1", all_term_labels))
  
  # remove treatment
  # interaction separate by :
  term.labels <- unlist(strsplit(term.labels, ":", fixed = TRUE))
  term.labels <- setdiff(term.labels, trt_var)
  
  dplyr::filter(
    data,
    .data$variable %in% c("y", term.labels) | .data$statistic == "N")
}

# Get study comparator treatment names
#
get_comparator <- function(dat, ref_trt, trt_var = "trt") {
  
  all_trt <- levels(as.factor(dat[[trt_var]]))
  all_trt[all_trt != ref_trt]
}

#' Get reference treatment
#'
#' @param ref_trt Reference treatment
#' @param trt Treatment
#' @param ipd_trial A dataframe of IPD
#' @param ald_trial A dataframe of ALD
#' @return String
get_ref_trt <- function(ref_trt, trt, ipd_trial, ald_trial) {
  
  if (!is.na(ref_trt)) {
    # check exists
    if (!ref_trt %in% names(ipd_trial) || !ref_trt %in% names(ald_trial)) {
      stop("Reference treatment not in IPD and ALD.", call. = FALSE)
    }
    
    return(ref_trt)
  }
  
  ref_trt <- intersect(
    unique(ipd_trial[[trt]]),
    unique(ald_trial[[trt]]))
  
  if (length(ref_trt) == 0) {
    stop("No common treatment in aggregate and individual level data.", call. = FALSE)
  }
  
  if (length(ref_trt) > 1) {
    stop("More than one common treatment in aggregate and individual level data.", call. = FALSE)
  }
  
  ref_trt
}

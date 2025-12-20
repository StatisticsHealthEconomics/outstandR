# prepare data functions ---

#' Prepare Individual Patient Data
#'
#' @param form Formula
#' @param data Individual patient data
#' @returns Model frame
#' @keywords internal
prep_ipd <- function(form, data) {
  # select data according to formula
  model.frame(form, data = data)
}

#' Prepare Aggregate Level Data
#' 
#' @param form Formula
#' @param data A dataframe of aggregate level data
#' @param trt_var Treatment variable name
#' @returns A data frame of filtered ALD by variables in formula.
#' @keywords internal
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

#' Get study comparator treatment names
#'
#' @param dat Data frame.
#' @param ref_trt Reference treatment. String.
#' @param trt_var Treatment variable. String default "trt".
#'
#' @returns Comparator string names
#' @keywords internal
#'
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
#' @return String of reference treatment name.
#' @keywords internal
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

#' Prepare Covariate Distributions
#' 
#' Resolves missing distributions and parameters by looking at the ALD.
#' Allows for partial specification (e.g., user specifies one variable, 
#' function auto-detects the rest).
#'
#' @keywords internal
prepare_covariate_distns <- function(formula, 
                                     ald, 
                                     trt_var, 
                                     marginal_distns, 
                                     marginal_params,
                                     verbose = FALSE) {
  
  # 1. Identify all covariates needed (excluding treatment)
  covariate_names <- get_covariate_names(formula)
  covariate_names <- setdiff(covariate_names, trt_var)
  
  # 2. Initialize empty containers if the user provided nothing (NA)
  if (length(marginal_distns) == 1 && is.na(marginal_distns)) {
    marginal_distns <- setNames(rep(NA_character_, length(covariate_names)), covariate_names)
  }
  
  if (length(marginal_params) == 1 && is.na(marginal_params)) {
    marginal_params <- vector("list", length(covariate_names))
    names(marginal_params) <- covariate_names
  }
  
  # 3. Handle unnamed vectors (assume order matches formula)
  if (is.null(names(marginal_distns)) && length(marginal_distns) > 0) {
    # If lengths don't match, we can't safely assume order, but if they do:
    if (length(marginal_distns) == length(covariate_names)) {
      names(marginal_distns) <- covariate_names
    }
    # Note: If user provides unnamed partial vector (e.g. length 1 for 3 vars),
    # it's hard to guess which var it is. You might want to error here or assume it's the first one.
    # For safety, let's leave it; accessing by name ["var"] will return NA if name missing.
  }
  
  # 4. Iterate through EVERY required covariate
  for (cov in covariate_names) {
    var_ald <- dplyr::filter(ald, .data$variable == cov)
    
    if (nrow(var_ald) == 0) {
      stop(paste("No ALD found for covariate:", cov), call. = FALSE)
    }
    
    # --- FILL MISSING DISTRIBUTIONS ---
    # logic: If the user didn't name this specific covariate, or passed NA
    if (is.na(marginal_distns[cov])) {
      
      if ("prop" %in% var_ald$statistic) {
        marginal_distns[cov] <- "binom"
      } else if (all(c("mean", "sd") %in% var_ald$statistic)) {
        marginal_distns[cov] <- "norm"
      } else {
        stop(paste("Cannot auto-detect distribution for", cov, 
                   "- provide it manually in 'marginal_distns'"), call. = FALSE)
      }
    }
    
    # --- FILL MISSING PARAMETERS ---
    if (is.null(marginal_params[[cov]])) {
      dist <- marginal_distns[cov]
      
      # Extract ALD stats
      m <- var_ald$value[var_ald$statistic == "mean"]
      s <- var_ald$value[var_ald$statistic == "sd"]
      p <- var_ald$value[var_ald$statistic == "prop"]
      
      if (dist == "norm") {
        if (length(m) == 0 || length(s) == 0) stop(paste("Need mean/sd for norm:", cov))
        marginal_params[[cov]] <- list(mean = m, sd = s)
        
      } else if (dist == "binom") {
        if (length(p) == 0) stop(paste("Need prop for binom:", cov))
        marginal_params[[cov]] <- list(size = 1, prob = p)
        
      } else if (dist == "gamma") {
        if (length(m) == 0 || length(s) == 0) stop(paste("Need mean/sd for gamma:", cov))
        marginal_params[[cov]] <- list(shape = m^2/s^2, rate = m/s^2)
        
      } else if (dist == "lnorm") {
        if (length(m) == 0 || length(s) == 0) stop(paste("Need mean/sd for lnorm:", cov))
        var_log <- log(1 + (s^2 / m^2))
        marginal_params[[cov]] <- list(meanlog = log(m) - 0.5 * var_log, sdlog = sqrt(var_log))
        
      } else if (dist == "beta") {
        if (length(m) == 0 || length(s) == 0) stop(paste("Need mean/sd for beta:", cov))
        term <- (m * (1 - m) / s^2) - 1
        marginal_params[[cov]] <- list(shape1 = m * term, shape2 = (1 - m) * term)
        
      } else {
        stop(paste("Auto-parameterization not implemented for '", dist, 
                   "'. Please supply marginal_params for: ", cov, sep=""), call. = FALSE)
      }
    }
    
    # Safety for binom size
    if (marginal_distns[cov] == "binom" && is.null(marginal_params[[cov]]$size)) {
      marginal_params[[cov]]$size <- 1
    }
  }
  
  list(distns = marginal_distns, params = marginal_params)
}

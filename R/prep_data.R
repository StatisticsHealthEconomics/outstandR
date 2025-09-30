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

#' Convert from long to wide format
#'
#' @param df A Dataframe of ALD
#'
#' @importFrom tidyr unite pivot_wider
#' @importFrom stringr str_replace
#' @importFrom dplyr filter arrange select bind_cols
#' 
#' @examples
#' df <- 
#'   data.frame(
#'     variable = c("age", "age", "y", "y", "y", "y", "y", "y", "y", "y"),
#'     statistic = c("mean", "sd", "sum", "bar", "sd", "N", "sum", "bar", "sd", "N"),
#'     trt = c(NA, NA, "B", "B", "B", "B", "C", "C", "C", "C"),
#'     value = c(1,1,1,1,1,1,1,1,1,1))
#' 
reshape_ald_to_wide <- function(df) {
  
  variable_df <- df |> 
    dplyr::filter(.data$variable != "y") |> 
    dplyr::arrange(.data$variable) |> 
    dplyr::select(-.data$study) |> 
    tidyr::unite("colname", .data$stat, .data$variable, sep = ".") |> 
    tidyr::pivot_wider(names_from = .data$colname,
                       values_from = .data$value)
  
  y_df <- df |> 
    dplyr::filter(.data$variable == "y") |> 
    tidyr::unite("colname", .data$variable, .data$trt, .data$stat, sep = ".") |> 
    dplyr::select(.data$colname, .data$value) |> 
    tidyr::pivot_wider(names_from = .data$colname,
                       values_from = .data$value)
  
  # Rename columns: y.X.N to N.X
  # i.e. swap study and stat, removes 'y.'
  colnames(y_df) <- colnames(y_df) |> 
    stringr::str_replace("^y\\.(.*?)\\.(.*?)$", "\\2.\\1")
  
  dplyr::bind_cols(variable_df, y_df)
}

#' Convert from wide to long format
#'
#' @param df A dataframe of ALD
#'
#' @importFrom tidyr pivot_longer separate
#' @importFrom dplyr select mutate bind_rows arrange
#' 
reshape_ald_to_long <- function(df) {
  
  # Separate the 'colname' into 'statistic', 'variable', and 'treatment' columns
  variable_cols <- df |> 
    dplyr::select(-starts_with("y"), -starts_with("N")) |> 
    tidyr::pivot_longer(cols = everything(),
                        names_to = "colname",
                        values_to = "value") |> 
    tidyr::separate(.data$colname,
                    into = c("statistic", "variable"), sep = "\\.") |> 
    dplyr::select(variable, statistic, value) |> 
    dplyr::mutate(trt = NA)  # Set study to NA for these rows
  
  N_cols <- df |> 
    dplyr::select(starts_with("N")) |> 
    tidyr::pivot_longer(cols = everything(),
                        names_to = "colname",
                        values_to = "value") |> 
    tidyr::separate(.data$colname,
                    into = c("statistic", "trt"), sep = "\\.") |> 
    dplyr::select(statistic, trt, value) |> 
    dplyr::mutate(variable = NA)
  
  # process the columns related to 'y' (stat, study, and variable)
  y_cols <- df |> 
    dplyr::select(starts_with("y")) |> 
    tidyr::pivot_longer(cols = everything(),
                        names_to = "colname",
                        values_to = "value") |> 
    tidyr::separate(.data$colname,
                    into = c("variable", "trt", "statistic"),
                    sep = "\\.") |> 
    dplyr::select(variable, statistic, trt, value) 
  
  dplyr::bind_rows(variable_cols, y_cols, N_cols) |> 
    dplyr::arrange(.data$variable, .data$statistic, .data$trt)
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

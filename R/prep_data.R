# prepare data functions

#
prep_ipd <- function(form, data) {
  # select data according to formula
  model.frame(form, data = data)
}

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
    variable %in% c("y", term.labels) | statistic == "N")
}

#' Convert from long to wide format
#'
#' @importFrom tidyr unite pivot_wider
#' @importFrom stringr str_replace
#' @importFrom dplyr filter arrange select bind_cols
#' 
#' @examples
#' df <- 
#'   data.frame(
#'     variable = c("age", "age", "y", "y", "y", "y", "y", "y", "y", "y"),
#'     stat = c("mean", "sd", "sum", "bar", "sd", "N", "sum", "bar", "sd", "N"),
#'     study = c(NA, NA, "B", "B", "B", "B", "C", "C", "C", "C"),
#'     value = c(1,1,1,1,1,1,1,1,1,1))
#' 
reshape_ald_to_wide <- function(df) {
  
  variable_df <- df |> 
    dplyr::filter(variable != "y") |> 
    dplyr::arrange(variable) |> 
    dplyr::select(-study) |> 
    tidyr::unite("colname", stat, variable, sep = ".") |> 
    tidyr::pivot_wider(names_from = colname,
                       values_from = value)
  
  y_df <- df |> 
    dplyr::filter(variable == "y") |> 
    tidyr::unite("colname", variable, study, stat, sep = ".") |> 
    dplyr::select(colname, value) |> 
    tidyr::pivot_wider(names_from = colname,
                       values_from = value)
  
  # Rename columns: y.X.N to N.X
  # i.e. swap study and stat, removes 'y.'
  colnames(y_df) <- colnames(y_df) |> 
    stringr::str_replace("^y\\.(.*?)\\.(.*?)$", "\\2.\\1")
  
  dplyr::bind_cols(variable_df, y_df)
}

##TODO: check this
#' Convert from wide to long format
#'
#' @importFrom tidyr pivot_longer separate
#' @importFrom dplyr select mutate bind_rows arrange
#' 
#' @examples
reshape_ald_to_long <- function(df) {
  
  # Separate the 'colname' into 'statistic', 'variable', and 'treatment' columns
  variable_cols <- df |> 
    dplyr::select(-starts_with("y")) |> 
    tidyr::pivot_longer(cols = everything(),
                        names_to = "colname",
                        values_to = "value") |> 
    tidyr::separate(colname, into = c("statistic", "variable"), sep = "\\.") |> 
    dplyr::select(variable, statistic, value) |> 
    dplyr::mutate(treatment = NA)  # Set study to NA for these rows
  
  # process the columns related to 'y' (stat, study, and variable)
  y_cols <- df |> 
    dplyr::select(starts_with("y")) |> 
    tidyr::pivot_longer(cols = everything(),
                        names_to = "colname",
                        values_to = "value") |> 
    tidyr::separate(colname,
                    into = c("variable", "treatment", "statistic"),
                    sep = "\\.") |> 
    dplyr::select(variable, statistic, treatment, value) 
  
  dplyr::bind_rows(variable_cols, y_cols) |> 
    dplyr::arrange(variable, statistic, treatment)
}


# Get study comparator treatment names
#
get_comparator <- function(dat, ref_trt, trt_var = "trt") {
  
  all_trt <- levels(as.factor(dat[[trt_var]]))
  all_trt[all_trt != ref_trt]
}

#' Get reference treatment
#'
get_ref_trt <- function(ref_trt, trt, ipd_trial, ald_trial) {
  
  if (!is.na(ref_trt)) {
    return()
  }
  
  ref_trt <- intersect(
    unique(ipd_trial[[trt]]),
    unique(ald_trial[[trt]]))
  
  if (length(ref_trt) != 1) {
    stop("More than one common treatment in aggregate and individual level data.", call. = FALSE)
  }  
  
  ref_trt
}

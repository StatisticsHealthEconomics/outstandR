
#' Convert aggregate data from long to wide format
#'
#' @param df A Dataframe of ALD
#' @returns Data frame in wide format
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

#' Convert aggregate data from wide to long format
#'
#' @param df A dataframe of ALD
#' @returns Data frame in long format
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
    dplyr::select(.data$variable, .data$statistic, .data$value) |> 
    dplyr::mutate(trt = NA)  # Set study to NA for these rows
  
  N_cols <- df |> 
    dplyr::select(starts_with("N")) |> 
    tidyr::pivot_longer(cols = everything(),
                        names_to = "colname",
                        values_to = "value") |> 
    tidyr::separate(.data$colname,
                    into = c("statistic", "trt"), sep = "\\.") |> 
    dplyr::select(.data$statistic, .data$trt, .data$value) |> 
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
    dplyr::select(.data$variable, .data$statistic, .data$trt, .data$value) 
  
  dplyr::bind_rows(variable_cols, y_cols, N_cols) |> 
    dplyr::arrange(.data$variable, .data$statistic, .data$trt)
}

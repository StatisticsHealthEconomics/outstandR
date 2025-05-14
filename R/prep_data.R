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
  
  mean_names <- paste0("mean.", term.labels)
  sd_names <- paste0("sd.", term.labels)  ##TODO: for maic do we need these?
  term_names <- c(mean_names, sd_names)
  
  # replace outcome variable name
  response_var <- all.vars(form)[1]
  response_names <- gsub(pattern = "y", replacement = response_var,
                         x = c("y.B.sum", "y.B.bar", "y.B.sd", "N.B",
                               "y.C.sum", "y.C.bar", "y.C.sd", "N.C")) 
  
  keep_names <- c(term_names, response_names)
  data_names <- names(data)
  
  data[data_names %in% keep_names]
}

#' Convert from long to wide format
#'
#' @import plyr
#' @import tidyr
#' @import stringr
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
    
    variable_df <- df %>%
      filter(variable != "y") %>%
      arrange(variable) |> 
      select(-study) |> 
      unite("colname", stat, variable, sep = ".") %>%
      pivot_wider(names_from = colname,
                  values_from = value)
    
    y_df <- df %>%
      filter(variable == "y") %>%
      unite("colname", variable, study, stat, sep = ".") %>%
      select(colname, value) %>%
      pivot_wider(names_from = colname,
                  values_from = value)
    
    # Rename columns: y.X.N to N.X
    # i.e. swap study and stat, removes 'y.'
    colnames(y_df) <- colnames(y_df) %>%
      str_replace("^y\\.(.*?)\\.(.*?)$", "\\2.\\1")
    
    bind_cols(variable_df, y_df)
}

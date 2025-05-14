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

##TODO: check this
#' Convert from wide to long format
#'
#' @import tidyr
#' @import stringr
#' 
#' @examples
reshape_ald_to_long <- function(df) {
  
  # Separate the 'colname' into 'stat', 'variable', and 'study' columns
  variable_cols <- df %>%
    select(-starts_with("y")) %>%
    pivot_longer(cols = everything(),
                 names_to = "colname",
                 values_to = "value") %>%
    separate(colname, into = c("stat", "variable"), sep = "\\.") %>%
    select(variable, stat, value) %>%
    mutate(study = NA)  # Set study to NA for these rows
  
  # process the columns related to 'y' (stat, study, and variable)
  y_cols <- df %>%
    select(starts_with("y")) %>%
    pivot_longer(cols = everything(),
                 names_to = "colname",
                 values_to = "value") %>%
    separate(colname, into = c("variable", "study", "stat"), sep = "\\.") %>%
    select(variable, stat, study, value) 
  
  bind_rows(variable_cols, y_cols) %>%
    arrange(variable, stat, study)
}


# Get study comparator treatment names

#
get_ald_comparator <- function(ald, ref_trt = "C") {
  
  pattern <- paste0("^y\\.(?!", ref_trt, ")")
  
  # filter for names starting with "y." but not with "y.C"
  y_non_ref <- grep(pattern, colnames(ald),
                    value = TRUE, perl = TRUE)
  
  # extract treatment of first match
  sub("^y\\.([A-Z]).*", "\\1", y_non_ref[1])
}

#
get_ipd_comparator <- function(ipd, ref_trt = "C", trt_var = "trt") {
  all_trt <- levels(ipd[[trt_var]])
  all_trt[all_trt != ref_trt]
}


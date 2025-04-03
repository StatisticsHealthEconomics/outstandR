# prepare data functions

#
prep_ipd <- function(form, data) {
  # select data according to formula
  model.frame(form, data = data)
}

#
prep_ald <- function(form, data) {
  
  all_term_labels <- attr(terms(form), "term.labels")
  
  # Remove duplicates (since age and I(age^2) will both appear)
  term.labels <- unique(gsub("I\\(([^)]+)\\^2\\)", "\\1", all_term_labels))
  
  mean_names <- paste0("mean.", term.labels)
  sd_names <- paste0("sd.", term.labels)  ##TODO: for maic do we need these?
  term_names <- c(mean_names, sd_names)
  
  # remove treatment labels
  term_names <- sort(term_names[!grepl(pattern = "trt", term_names)])
  
  # replace outcome variable name
  response_var <- all.vars(form)[1]
  response_names <- gsub(pattern = "y", replacement = response_var,
                         x = c("y.B.sum", "y.B.bar", "N.B", "y.C.sum", "y.C.bar", "N.C")) 
  
  keep_names <- sort(c(term_names, response_names))
  data_names <- sort(names(data))
  
  col_ids <- na.omit(match(keep_names, data_names))
  
  data[col_ids]
}

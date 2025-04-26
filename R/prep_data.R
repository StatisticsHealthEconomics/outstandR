# prepare data functions

#
prep_ipd <- function(form, data) {
  # select data according to formula
  model.frame(form, data = data)
}

#
prep_ald <- function(form, data) {
  
  all_term_labels <- attr(terms(form), "term.labels")
  trt_name <- get_treatment_name(form)
  
  # pull out original of transformed variable
  # remove duplicates (since X and I(X^2) will both appear)
  term.labels <- unique(gsub("I\\(([^)]+)\\^2\\)", "\\1", all_term_labels))
  
  # remove treatment
  # interaction separate by :
  term.labels <- unlist(strsplit(term.labels, ":", fixed = TRUE))
  term.labels <- setdiff(term.labels, trt_name)
  
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

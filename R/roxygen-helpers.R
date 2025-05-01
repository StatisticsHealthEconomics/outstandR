
#
reg_args <- function(include_formula = TRUE,
                     include_family = TRUE) {
  out <- character()
  
  if (include_formula) out <- c(out, "@param formula Linear regression `formula` object. Prognostic factors (PF) are main effects and effect modifiers (EM) are
                                interactions with the treatment variable, i.e., y ~ X1 + trt + trt:X2. For covariates as both PF and EM use `*` syntax.")
  
  if (include_family)  out <- c(out, "@param family A 'family' object specifying the distribution and link function (e.g., 'binomial').
                                See stats::family() for more details.")
  
  out
}

#
study_data_args <- function(include_ipd = TRUE,
                            include_ald = TRUE) {
  out <- character()
  
  if (include_ipd) out <- c(out, "@param ipd Individual-level patient data. Dataframe with one row per patient with outcome, treatment and covariate columns.")
  
  if (include_ald) out <- c(out, "@param ald Aggregate-level data. Single row matrix with summary statistics for each covariate and treatment outcomes.
                             The format is 'mean.*' and 'sd.*' for covariates and 'y.*.sum', 'y.*.bar', 'y.*.sd' for treatments B and C.
                             We assume a common distribution for each treatment arm.")
  
  out
}


#' Data set for individual level patient data
#' 
#' This data set contains simulated patient covariate and outcome values.
#' 
#' @name AC_IPD
#' @docType data
#' @format A data list including the variables needed for the ModStanR basic example. The variables are as follows:
#' \describe{
#' \item{list("X1")}{Numeric covariate}
#' \item{list("X2")}{Numeric covariate}
#' \item{list("X3")}{Numeric covariate}
#' \item{list("X4")}{Numeric covariate}
#' \item{list("trt")}{Integer treatment identifier}
#' \item{list("y")}{Integer binary outcome}
#' }
#' @references Remiro‐Azócar A, Heath A, Baio G (2022)
#' 
#' @keywords datasets
NULL

#' Data set for aggregate level patient data
#' 
#' This data set contains summaries of simulated patient covariate and outcome values.
#' 
#' @name BC_ALD
#' @docType data
#' @format A data list including the variables needed for the ModStanR basic example. The variables are as follows:
#' \describe{
#' \item{list("mean.X1")}{Numeric value}
#' \item{list("mean.X2")}{Numeric value}
#' \item{list("mean.X3")}{Numeric value}
#' \item{list("mean.X4")}{Numeric value}
#' \item{list("sd.X1")}{Numeric standard deviation value}
#' \item{list("sd.X2")}{Numeric standard deviation value}
#' \item{list("sd.X3")}{Numeric standard deviation value}
#' \item{list("sd.X4")}{Numeric standard deviation value}
#' \item{list("y.B.sum")}{Integer number of patients with event on B treatment}
#' \item{list("y.B.bar")}{Numeric proportion of patients on B treatment}
#' \item{list("N.B")}{Integer Total number of patients on B treatment}
#' \item{list("y.C.sum")}{Integer number of patients with event on C treatment}
#' \item{list("y.C.bar")}{Numeric proportion of patients on C treatment}
#' \item{list("N.C")}{Integer Total number of patients on C treatment}
#' }
#' @references Remiro‐Azócar A, Heath A, Baio G (2022)
#' 
#' @keywords datasets
NULL

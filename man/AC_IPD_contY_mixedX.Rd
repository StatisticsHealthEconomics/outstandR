% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{AC_IPD_contY_mixedX}
\alias{AC_IPD_contY_mixedX}
\title{Individual-level patient data for continuous outcome, mixed covariates}
\format{
\subsection{\code{y ~ X1 + X3 + X4 + trt + trt:(X2 + X3 + X4)}}{

\describe{
\item{id}{Numeric unique identifier}
\item{X1}{Numeric prognostic factor continuous covariate}
\item{X2}{Numeric prognostic factor and effect modifier binary covariate}
\item{X3}{Numeric prognostic factor and effect modifier continuous covariate}
\item{X4}{Numeric effect modifier binary covariate}
\item{trt}{Factor treatment identifier. Levels A, C}
\item{y}{Integer binary outcome}
\item{true_eta}{Numeric linear predictor}
}
}
}
\source{
Simulated data
}
\usage{
AC_IPD_contY_mixedX
}
\description{
This data set contains simulated patient covariate and outcome values.
Corresponds to ALD data set.
}
\references{
Remiro‐Azócar A, Heath A, Baio G (2022)
}
\keyword{datasets}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gcomp_ml.R
\name{gcomp_ml.boot}
\alias{gcomp_ml.boot}
\title{G-computation maximum likelihood bootstrap}
\usage{
gcomp_ml.boot(data, indices, R, formula = NULL, family, N = 1000, ald)
}
\arguments{
\item{data}{Trial data}

\item{indices}{Indices sampled from rows of \code{data}}

\item{formula}{Linear regression \code{formula} object}

\item{N}{Synthetic sample size for g-computation}
}
\value{
Mean difference in expected log-odds
}
\description{
Using bootstrap resampling, calculates the log odds ratio.
}
\seealso{
\code{\link[=strategy_gcomp_ml]{strategy_gcomp_ml()}}, \code{\link[=gcomp_ml_log_odds_ratio]{gcomp_ml_log_odds_ratio()}}
}
\keyword{internal}

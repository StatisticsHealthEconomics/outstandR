% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mim.R
\name{calc_mim}
\alias{calc_mim}
\title{Multiple imputation marginalization (MIM)}
\usage{
calc_mim(
  strategy,
  ipd,
  ald,
  M = 1000,
  n.chains = 2,
  warmup = 1000,
  iters = 4000
)
}
\arguments{
\item{strategy}{Strategy object}

\item{ipd}{Individual-level data}

\item{ald}{Aggregate-level data}

\item{M}{Number of syntheses used in analysis stage (high for low Monte Carlo error)}

\item{warmup}{Number of warm-up iterations}

\item{iters}{Number of total iterations}

\item{n.chain}{Number of chains}
}
\description{
Multiple imputation marginalization (MIM)
}
\keyword{internal}

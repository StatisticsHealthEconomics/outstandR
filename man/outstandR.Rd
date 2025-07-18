% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outstandR.R
\name{outstandR}
\alias{outstandR}
\title{Calculate the difference between treatments using all evidence}
\usage{
outstandR(
  ipd_trial,
  ald_trial,
  strategy,
  ref_trt = NA,
  CI = 0.95,
  scale = NULL,
  ...
)
}
\arguments{
\item{ipd_trial}{Individual-level patient data. For example, suppose between studies \emph{A} and \emph{C}.
In a long format and must contain a treatment column and outcome column consistent with the formula object.
The labels in the treatment are used internally so there must be a common treatment with the aggregate-level data trial.}

\item{ald_trial}{Aggregate-level data. For example, suppose between studies \emph{B} and \emph{C}. The column names are
\itemize{
\item \code{variable}: Covariate name. In the case of treatment arm sample size this is \code{NA}
\item \code{statistic}: Summary statistic name from "mean", standard deviation "sd", probability "prop", or "sum"
\item \code{value}: Numerical value of summary statistic
\item \code{trt}: Treatment label. Because we assume a common covariate distribution between treatment arms this is \code{NA}
}}

\item{strategy}{Computation strategy function. These can be
\code{strategy_maic()}, \code{strategy_stc()}, \code{strategy_gcomp_ml()} and \code{strategy_gcomp_stan()}.}

\item{ref_trt}{Reference / common / anchoring treatment name.}

\item{CI}{Confidence interval; between 0,1}

\item{scale}{Relative treatment effect scale. If \code{NULL}, the scale is automatically determined from the model.
Choose from "log-odds", "log_relative_risk", "risk_difference", "delta_z", "mean_difference", "rate_difference" depending on the data type.}

\item{...}{Additional arguments. Currently, can pass named arguments to \code{rstanarm::stan_glm()} via \code{strategy_gcomp_stan()}.}
}
\value{
List of length 3 of statistics as a \code{outstandR} class object.
Containing statistics between each pair of treatments.
These are the mean, variances and confidence intervals,
for contrasts and absolute values.
}
\description{
This is the main, top-level wrapper for \code{{outstandR}}.
Methods taken from
\insertCite{RemiroAzocar2022}{outstandR}.
}
\examples{
data(AC_IPD)  # AC patient-level data
data(BC_ALD)  # BC aggregate-level data

# linear formula
lin_form <- as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")

# matching-adjusted indirect comparison
outstandR_maic <- outstandR(AC_IPD, BC_ALD,
                            strategy = strategy_maic(formula = lin_form))

# simulated treatment comparison
outstandR_stc <- outstandR(AC_IPD, BC_ALD,
                           strategy = strategy_stc(lin_form))

# G-computation with maximum likelihood
# outstandR_gcomp_ml <- outstandR(AC_IPD, BC_ALD,
                                  strategy = strategy_gcomp_ml(lin_form))

# G-computation with Bayesian inference
outstandR_gcomp_stan <- outstandR(AC_IPD, BC_ALD,
                                  strategy = strategy_gcomp_stan(lin_form))

# Multiple imputation marginalization
outstandR_mim <- outstandR(AC_IPD, BC_ALD,
                           strategy = strategy_mim(lin_form))

}
\references{
\insertRef{RemiroAzocar2022}{outstandR}
}
\seealso{
\code{\link[=strategy_maic]{strategy_maic()}} \code{\link[=strategy_stc]{strategy_stc()}} \code{\link[=strategy_gcomp_ml]{strategy_gcomp_ml()}} \code{\link[=strategy_gcomp_stan]{strategy_gcomp_stan()}}
}

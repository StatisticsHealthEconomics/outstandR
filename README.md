# outstandR: *Out*come regression *stand*ardisation
<img src="man/figures/logo.png" align="right" />
================

<!-- <img align="right" src="mime.png" width="100"> -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/StatisticsHealthEconomics/outstandR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/StatisticsHealthEconomics/outstandR/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![outstandR status
badge](https://statisticshealtheconomics.r-universe.dev/outstandR/badges/version)](https://statisticshealtheconomics.r-universe.dev/outstandR)

<!-- badges: end -->

> Indirect treatment comparison with limited subject-level data

## Overview

`{outstandR}` is an R package designed to facilitate **outcome regression
standardisation** using model-based approaches, particularly focusing on
G-estimation. The package provides tools to apply standardisation
techniques for indirect treatment comparisons, especially in scenarios
with limited individual patient data.

## Who is this package for?

The target audience of `{outstandR}` is those who want to perform
model-based standardization in the specific context of two-study
indirect treatment comparisons with limited subject-level data. This is
model-based standardization with two additional steps:

1.  Covariate simulation (to overcome limited subject-level data for one
    of the studies)
2.  Indirect comparison across studies

## Installation

Install the [development version from
GitHub](https://github.com/StatisticsHealthEconomics/) using R-universe:

``` r
install.packages("outstandR", repos = c("https://statisticshealtheconomics.r-universe.dev", "https://cloud.r-project.org"))
```

Alternatively, you may wish to download directly from the repo with
`remotes::install_github("StatisticsHealthEconomics/outstandR")`.

## Background

Population adjustment methods are increasingly used to compare marginal
treatment effects when there are cross-trial differences in effect
modifiers and limited patient-level data.

The `{outstandR}` package allows the implementation of a range of
methods for this situation including the following:

-   *Matching-Adjusted Indirect Comparison (MAIC)* is based on
    propensity score weighting, which is sensitive to poor covariate
    overlap and cannot extrapolate beyond the observed covariate space.
    It reweights the individual patient-level data (IPD) to match the
    aggregate characteristics of the comparator trial, thereby aligning
    the populations.

-   *Simulated Treatment Comparison (STC)* relies on outcome regression
    models fitted to IPD, conditioning on covariates to estimate the
    effect of treatment. These estimates are then applied to the
    aggregate-level comparator population. Like MAIC, STC is limited by
    its conditional nature and can produce biased marginal estimates if
    not properly marginalized.

-   *Parametric G-computation with maximum likelihood*: This method fits
    an outcome model to the IPD using maximum likelihood estimation,
    then uses that model to predict outcomes in the comparator
    population. It allows extrapolation beyond the observed covariate
    space but requires correct specification of the outcome model to
    avoid bias.

-   *Parametric G-computation with Bayesian inference*: Similar to the
    maximum likelihood version, this approach fits an outcome model but
    within a Bayesian framework. It allows coherent propagation of
    uncertainty through prior distributions and posterior inference,
    enabling probabilistic sensitivity analysis and full uncertainty
    quantification.

-   *Multiple imputation marginalization method based on parametric G-computation*: Current
    outcome regression-based alternatives can extrapolate but target a
    conditional treatment effect that is incompatible in the indirect
    comparison. When adjusting for covariates, one must integrate or
    average the conditional estimate over the relevant population to
    recover a compatible marginal treatment effect. This can be easily
    applied where the outcome regression is a generalized linear model
    or a Cox model. The approach views the covariate adjustment
    regression as a nuisance model and separates its estimation from the
    evaluation of the marginal treatment effect of interest. The method
    can accommodate a Bayesian statistical framework, which naturally
    integrates the analysis into a probabilistic framework.

## General problem

Consider one trial, for which the company has IPD, comparing treatments
*A* and *C*, from herein call the *AC* trial. Also, consider a second
trial comparing treatments *B* and *C*, similarly called the *BC* trial.
For this trial only published aggregate data are available. We wish to
estimate a comparison of the effects of treatments *A* and *B* on an
appropriate scale in some target population *P*, denoted by the
parameter $d_{AB(P)}$. We make use of bracketed subscripts to denote a
specific population. Within the *BC* population there are parameters
$\mu_{B(BC)}$ and $\mu_{C(BC)}$ representing the expected outcome on
each treatment (including parameters for treatments not studied in the
*BC* trial, e.g. treatment *A*). The *BC* trial provides estimators
$\bar Y_{B(BC)}$ and $\bar Y_{C(BC)}$ of $\mu_{B(BC)}$, $\mu_{C(BC)}$,
respectively, which are the summary outcomes. It is the same situation
for the *AC* trial.

For a suitable scale, for example a log-odds ratio, or risk difference,
we form estimators $\Delta_{BC(BC)}$ and $\Delta_{AC(AC)}$ of the trial
level (or marginal) relative treatment effects. We shall assume that
this is always represented as a difference so, for example, for the risk
ratio this is on the log scale.

$$
\Delta_{AC{(AC)}} = g(\bar{Y}_{C{(AC)}}) - g(\bar{Y}_{A{(AC)}})
$$

and similarly for $\Delta_{BC(BC)}$. If we assume that there is no difference in effect modifiers between trials, then the estimator of the relative treatment effect $d_{AB(BC)}$ is

$$
\Delta_{AB(BC)} = \Delta_{BC(BC)} - \Delta_{AC(AC)}.
$$

However, when distributions of the effect modifiers are different between trial populations, the relative treatment effect estimated from each trial cannot simply be combined as above. The purpose of population-adjustment in ITC is to include an estimate for $\Delta_{AC(BC)}$.

## References

This R package contains code originally written for the papers:

> Remiro-Azócar, A., Heath, A. & Baio, G. (2022) *Parametric
> G-computation for Compatible Indirect Treatment Comparisons with
> Limited Individual Patient Data.* Res Synth Methods;1–31.

and

> Remiro-Azócar, A., Heath, A., & Baio, G. (2023) *Model-based
> standardization using multiple imputation. BMC Medical Research
> Methodology*, 1–15. <https://doi.org/10.1186/s12874-024-02157-x>

## Contributing

We welcome contributions! Please submit contributions through
`Pull Requests`, following the [contributing
guidelines](https://github.com/n8thangreen/BCEA/blob/dev/CONTRIBUTING.md).
To report issues and/or seek support, please file a new ticket in the
[issue](https://github.com/StatisticsHealthEconomics/outstandR/issues)
tracker.

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/n8thangreen/BCEA/blob/dev/CONDUCT.md). By
participating in this project you agree to abide by its terms.

> [!NOTE]
> This package is licensed under the GPLv3. For more
> information, see [LICENSE](https://www.gnu.org/licenses/gpl-3.0).

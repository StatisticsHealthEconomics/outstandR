<div id="main" class="col-md-9" role="main">

# outstandR: *Out*come regression *stand*ardisation

<div class="section level1">

  

> Indirect treatment comparison with limited subject-level data

<div class="section level2">

## Overview

[outstandR](https://StatisticsHealthEconomics.github.io/outstandR/) is
an R package designed to facilitate a unified framework for
**Population-Adjusted Indirect Comparison (PAIC)** in Health Technology
Assessment (HTA). It facilitates outcome regression standardisation
using model-based approaches when head-to-head clinical trials are
absent and trial populations differ in effect modifiers.

By streamlining the workflow of covariate simulation, model
standardisation, and contrast estimation,
[outstandR](https://StatisticsHealthEconomics.github.io/outstandR/)
enables robust and compatible evidence synthesis, particularly in
scenarios where individual patient data (IPD) is limited to only one of
the trials being compared.

</div>

<div class="section level2">

## Who is this package for?

The target audience of
[outstandR](https://StatisticsHealthEconomics.github.io/outstandR/)
includes statisticians and health economists performing indirect
treatment comparisons requiring cross-trial population adjustment. It
simplifies the implementation of model-based standardization with two
core steps:

1.  Covariate simulation (to overcome limited subject-level data for
    aggregate-level data studies).
2.  Indirect comparison across studies targeting compatible marginal
    treatment effects.

</div>

<div class="section level2">

## Installation

Install the released version from CRAN:

<div id="cb1" class="sourceCode">

``` r
install.packages("outstandR")
```

</div>

Or install the [development version from
GitHub](https://github.com/StatisticsHealthEconomics/) using R-universe:

<div id="cb2" class="sourceCode">

``` r
install.packages("outstandR", 
  repos = c("https://statisticshealtheconomics.r-universe.dev", 
            "https://cloud.r-project.org"))
```

</div>

Alternatively, you may wish to download directly from the repo with

<div id="cb3" class="sourceCode">

``` r
remotes::install_github("StatisticsHealthEconomics/outstandR")
```

</div>

</div>

<div class="section level2">

## Background

Population adjustment methods are increasingly used to compare marginal
treatment effects when there are cross-trial differences in effect
modifiers and limited patient-level data.

The [outstandR](https://StatisticsHealthEconomics.github.io/outstandR/)
package allows the implementation of a range of methods for this
situation including the following:

-   *Matching-Adjusted Indirect Comparison (MAIC)* is based on
    weighting, which is sensitive to poor covariate overlap and cannot
    extrapolate beyond the observed covariate space. It reweights the
    individual patient-level data (IPD) to match the aggregate
    characteristics of the comparator trial, thereby aligning the
    populations.

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

-   *Multiple imputation marginalization method based on parametric
    G-computation*: Current outcome regression-based alternatives can
    extrapolate but target a conditional treatment effect that is
    incompatible in the indirect comparison. When adjusting for
    covariates, one must integrate or average the conditional estimate
    over the relevant population to recover a compatible marginal
    treatment effect. This can be easily applied where the outcome
    regression is a generalized linear model or a Cox model. The
    approach views the covariate adjustment regression as a nuisance
    model and separates its estimation from the evaluation of the
    marginal treatment effect of interest. The method can accommodate a
    Bayesian statistical framework, which naturally integrates the
    analysis into a probabilistic framework.

-   *Simulated Treatment Comparison (STC) \[DEPRECATED\]:* While
    previously common, standard STC can produce biased marginal
    estimates due to aggregation bias and non-collapsibility when using
    non-linear link functions. *Note: STC is formally deprecated as of
    v2.0.0 in favour of the more robust G-computation approaches.*

</div>

<div class="section level2">

## General problem

Consider one trial, for which the company has IPD, comparing treatments
*A* and *C*, from herein call the *AC* trial. Also, consider a second
trial comparing treatments *B* and *C*, similarly called the *BC* trial.
For this trial only published aggregate data are available. We wish to
estimate a comparison of the effects of treatments *A* and *B* on an
appropriate scale in some target population *P*, denoted by the
parameter *d*_(*A**B*(*P*)). We make use of bracketed subscripts to
denote a specific population. Within the *BC* population there are
parameters *μ*_(*B*(*B**C*)) and *μ*_(*C*(*B**C*)) representing the
expected outcome on each treatment (including parameters for treatments
not studied in the *BC* trial, e.g. treatment *A*). The *BC* trial
provides estimators *Ȳ*_(*B*(*B**C*)) and *Ȳ*_(*C*(*B**C*)) of
*μ*_(*B*(*B**C*)), *μ*_(*C*(*B**C*)), respectively, which are the
summary outcomes. It is the same situation for the *AC* trial.

For a suitable scale, for example a log-odds ratio, or risk difference,
we form estimators *Δ̂**B**C*(*B**C*) and *Δ̂*_(*A**C*(*A**C*)) of the
trial level (or marginal) relative treatment effects. We shall assume
that this is always represented as a difference so, for example, for the
risk ratio this is on the log scale.

*Δ̂*_(*A**C*(*A**C*)) = *g*(*Ȳ*_(*C*(*A**C*))) − *g*(*Ȳ*_(*A*(*A**C*)))

and similarly for *Δ̂*_(*B**C*(*B**C*)). If we assume that there is no
difference in effect modifiers between trials, then the estimator of the
relative treatment effect *d*_(*A**B*(*B**C*)) is

*Δ̂*_(*A**B*(*B**C*)) = *Δ̂*_(*B**C*(*B**C*)) − *Δ̂*_(*A**C*(*B**C*)).

However, when distributions of the effect modifiers are different
between trial populations, the relative treatment effect estimated from
each trial cannot simply be combined as above. The purpose of
population-adjustment in ITC is to include an estimate for
*Δ̂*_(*A**C*(*B**C*)).

</div>

<div class="section level2">

## References

This R package contains code originally written for the papers:

> Remiro-Azócar, A., Heath, A. & Baio, G. (2022) *Parametric
> G-computation for Compatible Indirect Treatment Comparisons with
> Limited Individual Patient Data.* Res Synth Methods;1–31.

and

> Remiro-Azócar, A., Heath, A., & Baio, G. (2023) *Model-based
> standardization using multiple imputation. BMC Medical Research
> Methodology*, 1–15. <https://doi.org/10.1186/s12874-024-02157-x>

</div>

<div class="section level2">

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

> \[!NOTE\] This package is licensed under the GPLv3. For more
> information, see [LICENSE](https://www.gnu.org/licenses/gpl-3.0).

</div>

</div>

</div>

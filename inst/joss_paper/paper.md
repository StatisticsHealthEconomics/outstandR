---
title: 'outstandR: An R Package for Indirect Treatment Comparison with Limited Subject-Level Data'
tags:
  - R
  - health economics
  - outcome standardisation
  - Bayesian inference
authors:
  - name: Nathan Green
    orcid: 0000-0003-2745-1736
    equal-contrib: true
    affiliation: 1
  - name: Antonio Remiro-Acócar
    orcid: 
    equal-contrib: true
    affiliation: 2
  - name: Gianluca Baio
    corresponding: true
    affiliation: 1
affiliations:
  - name: University College London (UCL), UK
    index: 1
    ror: 02jx3x895
  - name: Novonodisk
    index: 2
    ror: 
date: 1 August 2025
bibliography: paper.bib
editor_options: 
  markdown: 
    wrap: 72
---

# Summary

The goal of `outstandR` is to provide functionality to perform
population adjustment methods are increasingly used to compare marginal
treatment effects when there are cross-trial differences in effect
modifiers and limited patient-level data [@RemiroAzocar2022a]. This
presents a significant challenge in evidence synthesis, particularly in
health technology assessment. Existing methods often face limitations
such as sensitivity to poor covariate overlap or inability to
extrapolate beyond observed covariate spaces. The `outstandR` package
addresses these challenges by offering a robust framework for out- come
regression standardisation, focusing on G-estimation. It enables
researchers to perform model-based standardization with two additional
crucial steps: covariate simulation (to over- come limited subject-level
data for one of the studies) and indirect comparison across studies
[@RemiroAzocar2022a]. The target audience of `outstandR` is those who
want to perform model-based standardization in the specific context of
two-study indirect treatment comparisons with limited subject-level
data.

# Statement of need

Indirect treatment comparisons (ITCs) are a cornerstone of modern
evidence synthesis, especially in health technology assessment (HTA)
where decision-makers must compare novel treatments against a range of
competitors. A common and challenging scenario arises when individual
patient data (IPD) is available for one trial, but only aggregate-level
data (ALD) is available for the comparator trial. Naively comparing
these studies can introduce significant bias because of differences in
patient populations. There is a clear need for a unified and robust
software tool that implements a range of modern population adjustment
techniques to address these challenges systematically.

# Method

We developed the `outstandR` R package to provide a comprehensive
framework for performing anchored ITCs using a suite of population
adjustment methods, with a focus on robust G-computation techniques.
`outstandR` streamlines the entire analytical workflow: from fitting
outcome models on IPD and standardizing them to the ALD population, to
performing the final indirect comparison. By implementing multiple
methods — including Matching-Adjusted Indirect Comparison (MAIC),
Simulated Treatment Comparison (STC) [@phillippo2016methods], parametric
G-computation (both frequentist and Bayesian) [@RemiroAzocar2022a], and
the Multiple Imputation Marginalization (MIM) method
[@RemiroAzocar2022b] — within a single interface, our package empowers
researchers to conduct sensitivity analyses and select the most
appropriate approach for their data. `outstandR` lowers the technical
barrier to entry for these complex analyses, promoting more reliable and
transparent evidence synthesis for healthcare decision-making.

Related R packages we are aware of are more general-purpose tools
implementing specific methods. The `marginaleffects` [@aari_marginaleffects_2024] package is not
designed for population adjustment between studies. `stdReg2` focuses on
standardising outcomes with a single data set [@g_sofer_2023_10022204].
For G-formula
implementations, `gfoRmula` can estimate effects in the presence of
time-varying treatments and confounders [@sjolander2023gformula]. It is designed for estimating
causal effects from longitudinal data with one study. `gFormulaMI`
employs multiple imputation using the `mice` package [@Sterne2023]. Finally,
`maicplus` is a specialist ITC package but focused only on the MAIC
approach [@maicplus_2024].

Our analysis performs an anchored indirect treatment comparison (ITC) to
estimate the relative effect of two treatments, A and B. This comparison
uses individual patient data (IPD) from a trial comparing treatments A
and C (the AC trial) and aggregate level data (ALD) from a trial
comparing treatments B and C (the BC trial). To account for potential
differences in the patient populations between the two trials, which can
bias the ITC, we use population adjustment methods. These methods
standardize the results from the AC trial to the baseline
characteristics of the BC trial population.

The general procedure involves two main stages. First, we fit an outcome
regression model using the IPD from the AC trial. This model describes
the relationship between the outcome, treatment assignment, and a set of
baseline covariates, including both prognostic factors and treatment
effect modifiers. Second, we use this fitted model to predict the
outcomes for treatment A and C in a target population that reflects the
aggregate baseline characteristics of the BC trial. For G-computation,
this step can involve simulating a large synthetic dataset that mirrors
the covariate distributions of the BC population. The resulting adjusted
treatment effect, $\Delta$AC(BC), is then indirectly compared against
the observed effect from the ALD, $\Delta$BC(BC), to yield the final
adjusted estimate for A versus B, $\Delta$AB(BC). Uncertainty is
quantified using non-parametric bootstrapping for frequentist methods or
by propagating parameter uncertainty from the full posterior
distribution in Bayesian implementations.

# References

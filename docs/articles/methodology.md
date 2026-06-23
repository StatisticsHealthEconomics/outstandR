<div id="main" class="col-md-9" role="main">

# Methodological Background

<div class="section level2">

## Introduction

Indirect treatment comparisons (ITCs) are vital when head-to-head
clinical trials are absent. Unadjusted comparisons are biased when trial
populations differ in effect modifiers; thus, population adjustment is
required. `outstandR` provides a unified framework for these methods.

</div>

<div class="section level2">

## Matching-Adjusted Indirect Comparison (MAIC)

MAIC is a method-of-moments weighting approach designed to match the
aggregate characteristics of the comparator trial. Weights *w*_(*i*) for
each patient *i* in the IPD trial are derived via a logistic regression
formulation:

*w*_(*i*) = exp (*β*^(⊤)*X*_(*i*))

The parameters *β* are found by minimizing a convex objective function
to match target covariate means:

$$\\min\\limits\_{\\beta \\in \\mathbb{R}^{p}}\\sum\\limits\_{i = 1}^{n\_{A}}\\exp\\left( \\beta^{\\top}\\left( X\_{i} - {\\overline{X}}\_{B} \\right) \\right)$$

</div>

<div class="section level2">

## Outcome Regression Models

For an individual *i* with outcome *Y*_(*i*), treatment
*T*_(*i*) ∈ {0, 1}, and baseline covariates *X*_(*i*), the general
outcome model is:

*g*(𝔼\[*Y*_(*i*)\|*X*_(*i*),*T*_(*i*)\]) = *α* + *X*_(*i*)^(⊤)*β*₀ + (*β*_(*t**r**t*)+*X*_(*i*)^(⊤)*β*₁)𝕀(*T*_(*i*)=1)

</div>

<div class="section level2">

## Parametric G-computation (Maximum Likelihood)

G-computation standardizes outcomes across treatment regimens. We
estimate the marginal mean outcome by simulating a pseudo-population of
size *N* reflecting the target ALD trial, predicting outcomes, and
averaging:

$$\\widehat{\\mathbb{E}}\\left\\lbrack Y^{T}\|\\widehat{\\theta} \\right\\rbrack = \\frac{1}{N}\\sum\\limits\_{i = 1}^{N}\\widehat{Y\_{i}}\\left( T,\\widehat{\\theta} \\right)$$

</div>

<div class="section level2">

## Parametric G-computation (Bayesian Inference)

Bayesian G-computation estimates the full posterior distribution of the
model parameters, offering robust uncertainty quantification. For *M*
posterior samples *θ*^((*m*)), the marginal mean is computed as:

$$\\widehat{\\mathbb{E}}\\left\\lbrack Y^{T}\|\\theta^{(m)} \\right\\rbrack = \\frac{1}{N}\\sum\\limits\_{i = 1}^{N}{\\widehat{Y}}\_{i}\\left( T,\\theta^{(m)} \\right)$$

</div>

<div class="section level2">

## Multiple Imputation Marginalization (MIM)

MIM conceptualizes unobserved potential outcomes as a missing data
problem. Using *M* posterior draws, MIM imputes counterfactuals and
utilizes modified Rubin’s rules to calculate total variance by
subtracting within-imputation variance ($\\overline{U}$) from
between-imputation variance (*B*):

$$Var\\left( \\overline{Q} \\right) = B - \\overline{U}$$

</div>

</div>

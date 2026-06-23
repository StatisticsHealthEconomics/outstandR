<div id="main" class="col-md-9" role="main">

# Simulated treatment comparison (STC)

<div class="ref-description section level2">

**\[deprecated\]**

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
strategy_stc(
  formula = NULL,
  family = gaussian(link = "identity"),
  trt_var = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   formula:

    Model formula.

-   family:

    Model family.

-   trt_var:

    Treatment variable name.

</div>

<div class="section level2">

## Value

`stc` class object

</div>

<div class="section level2">

## Details

`strategy_stc()` was deprecated in outstandR version 1.X.X. We recommend
using G-computation (`strategy_gcomp_ml()`) as a more robust alternative
for this type of analysis.

Outcome regression-based method which targets a conditional treatment
effect. STC is a modification of the covariate adjustment method. An
outcome model is fitted using IPD in the *AB* trial. For example,

$$ g(\\mu\_{t(AB)}(X)) = \\beta_0 + \\beta_1^T X + (\\beta_B +
\\beta_2^T X^{EM}) I(t=B) $$ where \\(\\beta_0\\) is an intercept term,
\\(\\beta_1\\) is a vector of coefficients for prognostic variables,
\\(\\beta_B\\) is the relative effect of treatment *B* compared to *A*
at \\(X=0\\), \\(\\beta_2\\) is a vector of coefficients for effect
modifiers \\(X^{EM}\\) subvector of the full covariate vector \\(X\\)),
and \\(\\mu\_{t(AB)}(X)\\) is the expected outcome of an individual
assigned treatment \\(t\\) with covariate values \\(X\\) which is
transformed onto a chosen linear predictor scale with link function
\\(g(\\cdot)\\).

</div>

</div>

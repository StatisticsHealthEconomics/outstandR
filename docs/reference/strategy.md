<div id="main" class="col-md-9" role="main">

# New strategy objects

<div class="ref-description section level2">

Create a type of strategy class for each modelling approach.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
strategy_maic(
  formula = NULL,
  family = gaussian(link = "identity"),
  trt_var = NULL,
  n_boot = 1000L,
  moments = 1,
  int = FALSE,
  verbatim = TRUE
)

strategy_gcomp_ml(
  formula = NULL,
  family = gaussian(link = "identity"),
  trt_var = NULL,
  rho = NA,
  marginal_distns = NA,
  marginal_params = NA,
  n_boot = 1000L,
  N = 1000L
)

strategy_gcomp_bayes(
  formula = NULL,
  family = gaussian(link = "identity"),
  trt_var = NULL,
  rho = NA,
  marginal_distns = NA,
  marginal_params = NA,
  N = 1000L
)

strategy_mim(
  formula = NULL,
  family = gaussian(link = "identity"),
  trt_var = NULL,
  rho = NA,
  marginal_distns = NA,
  marginal_params = NA,
  N = 1000L
)

new_strategy(strategy, ...)
```

</div>

</div>

<div class="section level2">

## Arguments

-   formula:

    Linear regression `formula` object. Prognostic factors (PF) are main
    effects and effect modifiers (EM) are interactions with the
    treatment variable, e.g., y \~ X1 + trt + trt:X2. For covariates as
    both PF and EM use `*` syntax.

-   family:

    A 'family' object specifying the distribution and link function
    (e.g., 'binomial'). See stats::family() for more details.

-   trt_var:

    Treatment variable name; string

-   n_boot:

    The number of resamples used for the non-parametric bootstrap;
    integer

-   moments:

    Integer. The number of moments of the covariates to balance. Setting
    `moments = 2` includes both the original variables and their squared
    terms, which effectively balances both the means and the variances.
    Default to 1.

-   int:

    Logical. If `TRUE`, includes two-way interactions between all
    covariates in the balancing model to effectively balance their
    covariances. Default `FALSE`

-   verbatim:

    Logical. Output to console during running.

-   rho:

    A named square matrix of covariate correlations; default NA

-   marginal_distns:

    Marginal distributions names; vector default NA. Available
    distributions are given in stats::Distributions. See
    `copula::Mvdc()` for details

-   marginal_params:

    Marginal distributions parameters; list of lists, default NA. See
    `copula::Mvdc()` for details

-   N:

    Synthetic sample size for g-computation

-   strategy:

    Class name from `strategy_maic`, `strategy_stc`,
    `strategy_gcomp_ml`, `strategy_gcomp_bayes`, `strategy_mim`

-   ...:

    Additional arguments

</div>

<div class="section level2">

## Value

`maic` class object

`gcomp_ml` class object

`gcomp_bayes` class object

`mim` class object

Strategy list object

</div>

<div class="section level2">

## Note

While current implementations focus on binary, continuous, and count
outcomes, support for survival data (using the `survival` package) is
under active development and scheduled for a future version.

</div>

<div class="section level2">

## Matching-adjusted indirect comparison (MAIC)

MAIC is a form of non-parametric likelihood reweighting method which
allows the propensity score logistic regression model to be estimated
without IPD in the *AC* population. The mean outcomes
\\(\\mu\_{t(AC)}\\) on treatment \\(t = A,B\\) in the *AC* target
population are estimated by taking a weighted average of the outcomes
\\(Y\\) of the \\(N\\) individuals in arm \\(t\\) of the *AB*
population.

Used to compare marginal treatment effects where there are cross-trial
differences in effect modifiers and limited patient-level data.

$$ \\hat{Y}\_{} = \\frac{\\sum\_{i=1}^{N} Y\_{it(AB)}
w\_{it}}{\\sum\_{i=1}^{N} w\_{it}} $$ where the weight \\(w\_{it}\\)
assigned to the \\(i\\)-th individual receiving treatment \\(t\\) is
equal to the odds of being enrolled in the *AC* trial vs the *AB* trial.

</div>

<div class="section level2">

## G-computation maximum likelihood

G-computation marginalizes the conditional estimates by separating the
regression modelling from the estimation of the marginal treatment
effect for *A* versus *C*. For example, a regression model of the
observed outcome \\(y\\) on the covariates \\(x\\) and treatment \\(z\\)
is fitted to the *AC* IPD:

$$ g(\\mu_n) = \\beta_0 + \\boldsymbol{x}\_n \\boldsymbol{\\beta_1} +
(\\beta_z + \\boldsymbol{x_n^{EM}} \\boldsymbol{\\beta_2})
\\mbox{I}(z_n=1) $$ In the context of G-computation, this regression
model is called the “Q-model". Having fitted the Q-model, the regression
coefficients are treated as nuisance parameters. The parameters are
applied to the simulated covariates \\(x\*\\) to predict hypothetical
outcomes for each subject under both possible treatments. Namely, a pair
of predicted outcomes, also called potential outcomes, under *A* and
under *C*, is generated for each subject.

By plugging treatment *C* into the regression fit for every simulated
observation, we predict the marginal outcome mean in the hypothetical
scenario in which all units are under treatment *C*:

$$ \\hat{\\mu}\_0 = \\int\_{x^\*} g^{-1} (\\hat{\\beta}\_0 + x^\*
\\hat{\\beta}\_1 ) p(x^\*) dx^\* $$ To estimate the marginal or
population-average treatment effect for *A* versus *C* in the linear
predictor scale, one back-transforms to this scale the average
predictions, taken over all subjects on the natural outcome scale, and
calculates the difference between the average linear predictions:

$$ \\hat{\\Delta}^{(2)}\_{10} = g(\\hat{\\mu}\_1) - g(\\hat{\\mu}\_0) $$

</div>

<div class="section level2">

## G-computation Bayesian

The difference between Bayesian G-computation and its maximum-likelihood
counterpart is in the estimated distribution of the predicted outcomes.
The Bayesian approach also marginalizes, integrates or standardizes over
the joint posterior distribution of the conditional nuisance parameters
of the outcome regression, as well as the joint covariate distribution.

Draw a vector of size \\(N\*\\) of predicted outcomes \\(y\*z\\) under
each set intervention \\(z\* \\in \\{0, 1\\}\\) from its posterior
predictive distribution under the specific treatment. This is defined as
\\(p(y\*\_{z\*} \| \\mathcal{D}\_{AC}) = \\int\_{\\beta} p(y\*\_{z\*} \|
\\beta) p(\\beta \| \\mathcal{D}\_{AC}) d\\beta\\) where \\(p(\\beta \|
\\mathcal{D}\_{AC})\\) is the posterior distribution of the outcome
regression coefficients \\(\\beta\\), which encode the predictor-outcome
relationships observed in the *AC* trial IPD.

This is given by:

$$ p(y\*\_{z\*} \\mid \\mathcal{D}\_{AC}) = \\int\_{x\*} p(y\* \\mid
z\*, x\*, \\mathcal{D}\_{AC}) p(x\* \\mid \\mathcal{D}\_{AC}) dx\* $$

$$ = \\int\_{x\*} \\int\_{\\beta} p(y\* \\mid z\*, x\*, \\beta) p(x\*
\\mid \\beta) p(\\beta \\mid \\mathcal{D}\_{AC}) d\\beta dx\* $$ In
practice, the integrals above can be approximated numerically, using
full Bayesian estimation via Markov chain Monte Carlo (MCMC) sampling.

</div>

<div class="section level2">

## Multiple imputation marginalization (MIM)

MIM targets a marginal treatment effect by using parametric
G-computation within a multiple imputation framework. This approach
views the covariate adjustment regression as a nuisance model and
separates its estimation from the evaluation of the marginal treatment
effect of interest. It is particularly useful for ensuring compatibility
in indirect comparisons when adjusting for effect modifiers.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`strategy_gcomp_bayes()`

`strategy_gcomp_ml()` `copula::Mvdc()`

</div>

</div>

</div>

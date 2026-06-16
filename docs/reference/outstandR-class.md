# outstandR class

The `outstandR` class contains the results from running a model with the
function
[`outstandR()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/outstandR.md).

## Details

Objects of class `outstandR` have the following

- contrasts:

  A list containing statistics for relative treatment effects:

  - `means`: Estimated relative effects (e.g., log-odds ratios, risk
    differences).

  - `variances`: Variance-covariance matrix of the relative effects.

  - `contrast_ci`: Confidence intervals for the relative effects.

- absolute:

  A list containing statistics for absolute treatment outcomes:

  - `means`: Estimated absolute outcomes (e.g., probabilities, mean
    response).

  - `variances`: Variance-covariance matrix of the absolute outcomes.

  - `ci`: Confidence intervals for the absolute outcomes.

- CI:

  The confidence level used (e.g., 0.95).

- ref_trt:

  The name of the reference treatment.

- scale:

  The scale of the outcome (e.g., "log odds", "probability").

- model:

  A list containing details of the underlying statistical model.
  Contents vary by strategy:

  - `family`: The error distribution and link function.

  - `fit`: The underlying model object (e.g., for STC, G-Comp ML, or
    Bayesian G-Comp).

  - `weights`, `ESS`: (MAIC only) The estimated weights and Effective
    Sample Size.

  - `stan_args`: (Bayesian G-Comp, MIM) Arguments passed to Stan.

  - `rho`: (G-Comp ML, MIM, Bayesian G-Comp) Correlation coefficient.

  - `N`: (G-Comp ML, MIM, Bayesian G-Comp) Number of iterations.

  - `nu`, `hats.v`, `M`: (MIM only) Imputation parameters and matrices.

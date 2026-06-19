# Calculate Naive (Unadjusted) Bucher Indirect Comparison

Calculates the unadjusted indirect comparison between two active
treatments via a common comparator using the standard Bucher method.
Estimates are calculated on the link scale (e.g., log-odds) and
back-transformed to the requested scale.

## Usage

``` r
calc_bucher_naive(
  ipd_trial,
  ald_trial,
  outcome_model,
  family,
  ref_trt = "C",
  scale = NULL
)
```

## Arguments

- ipd_trial:

  Data frame containing the Individual Patient Data.

- ald_trial:

  Data frame containing the Aggregate Level Data.

- outcome_model:

  Formula specifying the outcome model (used to identify variables).

- family:

  A family object (e.g., `binomial(link="logit")`).

- ref_trt:

  Character string identifying the common comparator treatment.

- scale:

  Character string specifying the desired output scale (e.g.,
  "risk_difference").

## Value

A list containing the naive point estimate, variance, and standard
error.

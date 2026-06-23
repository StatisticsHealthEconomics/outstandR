<div id="main" class="col-md-9" role="main">

# Package index

<div class="section level2">

## Deprecated

<div class="section-desc">

Functions that will be removed in future versions.

</div>

</div>

<div class="section level2">

-   `strategy_stc()` **\[deprecated\]** : Simulated treatment comparison
    (STC)

</div>

<div class="section level2">

## Main

<div class="section-desc">

Main analysis functions.

</div>

</div>

<div class="section level2">

-   `outstandR-class` : outstandR class
-   `outstandR()` : Calculate the difference between treatments using
    all evidence
-   `strategy_maic()` `strategy_gcomp_ml()` `strategy_gcomp_bayes()`
    `strategy_mim()` `new_strategy()` : New strategy objects

</div>

<div class="section level2">

## Constituent statistics

<div class="section-desc">

Helper functions used by the main functions.

</div>

</div>

<div class="section level2">

-   `marginal_treatment_effect()` : Marginal treatment effect from
    reported event counts
-   `marginal_variance()` : Marginal effect variance using the delta
    method
-   `calculate_trial_variance()` : Calculate trial variance
-   `calculate_ate()` : Calculate Average Treatment Effect
-   `calculate_trial_mean()` : Calculate Trial Mean Wrapper
-   `calculate_trial_mean_binary()` : Calculate Trial Mean Binary Data
-   `calculate_trial_mean_continuous()` : Calculate Trial Mean
    Continuous Data
-   `calculate_trial_mean_count()` : Calculate Trial Mean Count Data
-   `calculate_trial_variance_binary()` : Calculate trial variance
    binary
-   `calculate_trial_variance_continuous()` : Calculate trial variance
    continuous
-   `calculate_trial_variance_count()` : Calculate trial variance count
-   `estimate_var_sandwich()` : Estimate Variance Sandwich Estimator
-   `calc_ALD_stats()` : Aggregate-level data mean and variance
    statistics
-   `calc_bucher_naive()` : Calculate Naive (Unadjusted) Bucher Indirect
    Comparison
-   `calc_gcomp_bayes()` : Bayesian G-computation using Stan
-   `calc_gcomp_ml()` : G-computation Maximum Likelihood Bootstrap
-   `calc_IPD_stats()` : Calculate individual-level patient data
    statistics
-   `get_treatment_effect()` : Get treatment effect scale corresponding
    to a link function

</div>

<div class="section level2">

## Miscellaneous

</div>

<div class="section level2">

-   `print(<outstandR>)` : Print a Summary of a outstandR Object
-   `plot(<outstandR>)` : Default Plot Method for outstandR Objects
-   `summary(<outstandR>)` `print(<summary.outstandR>)` : Summary method
    for outstandR
-   `reshape_ald_to_long()` : Convert aggregate data from wide to long
    format
-   `reshape_ald_to_wide()` : Convert aggregate data from long to wide
    format

</div>

<div class="section level2">

## Data

</div>

<div class="section level2">

-   `AC_IPD_binY_contX` : Individual-level patient data for binary
    outcome, continuous covariates
-   `BC_ALD_binY_contX` : Aggregate level patient data for binary
    outcome, continuous covariates
-   `AC_IPD_contY_mixedX` : Individual-level patient data for continuous
    outcome, mixed covariates
-   `AC_IPD_countY_contX` : Individual-level patient data for count
    outcome, continuous covariates
-   `BC_ALD_contY_mixedX` : Aggregate level patient data for continuous
    outcome, mixed covariates
-   `BC_ALD_countY_contX` : Aggregate level patient data for count
    outcome, continuous covariates

</div>

<div class="section level2">

## Classes

</div>

<div class="section level2">

-   `strategy-class` : Strategy class and subclasses

</div>

</div>

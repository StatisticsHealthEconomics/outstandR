# Package index

## Deprecated

Functions that will be removed in future versions.

- [`strategy_maic()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  [`strategy_stc()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  [`strategy_gcomp_ml()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  [`strategy_gcomp_bayes()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  [`strategy_mim()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  [`new_strategy()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  : New strategy objects

## Main

Main analysis functions.

- [`outstandR-class`](https://StatisticsHealthEconomics.github.io/outstandR/reference/outstandR-class.md)
  : outstandR class
- [`outstandR()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/outstandR.md)
  : Calculate the difference between treatments using all evidence
- [`strategy_maic()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  [`strategy_stc()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  [`strategy_gcomp_ml()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  [`strategy_gcomp_bayes()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  [`strategy_mim()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  [`new_strategy()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  : New strategy objects

## Constituent statistics

Helper functions used by the main functions.

- [`marginal_treatment_effect()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/marginal_treatment_effect.md)
  : Marginal treatment effect from reported event counts
- [`marginal_variance()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/marginal_variance.md)
  : Marginal effect variance using the delta method
- [`calculate_trial_variance()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calculate_trial_variance.md)
  : Calculate trial variance
- [`calculate_ate()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calculate_ate.md)
  : Calculate Average Treatment Effect
- [`calculate_trial_mean()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calculate_trial_mean.md)
  : Calculate Trial Mean Wrapper
- [`calculate_trial_mean_binary()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calculate_trial_mean_binary.md)
  : Calculate Trial Mean Binary Data
- [`calculate_trial_mean_continuous()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calculate_trial_mean_continuous.md)
  : Calculate Trial Mean Continuous Data
- [`calculate_trial_mean_count()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calculate_trial_mean_count.md)
  : Calculate Trial Mean Count Data
- [`calculate_trial_variance_binary()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calculate_trial_variance_binary.md)
  : Calculate trial variance binary
- [`calculate_trial_variance_continuous()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calculate_trial_variance_continuous.md)
  : Calculate trial variance continuous
- [`calculate_trial_variance_count()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calculate_trial_variance_count.md)
  : Calculate trial variance count
- [`estimate_var_sandwich()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/estimate_var_sandwich.md)
  : Estimate Variance Sandwich Estimator
- [`calc_ALD_stats()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calc_ALD_stats.md)
  : Aggregate-level data mean and variance statistics
- [`calc_gcomp_bayes()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calc_gcomp_bayes.md)
  : Bayesian G-computation using Stan
- [`calc_gcomp_ml()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calc_gcomp_ml.md)
  : G-computation Maximum Likelihood Bootstrap
- [`calc_IPD_stats()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/calc_IPD_stats.md)
  : Calculate individual-level patient data statistics
- [`get_treatment_effect()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/get_treatment_effect.md)
  : Get treatment effect scale corresponding to a link function

## Miscellaneous

- [`print(`*`<outstandR>`*`)`](https://StatisticsHealthEconomics.github.io/outstandR/reference/print.outstandR.md)
  : Print a Summary of a outstandR Object
- [`plot(`*`<outstandR>`*`)`](https://StatisticsHealthEconomics.github.io/outstandR/reference/plot.outstandR.md)
  : Default Plot Method for outstandR Objects
- [`summary(`*`<outstandR>`*`)`](https://StatisticsHealthEconomics.github.io/outstandR/reference/summary.outstandR.md)
  [`print(`*`<summary.outstandR>`*`)`](https://StatisticsHealthEconomics.github.io/outstandR/reference/summary.outstandR.md)
  : Summary method for outstandR
- [`reshape_ald_to_long()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/reshape_ald_to_long.md)
  : Convert aggregate data from wide to long format
- [`reshape_ald_to_wide()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/reshape_ald_to_wide.md)
  : Convert aggregate data from long to wide format

## Data

- [`AC_IPD_binY_contX`](https://StatisticsHealthEconomics.github.io/outstandR/reference/AC_IPD_binY_contX.md)
  : Individual-level patient data for binary outcome, continuous
  covariates
- [`BC_ALD_binY_contX`](https://StatisticsHealthEconomics.github.io/outstandR/reference/BC_ALD_binY_contX.md)
  : Aggregate level patient data for binary outcome, continuous
  covariates
- [`AC_IPD_contY_mixedX`](https://StatisticsHealthEconomics.github.io/outstandR/reference/AC_IPD_contY_mixedX.md)
  : Individual-level patient data for continuous outcome, mixed
  covariates
- [`AC_IPD_countY_contX`](https://StatisticsHealthEconomics.github.io/outstandR/reference/AC_IPD_countY_contX.md)
  : Individual-level patient data for count outcome, continuous
  covariates
- [`BC_ALD_contY_mixedX`](https://StatisticsHealthEconomics.github.io/outstandR/reference/BC_ALD_contY_mixedX.md)
  : Aggregate level patient data for continuous outcome, mixed
  covariates
- [`BC_ALD_countY_contX`](https://StatisticsHealthEconomics.github.io/outstandR/reference/BC_ALD_countY_contX.md)
  : Aggregate level patient data for count outcome, continuous
  covariates

## Classes

- [`strategy-class`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy-class.md)
  : Strategy class and subclasses

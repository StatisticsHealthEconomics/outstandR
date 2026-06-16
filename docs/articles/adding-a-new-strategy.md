# Adding a New Strategy to \`outstandR\`

The `outstandR` package uses R’s S3 object-oriented system to manage
different statistical approaches (strategies) for population adjustment.
This modular design makes it straightforward to implement new methods
without rewriting the package’s core infrastructure.

To add a new statistical method (e.g., `custom_method`), you need to
complete the following 5 steps:

1.  Create a strategy constructor function.
2.  Document the new class fields.
3.  Write the core calculation function.
4.  Register the S3 dispatch method via the package factory.
5.  Create a print method.

------------------------------------------------------------------------

## Step 1: Create the Strategy Constructor

**File:** `R/strategy_.R`

First, create a user-facing constructor function that validates user
inputs and instantiates the strategy object. By convention, this
function should be named `strategy_<method_name>()`.

It must call the internal
[`new_strategy()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
function, passing the class name and the validated arguments.

``` r

#' @rdname strategy
#' 
#' @section Custom Method:
#' Describe the statistical theory behind the new method here.
#' 
#' @param trt_var Treatment variable name
#' @param my_custom_arg Description of your custom argument
#' 
#' @return `custom_method` class object
#' @export
strategy_custom_method <- function(formula = NULL,
                                   family = gaussian(link = "identity"),
                                   trt_var = NULL,
                                   my_custom_arg = 100) {
  
  # 1. Validate inputs using built-in helpers
  check_formula(formula, trt_var)
  check_family(family)
  
  if (my_custom_arg <= 0) {
    stop("my_custom_arg must be greater than 0.")
  }
  
  # 2. Bundle arguments
  args <- list(formula = formula,
               family = family,
               trt_var = get_treatment_name(formula, trt_var),
               my_custom_arg = my_custom_arg)
  
  # 3. Create object
  do.call(new_strategy, c(strategy = "custom_method", args))
}
```

## Step 2: Document the Class Fields

**File:** `R/strategy-class.R`

Update the `strategy-class` documentation to inform users about the
internal fields your new subclass carries. Add a new `\item` block under
the details section:

``` r

#'   \item{custom_method subclass}{Additional fields for the Custom Method:
#'     \itemize{
#'       \item `my_custom_arg`: Describes what this internal field does.
#'     }
#'   }
```

## Step 3: Write the Core Calculation Function

**File:** Create a new file (e.g., `R/custom_method.R`) or add to
existing calculation files.

Write the function that actually performs the maths. It will receive two
arguments: `strategy` (the object you created in Step 1) and
`analysis_params` (a list containing the pre-processed `ipd`, `ald`,
treatment labels, etc.).

**Crucial Requirement:** Your function MUST return a list containing
exactly two elements: `means` and `model`.

- `means$A`: Estimates for the comparator treatment group (can be a
  scalar, or a vector of bootstrap/MCMC draws).
- `means$C`: Estimates for the reference treatment group (can be a
  scalar, or a vector of bootstrap/MCMC draws).
- `model`: A list of any underlying models or diagnostics (e.g.,
  weights, glm fit objects) you want returned to the user.

``` r

#' Calculate Custom Method
#' @keywords internal
calc_custom_method <- function(strategy, analysis_params) {
  
  # Extract data and parameters
  ipd <- analysis_params$ipd
  ald <- analysis_params$ald
  formula <- strategy$formula
  custom_val <- strategy$my_custom_arg
  
  # ... [PERFORM YOUR STATISTICAL MODELLING HERE] ...
  # e.g., fit a model, run bootstraps, calculate weights
  
  # Construct the strict return object
  list(
    means = list(
      A = estimated_means_for_A, 
      C = estimated_means_for_C
    ),
    model = list(
      fit = my_fitted_model,
      custom_val = custom_val
    )
  )
}
```

## Step 4: Register the S3 Dispatch Method

**File:** `R/calc_IPD_stats.R`

To plug your core function into `outstandR`’s top-level execution, you
must register it as an S3 method for `calc_IPD_stats`.

Because `outstandR` handles the relative treatment effect
transformations (ATE) and variance estimation automatically, you simply
pass your core function to the
[`IPD_stat_factory()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/IPD_stat_factory.md).

``` r

#' @rdname calc_IPD_stats
#' @section Custom Method statistics:
#' Description of how the custom method calculates statistics.
#' @export
calc_IPD_stats.custom_method <- IPD_stat_factory(calc_custom_method)
```

## Step 5: Add a Print Method

**File:** `R/print-strategies.R`

Finally, create a nice console printout for when a user types the name
of their instantiated strategy. Extend the generic `print.strategy`
using [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) to print
the common fields, then print your specific parameters.

``` r

#' @export
#' @method print custom_method
print.custom_method <- function(x, ...) {
  # Print the standard strategy info (formula, family)
  NextMethod() 
  
  # Print custom subclass info
  cat(pillar::style_subtle("  Parameters:"), "\n")
  cat("    Custom Argument:     ", x$my_custom_arg, "\n")
  
  invisible(x)
}
```

## (Optional) Step 6: Passing Random Seeds

If your new method utilizes a random number generator (e.g., MCMC
sampling, pseudo-population simulation), you might need to catch the
top-level `seed` argument.

**File:** `R/outstandR.R`

Create an `add_seed` S3 method for your strategy:

``` r

#' @export
add_seed.custom_method <- function(strategy, fn_args, seed) {
  if (!is.null(seed)) {
    fn_args$seed <- seed
  }
  return(fn_args)
}
```

Here is a section you can append to your vignette or developer guide. It
provides a copy-pasteable utility script that future contributors can
use to automate the boilerplate setup for a new strategy.

------------------------------------------------------------------------

## Step 7: Automating Boilerplate with a Scaffolding Script (Optional)

If you or your contributors find yourselves adding new strategies
frequently, you can automate the creation of the boilerplate files using
the `usethis` package.

We recommend saving the following helper function in a script file
(e.g., `scripts/scaffold_strategy.R`). This script will automatically
create the necessary calculation file and open the registry files for
you to edit.

### How to Use It

When you are ready to build a new method, simply source the script and
run the function with the name of your new method:

``` r

# Load the helper function
source("scripts/scaffold_strategy.R")

# Scaffold a new method called "gcomp_rf" (e.g., Random Forest G-computation)
scaffold_strategy("gcomp_rf")
```

This will instantly set up your workspace, leaving you to focus solely
on the statistical implementation rather than package wiring!

Here is a section you can add to your Markdown document (perhaps at the
end, right before or after the scaffolding script) that outlines best
practices for writing a robust, user-friendly strategy in `outstandR`.

------------------------------------------------------------------------

## Best Practices for Writing a Good Strategy

When contributing a new statistical approach to `outstandR`, following
these guidelines ensures your method is robust, user-friendly, and
maintainable.

### 1. Fail Fast and Validate Early

Catch user errors in the constructor function
(`strategy_custom_method()`) before any heavy computation begins.

- Use the package’s built-in
  [`check_formula()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/check_formula.md)
  and `check_family()` helpers.
- If your method requires specific types of data (e.g., only binary
  outcomes), enforce it immediately:

``` r

if (family$family != "binomial") {
  cli::cli_abort("Custom Method currently only supports binary outcomes (`family = binomial()`).")
}
```

### 2. Return a Rich `model` Object

The core calculation function must return a list with `means` and
`model`. While `means` is strictly structured, the `model` element is a
free-form list. Use this to give the user as much diagnostic power as
possible.

- **Good things to include:** The raw fitted model objects (e.g., `glm`
  or `stanfit`), calculated weights, convergence status flags, and
  effective sample sizes (ESS).
- This allows users to extract the `model` from the final `outstandR`
  object and check diagnostics themselves.

### 3. Handle Convergence Issues Gracefully

Many population adjustment methods rely on optimization or iterative
fitting (like GLMs or MAIC). If an internal model fails to converge,
your function shouldn’t crash the entire R session—especially during a
bootstrap loop.

- Wrap volatile model fitting in
  [`tryCatch()`](https://rdrr.io/r/base/conditions.html).
- If a fit fails, return `NA` for the means and pass a warning to the
  user.

``` r

fit <- tryCatch({
  glm(formula, data = ipd, family = family)
}, warning = function(w) {
  # Handle warnings (e.g., probabilities numerically 0 or 1)
  glm(formula, data = ipd, family = family)
}, error = function(e) {
  NULL # Return NULL on hard failure
})

if (is.null(fit)) return(list(means = list(A = NA, C = NA), model = list(converged = FALSE)))
```

### 4. Respect the `verbose` Flag

`outstandR` provides a global `verbose` argument so users can silence
console output during heavy simulations.

- Extract `verbose <- isTRUE(analysis_params$verbose)` at the top of
  your calculation function.
- Use the `cli` package to print progress updates or warn users about
  computationally expensive steps, but *only* if `verbose` is true.

``` r

if (verbose) {
  cli::cli_h2("Custom Method Execution")
  cli::cli_alert_info("Optimizing weights. This may take a moment...")
}
```

### 5. Write Unit Tests for Your Method

Whenever you add a new strategy, add a corresponding test file in
`tests/testthat/` (e.g., `test-custom_method.R`).

- Test that the strategy constructor throws errors on bad inputs.
- Use the package’s built-in dummy data (`AC_IPD_binY_contX`, etc.) to
  run a fast, low-iteration version of your method to ensure the output
  list matches the required `list(means = ..., model = ...)` structure.

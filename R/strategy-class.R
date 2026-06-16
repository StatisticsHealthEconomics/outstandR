#' Strategy class and subclasses
#' 
#' The `strategy` class is a virtual class that defines the statistical approach 
#' for population adjustment in indirect treatment comparisons
#' These objects are constructors that validate hyperparameters and encapsulate 
#' modelling settings before execution by [outstandR()]
#' 
#' @rdname strategy-class
#' @name strategy-class
#' 
#' @details Objects of class `strategy` have a common structure but carry 
#' different subclasses to trigger specific S3 method dispatch
#' 
#' \describe{
#'   \item{General fields}{Shared by all strategies:
#'     \itemize{
#'       \item `formula`: The linear regression formula for the outcome model
#'       \item `family`: A base R family object specifying the distribution and link
#'       \item `trt_var`: The name of the treatment variable.
#'     }
#'   }
#'   \item{maic subclass}{Additional fields for Matching-Adjusted Indirect Comparison:
#'     \itemize{
#'       \item `n_boot`: Number of bootstrap resamples for variance estimation.
#'     }
#'   }
#'   \item{stc subclass}{Additional fields for Simulated Treatment Comparison:
#'     \itemize{
#'       \item `N`: Synthetic sample size for the target population.
#'     }
#'   }
#'   \item{gcomp_ml subclass}{Additional fields for Maximum Likelihood G-computation:
#'     \itemize{
#'       \item `rho`: Named square matrix of covariate correlations.
#'       \item `marginal_distns`: Names of the marginal distributions for covariates.
#'       \item `marginal_params`: Parameters for the marginal distributions.
#'       \item `N`: Synthetic sample size for the pseudo-population.
#'       \item `n_boot`: Number of bootstrap resamples.
#'     }
#'   }
#'   \item{gcomp_bayes subclass}{Additional fields for Bayesian G-computation:
#'     \itemize{
#'       \item `rho`, `marginal_distns`, `marginal_params`, `N`: Same as `gcomp_ml`.
#'       \item `...`: Additional arguments passed to the Stan engine via `rstanarm::stan_glm()`.
#'     }
#'   }
#'   \item{mim subclass}{Additional fields for Multiple Imputation Marginalization:
#'     \itemize{
#'       \item `rho`: Correlation matrix.
#'       \item `N`: Number of iterations/simulated individuals.
#'     }
#'   }
#' }
#' 
NULL

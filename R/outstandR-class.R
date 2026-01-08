#' outstandR class
#' 
#' The `outstandR` class contains the results from running a
#' model with the function [outstandR()].
#' 
#' @rdname outstandR-class
#' @name outstandR-class
#' 
#' @details Objects of class `outstandR` have the following
#' 
#' \describe{
#'   \item{contrasts}{A list containing statistics for relative treatment effects:
#'     \itemize{
#'       \item `means`: Estimated relative effects (e.g., log-odds ratios, risk differences).
#'       \item `variances`: Variance-covariance matrix of the relative effects.
#'       \item `contrast_ci`: Confidence intervals for the relative effects.
#'     }
#'   }
#'   \item{absolute}{A list containing statistics for absolute treatment outcomes:
#'     \itemize{
#'       \item `means`: Estimated absolute outcomes (e.g., probabilities, mean response).
#'       \item `variances`: Variance-covariance matrix of the absolute outcomes.
#'       \item `ci`: Confidence intervals for the absolute outcomes.
#'     }
#'   }
#'   \item{CI}{The confidence level used (e.g., 0.95).}
#'   \item{ref_trt}{The name of the reference treatment.}
#'   \item{scale}{The scale of the outcome (e.g., "log odds", "probability").}
#'   \item{model}{A list containing details of the underlying statistical model. 
#'     Contents vary by strategy:
#'     \itemize{
#'       \item `family`: The error distribution and link function.
#'       \item `fit`: The underlying model object (e.g., for STC, G-Comp ML, or Bayesian G-Comp).
#'       \item `weights`, `ESS`: (MAIC only) The estimated weights and Effective Sample Size.
#'       \item `stan_args`: (Bayesian G-Comp, MIM) Arguments passed to Stan.
#'       \item `rho`: (G-Comp ML, MIM, Bayesian G-Comp) Correlation coefficient.
#'       \item `N`: (G-Comp ML, MIM, Bayesian G-Comp) Number of iterations.
#'       \item `nu`, `hats.v`, `M`: (MIM only) Imputation parameters and matrices.
#'     }
#'   }
#' }
#'
NULL


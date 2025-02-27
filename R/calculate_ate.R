
#' Calculate Average Treatment Effect
#'
#' @param mean_A,mean_C Mean of the outcome for the treatment and control
#' @param family family object of the model
#'
#' @returns ATE
#' @export
#'
calculate_ate <- function(mean_A, mean_C, family) {

  link <- family$link
  
  if (link == "logit") {
    ate <- qlogis(mean_A) - qlogis(mean_C)  # log-odds ratio
  } else if (link == "identity") {
    ate <- mean_A - mean_C         # mean difference
  } else if (link == "probit") {
    ate <- qnorm(mean_A) - qnorm(mean_C)    # probit link
  } else if (link == "cloglog") {
    ate <- log(-log(1 - mean_A)) - log(-log(1 - mean_C))
  } else if (link == "log") {      # Poisson log link
    ate <- log(mean_A) - log(mean_C)        # log ratio
  } else {
    stop("Unsupported link function. Choose from 'logit', 'identity', 'probit', 'cloglog', or 'log'.")
  }
  
  ate
}


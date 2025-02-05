
#' Calculate ATE
#'
#' @param ppv model prediction samples
#' @param family family object of the model
#'
#' @returns ATE
#' @export
#'
calculate_ate <- function(mean_A, mean_C, family) {

  link <- family$link
  
  if (link == "logit") {
    ate <- qlogis(mean_A) - qlogis(mean_C)
  } else if (link == "identity") {
    ate <- mean_A - mean_C
  } else if (link == "probit") {
    ate <- qnorm(mean_A) - qnorm(mean_C)
  } else if (link == "cloglog") {
    ate <- log(-log(1 - mean_A)) - log(-log(1 - mean_C))
  } else if (link == "log") {  # Poisson log link
    ate <- log(mean_A) - log(mean_C)
  } else {
    stop("Unsupported link function. Choose from 'logit', 'identity', 'probit', 'cloglog', or 'log'.")
  }
  
  ate
}


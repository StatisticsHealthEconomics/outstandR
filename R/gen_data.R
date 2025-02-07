
#' Generate simulated datasets of IPD covariates and binary outcome for a trial
#' 
#' @param N Total number of patients
#' @param b_trt `b` coefficient for active treatment vs. common comparator 
#' @param b_X `b` coefficients for each prognostic variable `X` 
#' @param b_EM `b` coefficients effect modifiers
#' @param b_0 Intercept coefficient
#' @param meanX Mean of each normally-distributed covariate `X` 
#' @param sdX Standard deviation of each covariate `X` 
#' @param event_rate Event rate 
#' @param corX Covariate correlation coefficient of `X` 
#' @param allocation Allocation to active treatment as proportion of total; 0 to 1
#' @return Data frame of `X`, `trt` and `y`
#' 
#' @importFrom MASS mvrnorm
#' @export
#' 
#' @examples
#' 
#' \dontrun{
#' x <- gen_data(
#'  N = 100,
#'  b_trt = log(0.17),
#'  b_X = -log(0.5),
#'  b_EM = -log(0.67),
#'  b_0 = -0.62,
#'  meanX = 0.6,
#'  sdX = 0.4,
#'  event_rate = 0.35, 
#'  corX = 0.2,
#'  allocation = 2/3) 
#' 
#' head(x)
#' }
gen_data <- function(N, b_trt, b_X, b_EM, b_0,
                     meanX_EM, sdX_EM, 
                     meanX, sdX, 
                     corX, allocation,
                     family = binomial("logit")) {
  ##TODO: what does event_rate do?

  nX <- length(meanX)
  nX_EM <- length(meanX_EM)
  n_c <- nX + nX_EM
    
  rho <- matrix(corX, nrow=n_c, ncol=n_c) # set correlation matrix
  diag(rho) <- rep(1, n_c)
  N_active <- round(N*allocation)  # number of patients under active treatment
  N_control <- N - N_active        # number of patients under control
  sd.vec <- c(sdX, sdX_EM)         # vector of standard deviations
  
  cov.mat <- cor2cov(rho, sd.vec)  # covariance matrix
  
  # simulate correlated continuous covariates using multivariate normal
  # patients under active treatment
  X_active <- 
    as.data.frame(
      MASS::mvrnorm(n = N_active,
                    mu = c(meanX, meanX_EM),
                    Sigma = cov.mat))
  
  # patients under control treatment
  X_control <-
    as.data.frame(
      MASS::mvrnorm(n = N_control,
                    mu = c(meanX, meanX_EM),
                    Sigma = cov.mat))  
  # all patients
  X <- rbind(X_active, X_control)
  colnames(X) <- paste0("X", 1:ncol(X))
  
  # treatment assignment (1: active; 0: control)
  trt <- c(rep(1, N_active),
           rep(0, N_control))
  
  Xnames <- colnames(X)
  PF_names <- Xnames[1:nX] 
  EM_names <- Xnames[nX + (1:nX_EM)]
  
  design_mat <- X |> 
    mutate(trt = trt,
           X0 = 1) |> 
    # add interaction terms
    mutate(across(all_of(EM_names), ~ . * trt)) |> 
    relocate(X0) |>
    relocate(trt, .after = last_col())
  
  # generate outcomes using regression
  # linear predictor
  betas <- c(b_0, rep(b_X, nX), rep(b_EM, nX_EM), b_trt)
  LP <- as.matrix(design_mat) %*% betas
  
  if (family$family == "binomial") {
    yprob <- family$linkinv(LP)           # binary outcome probability
    y <- rbinom(n=N, size=1, prob=yprob)  # binary outcome
  } else if (family$family == "gaussian") {
    y <- rnorm(n=N, mean=LP, sd=1)        # continuous outcome
  } else if (family$family == "poisson") {
    yrate <- family$linkinv(LP)
    y <- rpois(n=N, lambda = yrate)        # counts outcome
  }
  
  as.data.frame(cbind(X, trt, y))
}

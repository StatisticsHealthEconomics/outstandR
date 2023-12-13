
#' Generate simulated datasets for index and comparator trials
#' 
#' @param N Size
#' @param b_trt b coefficient for treatment 
#' @param b_X b coefficients for X 
#' @param b_EM b coefficients effect modifiers
#' @param b_0 Intercept coefficient
#' @param meanX mean X 
#' @param sdX standard deviation X 
#' @param event_rate Event rate 
#' @param corX Correlation of X 
#' @param allocation Allocation
#' @return Data frame of `X`, `trt` and `y`
#' 
#' @keywords internal
#' 
gen_data <- function(N, b_trt, b_X, b_EM, b_0,
                     meanX, sdX, event_rate, 
                     corX, allocation) {
  # 4 baseline covariates
  rho <- matrix(corX, nrow=4, ncol=4) # set correlation matrix
  diag(rho) <- rep(1, 4)
  N_active <- round(N*allocation)  # number of patients under active treatment
  N_control <- N - N_active        # number of patients under control
  sd.vec <- rep(sdX, 4)            # vector of standard deviations
  
  cov.mat <- cor2cov(rho, sd.vec)  # covariance matrix
  # simulate correlated continuous covariates using multivariate normal
  # patients under active treatment
  X_active <- 
    as.data.frame(
      MASS::mvrnorm(n = N_active,
                    mu = rep(meanX, 4),
                    Sigma = cov.mat))
  # patients under control treatment
  X_control <-
    as.data.frame(
      MASS::mvrnorm(n = N_control,
                    mu = rep(meanX, 4),
                    Sigma = cov.mat))  
  # all patients
  X <- rbind(X_active, X_control)
  colnames(X) <- c("X1","X2","X3","X4") 
  # treatment assignment (1: active; 0: control)
  trt <- c(rep(1,N_active), rep(0,N_control))
  
  # generate binary outcomes using logistic regression
  # linear predictor
  LP <-
    b_0 + b_X*X$X1 + b_X*X$X2 +
    b_X*X$X3 + b_X*X$X4 + b_trt*trt +
    b_EM*X$X1*trt + b_EM*X$X2*trt
  
  yprob <- 1/(1 + exp(-LP))               # binary outcome probability
  y <- rbinom(n=N, size=1, prob=yprob)  # binary outcome
  
  as.data.frame(cbind(X, trt, y))
}

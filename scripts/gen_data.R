# This file specifies the simulation setup and generates the simulation study
# data

# setwd("C:/Users/Antonio/Desktop/Gcomp_indirect_comparisons_simstudy") 

rm(list=ls())

# Load packages
# package for data manipulation
if(!require("dplyr")) {install.packages("dplyr"); library(dplyr)}
# package to sample/simulate the covariates from a multivariate normal
if(!require("MASS")) {install.packages("MASS"); library(MASS)}

set.seed(555) # set random seed for reproducibility

# Define simulation study parameters

N_sim <- 2000 # number of Monte Carlo replicates
allocation <- 2/3 # active treatment vs. placebo allocation ratio (2:1)

N_AC <- c(200,400,600) # number of subjects in the AC trial
N_BC <- 600 # number of subjects in the BC trial 
b_trt <- log(0.17) # conditional effect of active treatment vs. common comparator
b_X <- -log(0.5) # conditional effect of each prognostic variable
b_EM <- -log(0.67) # conditional interaction effect of each effect modifier
event_rate <- 0.35
meanX_AC <- c(0.45,0.3,0.15) # mean of each normally-distributed covariate in AC trial
meanX_BC <- 0.6 # mean of each normally-distributed covariate in BC
sdX <- 0.4 # standard deviation of each covariate (same for AC and BC)
corX <- 0.2 # covariate correlation coefficient  

# parameter combinations for each scenario
param.combinations <- expand.grid(N_AC=N_AC, meanX_AC=meanX_AC)
pc <- param.combinations

scenarios <- nrow(pc) # number of simulation scenarios
save(pc, N_sim, allocation, file="binary_settings.RData")

# compute intercept based on baseline event rate (for patients under C without covariates)
optim.function <- function(param, y_prob) {
  b_0 <- param
  fit <- sum(((1 / (1 + exp(-b_0))) - y_prob)^2)
  return(fit)
}

b_0 <- optim(par=0,fn=optim.function,y_prob=event_rate,
             method="Brent",lower=-2,upper=2)$par

# generate simulated datasets for index and comparator trials
gen.data <- function(N, b_trt, b_X, b_EM, b_0, meanX, sdX, event_rate, 
                     corX, allocation) {
  # 4 baseline covariates
  rho <- matrix(corX, nrow=4, ncol=4) # set correlation matrix
  diag(rho) <- rep(1, 4)
  N_active <- round(N*allocation) # number of patients under active treatment
  N_control <- N - N_active # number of patients under control
  sd.vec <-rep(sdX, 4) # vector of standard deviations
  cor2cov <- function(R, S) {
    # function to compute covariance matrix from correlation matrix R and vector
    # of standard deviations S. covariance matrix required as input for mvrnorm
    sweep(sweep(R, 1, S, "*"),2,S,"*")
  }
  cov.mat <- cor2cov(rho, sd.vec) # covariance matrix
  # simulate correlated continuous covariates using multivariate normal
  # patients under active treatment
  X_active <- as.data.frame(MASS::mvrnorm(n=N_active,mu=rep(meanX,4), Sigma=cov.mat))
  # patients under control treatment
  X_control <- as.data.frame(MASS::mvrnorm(n=N_control,mu=rep(meanX,4), Sigma=cov.mat))  
  # all patients
  X <- rbind(X_active, X_control)
  colnames(X) <- c("X1","X2","X3","X4") 
  # treatment assignment (1: active; 0: control)
  trt <- c(rep(1,N_active),rep(0,N_control)) 
  # generate binary outcomes using logistic regression
  LP <- b_0 + b_X*X$X1 + b_X*X$X2 + b_X*X$X3 + b_X*X$X4 + b_trt*trt + 
        b_EM*X$X1*trt + b_EM*X$X2*trt # linear predictor
  yprob <- 1/ (1+exp(-LP)) # binary outcome probability
  y <- rbinom(n=N, size=1, prob=yprob) # binary outcome
  return(as.data.frame(cbind(X, trt, y)))
}

for (i in 1:scenarios) {
  print(i)
  # simulate IPD covariates and outcome for A vs. C trial (S=1)
  IPD.AC <- replicate(n=N_sim, expr=gen.data(pc$N_AC[i], b_trt, b_X, b_EM, 
                                             b_0, pc$meanX_AC[i], sdX, event_rate, 
                                             corX, allocation),
                      simplify=FALSE)
  # simulate IPD covariates and outcome for B vs. C trial (S=2)
  IPD.BC <- replicate(n=N_sim, expr=gen.data(N_BC, b_trt, b_X, b_EM, 
                                             b_0, meanX_BC, sdX, event_rate, 
                                             corX, allocation),
                      simplify=FALSE)
  # Summarize BC IPD as ALD
  ALD.BC <- lapply(1:N_sim, function(j) {
    as.data.frame(cbind(
      # aggregate the data for the BC trial 
      summarise(IPD.BC[[j]], mean.X1=mean(X1), mean.X2=mean(X2), mean.X3=mean(X3),
                mean.X4=mean(X4), sd.X1=sd(X1), sd.X2=sd(X2), sd.X3=sd(X3), sd.X4=sd(X4)),
      # summarize the outcomes for the BC trial (treatment B)
      filter(IPD.BC[[j]], trt == 1) %>%
        summarise(y.B.sum=sum(y), y.B.bar = mean(y), N.B = n()),
      # summarize the outcomes for the BC trial (treatment C)
      filter(IPD.BC[[j]], trt == 0) %>%
        summarise(y.C.sum=sum(y), y.C.bar = mean(y), N.C = n())))    
  } )
  file.id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i]) 
  save(IPD.AC, file=paste0("Data/IPD_AC_", file.id, ".RData"))
  save(IPD.BC, file=paste0("Data/IPD_BC_", file.id, ".RData"))
  save(ALD.BC, file=paste0("Data/ALD_BC_", file.id, ".RData"))  
}                       
                      
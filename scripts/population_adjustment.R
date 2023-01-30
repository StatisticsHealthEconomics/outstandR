# This file performs the indirect comparison methods on the simulated data

rm(list=ls())

# setwd("C:/Users/Antonio/Desktop/Gcomp_indirect_comparisons_simstudy") 
source('functions.R') # load MAIC function
load("binary_settings.RData") # load sim. study parameter combinations

# for detect cores
if(!require(parallel)) {install.packages("parallel"); library(parallel)}
# for parallel cluster
if(!require(doSNOW)) {install.packages("doSNOW"); library(doSNOW)}
# for fitting outcome regression and drawing predicted outcomes from posterior in Bayes G-comp
if(!require(rstanarm)) {install.packages("rstanarm"); library(rstanarm)}
# for non-parametric bootstrap in ML G-computation and MAIC
if(!require(boot)) {install.packages("boot"); library(boot)}
# for simulating BC trial covariates from multivariate Gaussian copula
if(!require(copula)) {install.packages("copula"); library(copula)}

set.seed(444) # set seed for reproducibility
scenarios <- nrow(pc) # number of simulation scenarios

# settings for the methods
resamples <- 1000 # Number of bootstrap resamples for MAIC and maximum-likelihood G-computation
N_star <- 1000 # size of simulated BC pseudo-population; high for small MC error
# MCMC info for Bayesian G-computation
n.chains <- 2 # number of Markov chains for MCMC
warmup <- 2000 # number of warmup aka discarded burn-in iterations per chain
iters <- 4000 # number of total iterations per chain (including warmup)

# simulated patient-level (AC) and aggregate-level (BC) datasets for all scenarios
IPD.AC.all <- vector(mode="list", scenarios)
ALD.BC.all <- vector(mode="list", scenarios)
# load data
for (i in 1:scenarios) {
  file.id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i]) 
  load(paste0("Data/IPD_AC_", file.id, ".RData")) # load AC patient-level data
  load(paste0("Data/ALD_BC_", file.id, ".RData")) # load BC aggregate-level data
  IPD.AC.all[[i]] <- IPD.AC
  ALD.BC.all[[i]] <- ALD.BC
}

replicates <- N_sim # number of Monte Carlo replicates per scenario

### Matching-adjusted indirect comparison (MAIC)
maic.wrapper <- function(data.AC, data.BC, resamples) { 
  # Inputs: data.AC - AC individual patient-level data; data.BC - BC aggregate-level data  
  # resamples - number of resamples for non-parametric bootstrap
  maic.boot <- function(data, indices) {
    dat <- data[indices,]
    x.EM <- dat[,c("X1","X2")] # AC effect modifiers 
    theta <- data.BC[c("mean.X1", "mean.X2")] # BC effect modifier means, assumed fixed
    # center the AC effect modifiers on the BC means
    x.EM$X1 <- x.EM$X1 - theta$mean.X1
    x.EM$X2 <- x.EM$X2 - theta$mean.X2
    # MAIC weights estimated using standard method of moments
    hat.w <- maic(X.EM=x.EM) # estimated weights
    # aess <- sum(hat.w)^2/sum(hat.w^2) # approximate effective sample size
    # fit weighted logistic regression model using glm
    outcome.fit <- glm(y~trt, family="quasibinomial", weights=hat.w, data=dat)
    # fitted treatment coefficient is marginal treatment effect for A vs. C
    hat.Delta.AC <- coef(outcome.fit)["trt"] 
    return(hat.Delta.AC)
  }
  # non-parametric bootstrap
  boot.object <- boot::boot(data=data.AC, statistic=maic.boot, R=resamples)
  # bootstrap mean of marginal A vs. C treatment effect estimate
  hat.Delta.AC <- mean(boot.object$t)
  # bootstrap variance of A vs. C treatment effect estimate   
  hat.var.Delta.AC <- var(boot.object$t)
  # B vs. C marginal treatment effect from reported event counts
  hat.Delta.BC <- with(data.BC, log(y.B.sum*(N.C-y.C.sum)/(y.C.sum*(N.B-y.B.sum))))
  # B vs. C marginal treatment effect variance using the delta method
  hat.var.Delta.BC <- with(data.BC, 1/y.C.sum+1/(N.C-y.C.sum)+1/y.B.sum+1/(N.B-y.B.sum))
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # A vs. B
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC  
  list(hat.Delta.AB, hat.var.Delta.AB)
}  

### Simulated treatment comparison (STC); conventional outcome regression approach
stc.wrapper <- function(data.AC, data.BC) {
  # Inputs: data.AC - AC individual patient-level data; data.BC - BC aggregate-level data
  # fit outcome regression model: regression of outcome on treatment and covariates
  # IPD effect modifiers centered at the mean BC values 
  # purely prognostic variables are included but not centered
  outcome.model <- glm(y~X3+X4+trt*I(X1-data.BC$mean.X1)+trt*I(X2-data.BC$mean.X2),
                       data=data.AC, family=binomial)
  # fitted treatment coefficient is relative conditional effect for A vs. C
  hat.Delta.AC <- coef(outcome.model)["trt"] 
  # estimated variance for A vs. C from model fit
  hat.var.Delta.AC <- vcov(outcome.model)["trt", "trt"] 
  # B vs. C marginal treatment effect estimated from reported event counts
  hat.Delta.BC <- with(data.BC, log(y.B.sum*(N.C-y.C.sum)/(y.C.sum*(N.B-y.B.sum))))
  # B vs. C marginal treatment effect variance using the delta method 
  hat.var.Delta.BC <- with(data.BC, 1/y.C.sum+1/(N.C-y.C.sum)+1/y.B.sum+1/(N.B-y.B.sum))
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # A vs. B
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC
  list(hat.Delta.AB, hat.var.Delta.AB)
}

### G-computation with maximum-likelihood estimation and bootstrapping
gcomp.ml.wrapper <- function(data.AC, data.BC, resamples, N_star) {
  # Inputs: data.AC - AC individual patient-level data; data.BC - BC aggregate-level data;
  # resamples - number of resamples for non-parametric bootstrap
  # N_star - size of simulated BC pseudo-population (high for small Monte Carlo error)
  # matrix of pairwise correlations between IPD covariates  
  rho <- cor(data.AC[,c("X1","X2","X3","X4")]) 
  #  covariate simulation for comparator trial using copula package
  cop <- normalCopula(param=c(rho[1,2],rho[1,3],rho[1,4],rho[2,3],rho[2,4],rho[3,4]), 
                      dim=4, dispstr="un") # AC IPD pairwise correlations
  # sample covariates from approximate joint distribution using copula
  mvd <- mvdc(copula=cop, margins=c("norm", "norm", "norm", "norm"), # Gaussian marginals
              # BC covariate means and standard deviations
              paramMargins=list(list(mean=data.BC$mean.X1, sd=data.BC$sd.X1),
                                list(mean=data.BC$mean.X2, sd=data.BC$sd.X2),       
                                list(mean=data.BC$mean.X3, sd=data.BC$sd.X3),
                                list(mean=data.BC$mean.X4, sd=data.BC$sd.X4)))
  # data frame of simulated covariates
  x_star <- as.data.frame(rMvdc(N_star, mvd))
  colnames(x_star) <- c("X1", "X2", "X3", "X4")
  gcomp.ml <- function(data, indices) {
    dat = data[indices,]
    # outcome logistic regression model fitted to original data using maximum likelihood
    outcome.model <- glm(y~X3+X4+trt*X1+trt*X2, data=dat, family=binomial)
    # counterfactual datasets
    data.trtA <- data.trtC <- x_star
    # intervene on received treatment while keeping set covariates fixed
    data.trtA$trt <- 1 # generate dataset where everyone receives treatment A
    data.trtC$trt <- 0 # generate dataset where all observations receive C
    # predict individual counterfactual event probs, conditional on treatment/covariates
    hat.mu.A.i <- predict(outcome.model, type="response", newdata=data.trtA)
    hat.mu.C.i <- predict(outcome.model, type="response", newdata=data.trtC)
    hat.mu.A <- mean(hat.mu.A.i) # (marginal) mean probability prediction under treatment A
    hat.mu.C <- mean(hat.mu.C.i) # (marginal) mean probability prediction under treatment C
    # estimate marginal A vs. C log-odds ratio (mean difference in expected log-odds)  
    # by transforming from probability to linear predictor scale 
    hat.Delta.AC <- log(hat.mu.A/(1-hat.mu.A)) - log(hat.mu.C/(1-hat.mu.C))    
    # hat.Delta.AC <- qlogis(hat.mu.A) - qlogis(hat.mu.C) 
    return(hat.Delta.AC)
  }  
  # non-parametric bootstrap
  boot.object <- boot::boot(data=data.AC, statistic=gcomp.ml, R=resamples)
  # bootstrap mean of marginal A vs. C treatment effect estimate
  hat.Delta.AC <- mean(boot.object$t)
  # bootstrap variance of A vs. C treatment effect estimate   
  hat.var.Delta.AC <- var(boot.object$t)
  # marginal log-odds ratio for B vs. C from reported event counts
  hat.Delta.BC <- with(data.BC,log(y.B.sum*(N.C-y.C.sum)/(y.C.sum*(N.B-y.B.sum))))
  # variance of B vs. C using delta method
  hat.var.Delta.BC <- with(data.BC,1/y.C.sum+1/(N.C-y.C.sum)+1/y.B.sum+1/(N.B-y.B.sum))
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # treatment effect for A vs. B
  # variance for A vs. B
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC 
  list(hat.Delta.AB, hat.var.Delta.AB)
}  

### Bayesian G-computation
gcomp.bayes.wrapper <- function(data.AC, data.BC, n.chains, warmup, iters,
                                N_star) {
  # Inputs: data.AC - AC individual patient-level data; data.BC - BC aggregate-level data;
  # n.chains, burnin, iters - MCMC info
  # N_star - size of simulated BC pseudo-population (high for small Monte Carlo error)
  # matrix of pairwise correlations between IPD covariates  
  rho <- cor(data.AC[,c("X1","X2","X3","X4")]) 
  #  covariate simulation for comparator trial using copula package
  cop <- normalCopula(param=c(rho[1,2],rho[1,3],rho[1,4],rho[2,3],rho[2,4],rho[3,4]), 
                      dim=4, dispstr="un") # AC IPD pairwise correlations
  # sample covariates from approximate joint distribution using copula
  mvd <- mvdc(copula=cop, margins=c("norm", "norm", "norm", "norm"), # Gaussian marginals
              # BC covariate means and standard deviations
              paramMargins=list(list(mean=data.BC$mean.X1, sd=data.BC$sd.X1),
                                list(mean=data.BC$mean.X2, sd=data.BC$sd.X2),       
                                list(mean=data.BC$mean.X3, sd=data.BC$sd.X3),
                                list(mean=data.BC$mean.X4, sd=data.BC$sd.X4)))
  # data frame of simulated covariates
  x_star <- as.data.frame(rMvdc(N_star, mvd))
  colnames(x_star) <- c("X1", "X2", "X3", "X4")  
  # outcome logistic regression fitted to original data using MCMC  
  outcome.model <- stan_glm(y~X3+X4+trt*X1+trt*X2, data=data.AC, 
                            family=binomial, algorithm="sampling",
                            iter=iters, warmup=warmup, chains=n.chains) # run Stan model
  # counterfactual datasets
  data.trtA <- data.trtC <- x_star
  # intervene on received treatment while keeping set covariates fixed
  data.trtA$trt <- 1 # generate dataset where everyone receives treatment A
  data.trtC$trt <- 0 # generate dataset where all observations receive C  
  # draw counterfactual binary responses from posterior predictive distribution
  # matrix of posterior predictive draws under treatment A
  y.star.A <- posterior_predict(outcome.model, newdata=data.trtA) 
  # matrix of posterior predictive draws under treatment C
  y.star.C <- posterior_predict(outcome.model, newdata=data.trtC)
  # compute marginal log-odds ratio for A vs. C for each MCMC sample
  # by transforming from probability to linear predictor scale  
  hat.delta.AC <- qlogis(rowMeans(y.star.A)) - qlogis(rowMeans(y.star.C)) 
  hat.Delta.AC <- mean(hat.delta.AC) # average over samples
  hat.var.Delta.AC <- var(hat.delta.AC) # sample variance
  # B vs. C from reported aggregate event counts, e.g. in contingency table
  hat.Delta.BC <- with(data.BC, log(y.B.sum*(N.C-y.C.sum)/(y.C.sum*(N.B-y.B.sum))))
  # B vs. C variance using the delta method 
  hat.var.Delta.BC <- with(data.BC, 1/y.C.sum+1/(N.C-y.C.sum)+1/y.B.sum+1/(N.B-y.B.sum))
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # treatment effect for A vs. B
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC # A vs. B variance
  list(hat.Delta.AB, hat.var.Delta.AB)
}

# set up cluster for parallel computing
num.cores <- detectCores()-1
cluster <- makeCluster(num.cores, type="SOCK", outfile="")
registerDoSNOW(cluster)
# progress bar
pb <- txtProgressBar(max=replicates, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

# combine lists in parallelisation
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# run population adjustment methods for all replicates/scenarios in parallel
for(i in 1:scenarios) {
  IPD.AC <- IPD.AC.all[[i]]
  ALD.BC <- ALD.BC.all[[i]]
  file.id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i])  
  ### Matching-adjusted indirect comparison (MAIC)
  maic.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                          .init=list(list(), list()), .options.snow=opts,
                          .packages=c("boot")) %dopar% {
                            results <- maic.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                    resamples)
                            return(results)
                          }
  close(pb)
  means <- unlist(maic.results[[1]])
  variances <- unlist(maic.results[[2]])
  save(means, file=paste0("Results/MAIC/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/MAIC/variances_", file.id, ".RData"))
  ## Simulated treatment comparison (STC); conventional outcome regression approach
  stc.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                         .init=list(list(),list()), .options.snow=opts) %dopar% {
                           results <- stc.wrapper(IPD.AC[[j]], ALD.BC[[j]])
                           return(results)
                         }
  close(pb)
  means <- unlist(stc.results[[1]])
  variances <- unlist(stc.results[[2]])
  save(means, file=paste0("Results/STC/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/STC/variances_", file.id, ".RData"))
  ### G-computation with maximum-likelihood estimation and bootstrapping
  gcomp.ml.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                              .init=list(list(),list()), .options.snow=opts,
                              .packages=c("boot", "copula")) %dopar% {
                                results <- gcomp.ml.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                            resamples, N_star)
                                return(results)
                              }
  means <- unlist(gcomp.ml.results[[1]])
  variances <- unlist(gcomp.ml.results[[2]])
  save(means, file=paste0("Results/GcompML/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/GcompML/variances_", file.id, ".RData"))
  ### Bayesian G-computation
  gcomp.bayes.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                                 .init=list(list(),list()), .options.snow=opts,
                                 .packages=c("copula", "rstanarm")) %dopar% {
                                  results <- gcomp.bayes.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                                 n.chains, warmup, iters, N_star)
                                 return(results)
                                 }
  means <- unlist(gcomp.bayes.results[[1]])
  variances <- unlist(gcomp.bayes.results[[2]])
  save(means, file=paste0("Results/GcompBayes/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/GcompBayes/variances_", file.id, ".RData"))
}

stopCluster(cluster)


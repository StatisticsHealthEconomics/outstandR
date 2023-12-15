# perform the indirect comparison methods on the simulated data

# load sim. study parameter combinations
load(here::here("data", "binary_settings_input_params.RData")) 

library(parallel)
library(doSNOW)
library(rstanarm)  # fit outcome regression and draw predicted outcomes from posterior in Bayes G-comp
library(boot)      # non-parametric bootstrap in ML G-computation and MAIC
library(copula)    # simulate BC trial covariates from multivariate Gaussian copula

set.seed(444)
n_scenarios <- nrow(pc) # number of simulation scenarios

# settings for the methods
resamples <- 1000 # Number of bootstrap resamples for MAIC and maximum-likelihood G-computation
N_star <- 1000    # size of simulated BC pseudo-population; high for small MC error

# MCMC info for Bayesian G-computation
n.chains <- 2     # number of Markov chains for MCMC
warmup <- 2000    # number of warmup aka discarded burn-in iterations per chain
iters <- 4000     # number of total iterations per chain (including warmup)

# simulated patient-level (AC) and aggregate-level (BC) datasets for all scenarios
IPD.AC.all <- vector(mode="list", n_scenarios)
ALD.BC.all <- vector(mode="list", n_scenarios)

# load data
for (i in seq_len(n_scenarios)) {
  file_id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i]) 
  load(paste0("data/IPD_AC_", file_id, ".RData"))  # AC patient-level data
  load(paste0("data/ALD_BC_", file_id, ".RData"))  # BC aggregate-level data
  IPD.AC.all[[i]] <- IPD.AC
  ALD.BC.all[[i]] <- ALD.BC
}

replicates <- N_sim  # number of Monte Carlo replicates per scenario
repl_seq <- 1:replicates

# set up cluster for parallel computing
num.cores <- detectCores() - 1
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


# run population adjustment methods
# for all replicates/scenarios in parallel
for(i in seq_len(n_scenarios)) {
  IPD.AC <- IPD.AC.all[[i]]
  ALD.BC <- ALD.BC.all[[i]]
  file_id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i])
  
  ### Matching-adjusted indirect comparison (MAIC)
  maic.results <- foreach(j=repl_seq, .combine='comb', .multicombine=TRUE,
                          .init=list(list(), list()), .options.snow=opts,
                          .packages=c("boot")) %dopar% {
                            results <- maic.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                    resamples)
                            return(results)
                          }
  close(pb)
  means <- unlist(maic.results[[1]])
  variances <- unlist(maic.results[[2]])
  
  save(means, file=paste0("Results/MAIC/means_", file_id, ".RData"))
  save(variances, file=paste0("Results/MAIC/variances_", file_id, ".RData"))
  
  ## Simulated treatment comparison (STC); conventional outcome regression approach
  stc.results <- foreach(j=repl_seq, .combine='comb', .multicombine=TRUE,
                         .init=list(list(),list()), .options.snow=opts) %dopar% {
                           results <- stc.wrapper(IPD.AC[[j]], ALD.BC[[j]])
                           return(results)
                         }
  close(pb)
  means <- unlist(stc.results[[1]])
  variances <- unlist(stc.results[[2]])
  
  save(means, file=paste0("Results/STC/means_", file_id, ".RData"))
  save(variances, file=paste0("Results/STC/variances_", file_id, ".RData"))
  
  ### G-computation with maximum-likelihood estimation and bootstrapping
  gcomp.ml.results <- foreach(j=repl_seq, .combine='comb', .multicombine=TRUE,
                              .init=list(list(),list()), .options.snow=opts,
                              .packages=c("boot", "copula")) %dopar% {
                                results <- gcomp.ml.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                            resamples, N_star)
                                return(results)
                              }
  means <- unlist(gcomp.ml.results[[1]])
  variances <- unlist(gcomp.ml.results[[2]])
  
  save(means, file=paste0("Results/GcompML/means_", file_id, ".RData"))
  save(variances, file=paste0("Results/GcompML/variances_", file_id, ".RData"))
  
  ### Bayesian G-computation
  gcomp.bayes.results <- foreach(j=repl_seq, .combine='comb', .multicombine=TRUE,
                                 .init=list(list(),list()), .options.snow=opts,
                                 .packages=c("copula", "rstanarm")) %dopar% {
                                  results <- gcomp.bayes.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                                 n.chains, warmup, iters, N_star)
                                 return(results)
                                 }
  means <- unlist(gcomp.bayes.results[[1]])
  variances <- unlist(gcomp.bayes.results[[2]])
  
  save(means, file=paste0("Results/GcompBayes/means_", file_id, ".RData"))
  save(variances, file=paste0("Results/GcompBayes/variances_", file_id, ".RData"))
}

stopCluster(cluster)


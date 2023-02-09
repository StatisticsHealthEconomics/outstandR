# specify the simulation set-up and generates the simulation study data


library(dplyr)
library(MASS) # sample/simulate covariates from a multivariate normal

set.seed(555) # set random seed for reproducibility

# Define simulation study parameters

N_sim <- 2000     # number of Monte Carlo replicates
allocation <- 2/3 # active treatment vs. placebo allocation ratio (2:1)

N_AC <- c(200,400,600) # number of subjects in the AC trial
N_BC <- 600            # number of subjects in the BC trial 
b_trt <- log(0.17)     # conditional effect of active treatment vs. common comparator
b_X <- -log(0.5)       # conditional effect of each prognostic variable
b_EM <- -log(0.67)     # conditional interaction effect of each effect modifier
event_rate <- 0.35
meanX_AC <- c(0.45,0.3,0.15) # mean of each normally-distributed covariate in AC trial
meanX_BC <- 0.6        # mean of each normally-distributed covariate in BC
sdX <- 0.4             # standard deviation of each covariate (same for AC and BC)
corX <- 0.2            # covariate correlation coefficient  

# parameter combinations for each scenario
param.combinations <- expand.grid(N_AC=N_AC, meanX_AC=meanX_AC)
pc <- param.combinations

scenarios <- nrow(pc) # number of simulation scenarios
save(pc, N_sim, allocation, file="data/binary_settings.RData")

# compute intercept based on baseline event rate (for patients under C without covariates)
optim.function <- function(param, y_prob) {
  b_0 <- param
  fit <- sum(((1 / (1 + exp(-b_0))) - y_prob)^2)
  return(fit)
}

b_0 <- optim(par=0, fn=optim.function,
             y_prob=event_rate,
             method="Brent",
             lower=-2, upper=2)$par

for (i in seq_len(scenarios)) {
  print(i)
  # simulate IPD covariates and outcome for A vs. C trial (S=1)
  IPD.AC <- replicate(n=N_sim, expr=gen_data(pc$N_AC[i], b_trt, b_X, b_EM, 
                                             b_0, pc$meanX_AC[i], sdX, event_rate, 
                                             corX, allocation),
                      simplify=FALSE)
  # simulate IPD covariates and outcome for B vs. C trial (S=2)
  IPD.BC <- replicate(n=N_sim, expr=gen_data(N_BC, b_trt, b_X, b_EM, 
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
                      
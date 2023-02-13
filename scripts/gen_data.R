# specify the simulation set-up and generates the simulation study data


library(dplyr)
library(MASS)  # sample/simulate covariates from multivariate normal

set.seed(555)

# Define simulation study parameters

n_sim <- 2000      # number of Monte Carlo replicates
allocation <- 2/3  # active treatment vs. placebo allocation ratio (2:1)

N_AC <- c(200, 400, 600) # number of subjects in the AC trial
N_BC <- 600            # number of subjects in the BC trial 
b_trt <- log(0.17)     # conditional effect of active treatment vs. common comparator
b_X <- -log(0.5)       # conditional effect of each prognostic variable
b_EM <- -log(0.67)     # conditional interaction effect of each effect modifier
event_rate <- 0.35
meanX_AC <- c(0.45, 0.3, 0.15) # mean of each normally-distributed covariate in AC trial
meanX_BC <- 0.6        # mean of each normally-distributed covariate in BC
sdX <- 0.4             # standard deviation of each covariate (same for AC and BC)
corX <- 0.2            # covariate correlation coefficient  

# parameter combinations for each scenario
pc_AC <- expand.grid(N = N_AC, meanX = meanX_AC)

scenarios <- nrow(pc_AC)    # number of simulation scenarios
save(pc_AC, n_sim, allocation,
     file = here::here("data", "binary_settings.RData"))

b_0 <- optim(par = 0,
             fn = function(param, y_prob)
               sum(((1 / (1 + exp(-param))) - y_prob)^2),
             y_prob = event_rate,
             method = "Brent",
             lower = -2, upper = 2)$par

gen_data_args <-
  tibble::lst(b_trt, b_X,
              b_EM, b_0,
              sdX, event_rate,
              corX, allocation)

# simulate IPD covariates and outcome for B vs. C trial (S=2)
IPD.BC <-
  replicate(n = n_sim,
            expr = do.call(gen_data,
                           c(N = N_BC, meanX = meanX_BC, gen_data_args)), 
            simplify = FALSE)

# summarize BC IPD as ALD
ALD.BC <- lapply(1:n_sim, function(j) {
  as.data.frame(cbind(
    # aggregate data for the BC trial 
    summarise(
      IPD.BC[[j]],
      mean.X1 = mean(X1),
      mean.X2 = mean(X2),
      mean.X3 = mean(X3),
      mean.X4 = mean(X4),
      sd.X1 = sd(X1),
      sd.X2 = sd(X2),
      sd.X3 = sd(X3),
      sd.X4 = sd(X4)),
    # summarize outcomes for the BC trial (treatment B)
    filter(IPD.BC[[j]], trt == 1) %>%
      summarise(y.B.sum = sum(y),
                y.B.bar = mean(y),
                N.B = n()),
    # summarize outcomes for the BC trial (treatment C)
    filter(IPD.BC[[j]], trt == 0) %>%
      summarise(y.C.sum = sum(y),
                y.C.bar = mean(y),
                N.C = n())))   
})

save(IPD.BC, file = "Data/IPD_BC.RData")
save(ALD.BC, file = "Data/ALD_BC.RData")  


# simulate IPD covariates and outcome for A vs. C trial (S=1)
for (i in seq_len(scenarios)) {
  
  params <- pc_AC[i, ]
  file_id <- glue::glue("N_AC{params$N}meanX_AC{params$meanX}")
  print(file_id)
  
  replicate(n = n_sim,
            expr =
              do.call(gen_data,
                      c(N = params$N, meanX = params$meanX, gen_data_args)),
            simplify = FALSE) |> 
    save(file = glue::glue("Data/IPD_AC_{file_id}.RData"))
}                       


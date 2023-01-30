## This file processes the results of the simulation study and computes
## and graphs the relevant performance metrics.

rm(list=ls())

# setwd("C:/Users/Antonio/Desktop/Gcomp_indirect_comparisons_simstudy") 
load(file="binary_settings.RData")
source("functions.R") # load functions to compute performance measures

# packages required for plotting simulation study results
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(gridExtra)) {install.packages("gridExtra"); library(gridExtra)}

scenario.settings <- pc # parameter combinations
scenarios <- nrow(pc) # number of scenarios

replicates <- N_sim # number of data replicates per scenario
Delta.AB <- 0 # true value of marginal A vs. B treatment effect is zero

# data frame to store performance metrics (4 methods, each row repeated 4 times)
simulation.metrics <- scenario.settings[rep(seq_len(nrow(scenario.settings)), 
                                            each=4), ]
metrics.names <- c("Method",
                   "Bias", "Bias.MCSE", "LCI", "LCI.MCSE", "UCI", "UCI.MCSE",
                   "VR", "VR.MCSE", "Cov", "Cov.MCSE", "ESE", "ESE.MCSE",
                   "MSE", "MSE.MCSE", "MAE", "MAE.MCSE", "number.NAs",
                   "negative.vars")
simulation.metrics[metrics.names] <- NA

# data frame to store all A vs. B marginal treatment effect point estimates
ate.table <- scenario.settings[rep(seq_len(nrow(scenario.settings)),
                                   each=4*replicates),]
ate.table["Method"] <- NA
ate.table["ATE"] <- NA

# function that computes performance metrics for a given method
process.metrics <- function(means, variances, truth) {
  # remove NAs (an issue for MAIC in Scenario 7, separation issues)
  NAs <- is.na(means)
  means <- means[!NAs]  
  variances <- variances[!NAs]
  number.NAs <- sum(NAs) # number that do not converge
  replicates <- length(means)
  bias.metric <- bias(means, truth)
  bias.metric.mcse <- bias.mcse(means)
  # # bias is equal to estimated marginal treatment effect (truth is zero)
  # ate <- mean(means)
  # ate.mcse <- mcse.estimate(means)
  mae.metric <- mae(means, truth)
  abs.err <- abs(means - truth)
  mae.mcse <- mcse.estimate(abs.err)
  mse.metric <- mse(means, truth) 
  mse.metric.mcse <- mse.mcse(means, truth) 
  # construct Wald-type interval estimates using normal distribution
  lci <- means + qnorm(0.025)*sqrt(variances)
  uci <- means + qnorm(0.975)*sqrt(variances)
  lci.mean <- mean(lci)
  lci.mcse <- mcse.estimate(lci)
  uci.mean <- mean(uci)
  uci.mcse <- mcse.estimate(uci)
  cov <- coverage(lci, uci, truth)
  cov.mcse <- coverage.mcse(cov, replicates)
  empse.metric <- empse(means)
  empse.metric.mcse <- empse.mcse(empse.metric, replicates)
  vr <- var.ratio(means, sqrt(variances))
  vr.mcse <- var.ratio.mcse(avg.se=mean(sqrt(variances)), 
                            emp.se=empse.metric,
                            var.avg.se=mcse.estimate(sqrt(variances))^2,
                            var.emp.se=empse.metric.mcse^2)
  list(bias.metric, bias.metric.mcse, lci.mean, lci.mcse, uci.mean, uci.mcse,
       vr, vr.mcse, cov, cov.mcse, empse.metric, empse.metric.mcse,
       mse.metric, mse.metric.mcse, mae.metric, mae.mcse, number.NAs)
} 

j <- 1 # row counter for simulation metrics
k <- 1 # row counter for ATEs

for (i in 1:scenarios) {
  file.id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i])   
  ### Matching-adjusted indirect comparison (MAIC)
  load(paste0("Results/MAIC/means_", file.id, ".RData"))
  load(paste0("Results/MAIC/variances_", file.id, ".RData"))
  simulation.metrics[j,3] <- "MAIC"
  maic.metrics <- process.metrics(means, variances, Delta.AB) 
  simulation.metrics[j,4:20] <- unlist(maic.metrics)
  ate.table[k:(k+replicates-1),3] <- "MAIC"
  ate.table[k:(k+replicates-1),4] <- means
  j <- j+1
  k <- k+replicates
  ### Simulated treatment comparison (STC); conventional outcome regression approach
  load(paste0("Results/STC/means_", file.id, ".RData"))
  load(paste0("Results/STC/variances_", file.id, ".RData"))  
  simulation.metrics[j,3] <- "STC"
  stc.metrics <- process.metrics(means, variances, Delta.AB)
  simulation.metrics[j,4:20] <- unlist(stc.metrics)
  ate.table[k:(k+replicates-1),3] <- "STC"
  ate.table[k:(k+replicates-1),4] <- means
  j <- j+1
  k <- k+replicates  
  ### G-computation with maximum-likelihood estimation and bootstrapping
  load(paste0("Results/GcompML/means_", file.id, ".RData"))
  load(paste0("Results/GcompML/variances_", file.id, ".RData"))    
  simulation.metrics[j,3] <- "G-comp (ML)"
  gcomp.ml.metrics <- process.metrics(means, variances, Delta.AB)
  simulation.metrics[j,4:20] <- unlist(gcomp.ml.metrics)
  ate.table[k:(k+replicates-1),3] <- "G-comp (ML)"
  ate.table[k:(k+replicates-1),4] <- means
  j <- j+1
  k <- k+replicates    
  ### Bayesian G-computation
  load(paste0("Results/GcompBayes/means_", file.id, ".RData"))
  load(paste0("Results/GcompBayes/variances_", file.id, ".RData"))
  simulation.metrics[j,3] <- "G-comp (Bayes)"
  gcomp.bayes.metrics <- process.metrics(means, variances, Delta.AB)
  simulation.metrics[j,4:20] <- unlist(gcomp.bayes.metrics)
  ate.table[k:(k+replicates-1),3] <- "G-comp (Bayes)"
  ate.table[k:(k+replicates-1),4] <- means
  j <- j+1
  k <- k+replicates      
}

# Save simulation study performance metrics
write.csv(simulation.metrics, "Analysis/scenarios.csv", row.names = FALSE)

## Plot results for a specific scenario
i <- 1 # change from 1-9 for different scenarios
scenario.ates <- subset(ate.table, meanX_AC==pc$meanX_AC[i]&N_AC==pc$N_AC[i])
scenario.metrics <- subset(simulation.metrics, 
                           meanX_AC==pc$meanX_AC[i]&N_AC==pc$N_AC[i])
display.table <- cbind(Method=scenario.metrics$Method, 
                       ATE=paste0(format(round(scenario.metrics$Bias,digits=3),nsmall=3)," (",
                                  format(round(scenario.metrics$Bias.MCSE,digits=3),nsmall=3),")"),
                       LCI=paste0(format(round(scenario.metrics$LCI,digits=3),nsmall=3)," (",
                                  format(round(scenario.metrics$LCI.MCSE,digits=3),nsmall=3),")"),
                       UCI=paste0(format(round(scenario.metrics$UCI,digits=3),nsmall=3)," (",
                                  format(round(scenario.metrics$UCI.MCSE,digits=3),nsmall=3),")"),
                       VR=paste0(format(round(scenario.metrics$VR,digits=3),nsmall=3)," (",
                                 format(round(scenario.metrics$VR.MCSE,digits=3),nsmall=3),")"),
                       Cov=paste0(format(round(scenario.metrics$Cov,digits=3),nsmall=3)," (",
                                  format(round(scenario.metrics$Cov.MCSE,digits=3),nsmall=3),")"),
                       ESE=paste0(format(round(scenario.metrics$ESE,digits=3),nsmall=3)," (",
                                  format(round(scenario.metrics$ESE.MCSE,digits=3),nsmall=3),")"),
                       MSE=paste0(format(round(scenario.metrics$MSE,digits=3),nsmall=3)," (",
                                  format(round(scenario.metrics$MSE.MCSE,digits=3),nsmall=3),")")) 
p1 <- ggplot(scenario.ates, aes(x=Method, y=ATE, fill=Method)) + 
  geom_boxplot(alpha=0.7) + geom_hline(yintercept=0, linetype="dashed", color = "red") +
  scale_x_discrete(limits=c("MAIC", "STC", "G-comp (ML)", 
                            "G-comp (Bayes)")) +
  scale_fill_brewer(palette="Dark2") + theme_classic() + 
  theme(legend.position="none", axis.text.x = element_text(color = "grey20", size = 14, 
                                                           face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, 
                                   face = "plain"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color= "grey20", size=14, face="plain")) +
  ## Additional options for outlier in scenario 7
  # coord_cartesian(ylim = c(-4, 4)) + geom_point(aes(x="MAIC", y=-4), colour="red") +
  # geom_text(x="MAIC", y=-4, label="-10.3",hjust = -0.3,color="red") +
  ylab("Marginal treatment effect")
p2 <- tableGrob(display.table, theme=ttheme_minimal())
grid.arrange(p1,p2,ncol=1) 
## ----include = FALSE------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup, warning=FALSE, message=FALSE----------------------------------
library(boot)      # non-parametric bootstrap in MAIC and ML G-computation
library(copula)    # simulating BC covariates from Gaussian copula
library(rstanarm)  # fit outcome regression, draw outcomes in Bayesian G-computation
library(mimR)


## ----load-data------------------------------------------------------------
set.seed(555)

AC.IPD <- read.csv(here::here("data", "AC_IPD.csv"))  # AC patient-level data
BC.ALD <- read.csv(here::here("data", "BC_ALD.csv"))  # BC aggregate-level data


## -------------------------------------------------------------------------
head(AC.IPD)


## -------------------------------------------------------------------------
BC.ALD


## ----hat_Delta_stats_maic-------------------------------------------------
hat_Delta_stats_maic <- hat_Delta_stats(AC.IPD, BC.ALD, strategy = strategy_maic())


## ----hat_Delta_stats_maic-print-------------------------------------------
hat_Delta_stats_maic


## ----hat_Delta_stats_stc--------------------------------------------------
hat_Delta_stats_stc <- hat_Delta_stats(AC.IPD, BC.ALD, strategy = strategy_stc())
hat_Delta_stats_stc


## ----hat_Delta_stats_gcomp_ml---------------------------------------------
# hat_Delta_stats_gcomp_ml <- hat_Delta_stats(AC.IPD, BC.ALD, strategy = strategy_gcomp_ml())


## ----hat_Delta_stats_gcomp_stan-------------------------------------------
hat_Delta_stats_gcomp_stan <- hat_Delta_stats(AC.IPD, BC.ALD, strategy = strategy_gcomp_stan())
hat_Delta_stats_gcomp_stan


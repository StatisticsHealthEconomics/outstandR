## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(boot)      # non-parametric bootstrap in MAIC and ML G-computation
library(copula)    # simulating BC covariates from Gaussian copula
library(rstanarm)  # fit outcome regression, draw outcomes in Bayesian G-computation
library(outstandR)

## ----load-data----------------------------------------------------------------
set.seed(555)

AC.IPD <- read.csv(here::here("raw-data", "AC_IPD.csv"))  # AC patient-level data
BC.ALD <- read.csv(here::here("raw-data", "BC_ALD.csv"))  # BC aggregate-level data

## -----------------------------------------------------------------------------
head(AC.IPD)

## -----------------------------------------------------------------------------
BC.ALD

## -----------------------------------------------------------------------------
lin_form <- as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")

## ----outstandR_maic------------------------------------------------------------
outstandR_maic <- outstandR(AC.IPD, BC.ALD, strategy = strategy_maic(formula = lin_form))

## ----outstandR_maic-print------------------------------------------------------
outstandR_maic

## ----outstandR_stc-------------------------------------------------------------
outstandR_stc <- outstandR(AC.IPD, BC.ALD, strategy = strategy_stc(formula = lin_form))
outstandR_stc

## ----outstandR_gcomp_ml--------------------------------------------------------
# outstandR_gcomp_ml <- outstandR(AC.IPD, BC.ALD, strategy = strategy_gcomp_ml(formula = lin_form))

## ----outstandR_gcomp_stan------------------------------------------------------
outstandR_gcomp_stan <- outstandR(AC.IPD, BC.ALD, strategy = strategy_gcomp_stan(formula = lin_form))
outstandR_gcomp_stan


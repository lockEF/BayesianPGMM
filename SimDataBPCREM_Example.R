# Load packages
library(MASS)
library(matrixStats)
library(rjags)
library(coda)

# Load BPCREM function
source("BayesPCREM.R")

# Load simulated data
load(file = "SimDataBPCREM.Rdata")

# Estimate the bivariate model (can take about 20-30 minutes)
Results <- BayesPCREM(data = SimDataBPCREM,
                      y_var = "y", time_var = "time", 
                      ind_id_var = "id", cross_id_var = "teacherid",
                      TracePlots=TRUE)

# Mean potential scale reduction factor to assess convergence
Results$mean.psrf

# DIC
Results$DIC

# Parameter estimates and their 95% credible intervals
Results$parameter.estimates

# Run time for each component
Results$run.times

# Traceplots
traceplot(Results$mcmcList) # from coda package
plot(Results$mcmcList) # base R function
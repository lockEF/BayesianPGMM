# Required libraries
library(MASS)          
library(matrixStats)   
library(rjags)
library(coda)

# Load BayesBPLMEM function
source("BayesBPLMEM.R")

# Load simulated data
load(file = "SimDataBPLMEM.rda")

# Estimate the bivariate model (can take about 45 minutes)
Results <- BayesBPLMEM(X = SimDataBPLMEM[[1]], 
                       Y = SimDataBPLMEM[[2]],
                       iters_adapt = 5000, 
                       iters_BurnIn = 10000, 
                       iters_sampling = 10000, 
                       Thin=15)

# Mean potential scale reduction factor to assess convergence
Results$Gelman.msrf

# DIC
Results$DIC

# Fixed-effects of the two outcome variables and their respective 95% credible intervals
Results$y1.mean
Results$y2.mean
Results$y1.mean.CI
Results$y2.mean.CI

# Error variances, covariance, and correlation and their respective 95% credible intervals
Results$error.var
Results$error.cov
Results$error.corr
Results$error.var.CI
Results$error.cov.CI
Results$error.corr.CI

# Variances, covariances, and correlations of random-effects and their respective 95% credible intervals
# Suffixes of covariances and correlations refer to the entries of the covariance/correlation matrix 
Results$random.var
Results$random.cov
Results$random.corr
Results$random.var.CI
Results$random.cov.CI
Results$random.corr.CI






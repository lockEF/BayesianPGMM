BayesBPLMEM <- function(X, Y, iters_adapt=5000, iters_BurnIn=100000, iters_sampling=50000, Thin=15, SaveChains=FALSE){
  #X[i,j] gives timepoint for subject i at measurement j
  #Y[i,j] gives value for subject i at measurement j
  
  #load module to compute DIC
  load.module("dic")
  
  ###Bivariate piecewise model definition for JAGS 
  bivariate_pw <- "model{
  for (i in 1:Nind) {
    for (j in 1:Ntime) {
      y[i,j,1:2] ~ dmnorm(mu_y[i,j,1:2], tau_error[1:2,1:2])
      mu_y[i,j,1] <- b[i,1] + b[i,2]*x[i,j] + b[i,3]*(max(0, x[i,j]-b[i,4])) 
      mu_y[i,j,2] <- b[i,5] + b[i,6]*x[i,j] + b[i,7]*(max(0, x[i,j]-b[i,8]))
    }
  }
  
  ##Level 1 precision and error variance
  tau_error[1:2,1:2] <- inverse(sigma2_error[,])
  sigma2_error[1,1] <- sigma2_e1
  sigma2_error[2,2] <- sigma2_e2
  sigma2_error[2,1] <- rho*sqrt(sigma2_e1)*sqrt(sigma2_e2)
  sigma2_error[1,2] <- sigma2_error[2,1] 
  sigma2_e1 ~ dunif(0,bound_e1)
  sigma2_e2 ~ dunif(0,bound_e2)
  rho ~ dunif(-1,1)
  
  ##Distribution of random-effects and error variance
  for (i in 1:Nind){
    for(k in 1:8){
      b[i,k] <- mub[k] + b.rand[i,k]
      b.rand[i,k] <- c[k]*b.raw[i,k]
    }
    b.raw[i,1:8] ~ dmnorm(mub.zero[1:8], taub.raw[1:8,1:8]) 
  }
  
  ##Priors for fixed-effects
  mub[1] ~ dnorm(0, 0.0001)
  mub[2] ~ dnorm(0, 0.0001)
  mub[3] ~ dnorm(0, 0.0001)
  mub[4] ~ dnorm(mean.cp, prec.cp)T(minT,maxT)
  mub[5] ~ dnorm(0, 0.0001)
  mub[6] ~ dnorm(0, 0.0001)
  mub[7] ~ dnorm(0, 0.0001)
  mub[8] ~ dnorm(mean.cp, prec.cp)T(minT,maxT)
  
  ##Priors for scaling constants
  for(k in 1:8){
    c[k] ~ dunif(0, 100)
  }
  
  ##Prior for covariance matrix of random-effects
  #Inverse Wishart prior for the raw covariance matrix
  taub.raw[1:8,1:8] ~ dwish(OmegaB[1:8,1:8], 9)
  sigma2_b.raw[1:8,1:8] <- inverse(taub.raw[,])
  
  ##Define elements of correlation and covariance matrices to recover
  for(k in 1:8){
    for(k.prime in 1:8){
      rho_b[k,k.prime] <- sigma2_b.raw[k,k.prime]/sqrt(sigma2_b.raw[k,k]*sigma2_b.raw[k.prime,k.prime])
      cov_b_mat[k,k.prime] <- c[k]*c[k.prime]*sigma2_b.raw[k,k.prime]  
    }
  }
  
  #Variances:
  for(k in 1:8){
    var_b[k] <- cov_b_mat[k,k]
  }
  
  #Covariances
  cov_b[1] <- cov_b_mat[2,1]
  cov_b[2] <- cov_b_mat[3,1]
  cov_b[3] <- cov_b_mat[3,2]
  cov_b[4] <- cov_b_mat[4,1]
  cov_b[5] <- cov_b_mat[4,2]
  cov_b[6] <- cov_b_mat[4,3]
  cov_b[7] <- cov_b_mat[5,1]
  cov_b[8] <- cov_b_mat[5,2]
  cov_b[9] <- cov_b_mat[5,3]
  cov_b[10] <- cov_b_mat[5,4]
  cov_b[11] <- cov_b_mat[6,1]
  cov_b[12] <- cov_b_mat[6,2]
  cov_b[13] <- cov_b_mat[6,3]
  cov_b[14] <- cov_b_mat[6,4]
  cov_b[15] <- cov_b_mat[6,5]
  cov_b[16] <- cov_b_mat[7,1]
  cov_b[17] <- cov_b_mat[7,2]
  cov_b[18] <- cov_b_mat[7,3]
  cov_b[19] <- cov_b_mat[7,4]
  cov_b[20] <- cov_b_mat[7,5]
  cov_b[21] <- cov_b_mat[7,6]
  cov_b[22] <- cov_b_mat[8,1]
  cov_b[23] <- cov_b_mat[8,2]
  cov_b[24] <- cov_b_mat[8,3]
  cov_b[25] <- cov_b_mat[8,4]
  cov_b[26] <- cov_b_mat[8,5]
  cov_b[27] <- cov_b_mat[8,6]
  cov_b[28] <- cov_b_mat[8,7]
  
  #Correlations
  cor_b[1] <- rho_b[2,1]
  cor_b[2] <- rho_b[3,1]
  cor_b[3] <- rho_b[3,2]
  cor_b[4] <- rho_b[4,1]
  cor_b[5] <- rho_b[4,2]
  cor_b[6] <- rho_b[4,3]
  cor_b[7] <- rho_b[5,1]
  cor_b[8] <- rho_b[5,2]
  cor_b[9] <- rho_b[5,3]
  cor_b[10] <- rho_b[5,4]
  cor_b[11] <- rho_b[6,1]
  cor_b[12] <- rho_b[6,2]
  cor_b[13] <- rho_b[6,3]
  cor_b[14] <- rho_b[6,4]
  cor_b[15] <- rho_b[6,5]
  cor_b[16] <- rho_b[7,1]
  cor_b[17] <- rho_b[7,2]
  cor_b[18] <- rho_b[7,3]
  cor_b[19] <- rho_b[7,4]
  cor_b[20] <- rho_b[7,5]
  cor_b[21] <- rho_b[7,6]
  cor_b[22] <- rho_b[8,1]
  cor_b[23] <- rho_b[8,2]
  cor_b[24] <- rho_b[8,3]
  cor_b[25] <- rho_b[8,4]
  cor_b[26] <- rho_b[8,5]
  cor_b[27] <- rho_b[8,6]
  cor_b[28] <- rho_b[8,7]
  
  #Define elements of error covariance (level 1) to recover
  sigma2_e[1] <- sigma2_error[1,1]
  sigma2_e[2] <- sigma2_error[2,2]
  sigma2_e[3] <- sigma2_error[2,1]
  
  #Collect important parameters to assess convergence
  ForConv[1:8] <- mub
  ForConv[9:11] <- sigma2_e
  ForConv[12:19] <- var_b
  ForConv[20:47] <- cov_b
  ForConv[48:75] <- cor_b
  ForConv[76] <- rho
  }"
  
  
  # Define values the models will require
  Nind <- dim(Y)[1] 
  Ntime <- dim(Y)[2]
  minCP <- unique(sort(X))[2]
  maxCP <- unique(sort(X, decreasing = TRUE))[2]
  minT <- min(X)  
  maxT <- max(X)
  mean.cp <- (maxT-minT)/2          # mean of changepoint
  prec.cp <- 1/((maxT-minT)/4)^2    # precision of changepoint
  bound_e1 <- min(colVars(Y[,,1], na.rm = TRUE))  # bound of error variance 1
  bound_e2 <- min(colVars(Y[,,2], na.rm = TRUE))  # bound of error variance 2
  mub.zero <- rep(0,8)  # Vector of zeros to center the distribution of raw random-effects
  OmegaB <- diag(8)     # Scale matrix for the Wishart distribution
  
  # Parameters for which I want to have the posterior distribution reported
  par_recov <- c("ForConv", "deviance", "pD")
  
  
  ### Estimate full bivariate model
  #Define JAGS model
  cat('Calibrating MCMC...\n')
  biv_spec <- textConnection(bivariate_pw)
  mod.jags <- jags.model(biv_spec,
                         data = list('y' = Y,
                                     'x' = X,
                                     'Nind' = Nind,
                                     'Ntime' = Ntime,
                                     'mean.cp' = mean.cp,
                                     'prec.cp' = prec.cp,
                                     'minT' = minT,
                                     'maxT' = maxT,
                                     'bound_e1' = bound_e1,
                                     'bound_e2' = bound_e2,
                                     'mub.zero' = mub.zero,
                                     'OmegaB' = OmegaB),
                         n.chains = 3,
                         n.adapt = iters_adapt)
  #Run model for burn-in period
  cat('Burn in of jags model...\n')
  update(mod.jags, iters_BurnIn)
  #Sample every Thin'th iteration for iters_sampling more iterations
  cat('Collecting samples...\n')
  fit_biv <- jags.samples(mod.jags, 
                          variable.names = par_recov, 
                          n.iter = iters_sampling, 
                          thin = Thin)
  
  # Elements to recover
  mcmcList <- as.mcmc.list(fit_biv$ForConv)
  gelman_msrf <- gelman.diag(mcmcList)
  mean_psrf <- mean(gelman_msrf$psrf[1:47])
  deviance <- mean(fit_biv$deviance)
  pD <- mean(fit_biv$pD)
  dic <- deviance + pD
  sum_mcmc <- summary(mcmcList)
  fixed_param1 <- sum_mcmc$statistics[1:4,1]
  fixed_param2 <- sum_mcmc$statistics[5:8,1]
  error_var <- sum_mcmc$statistics[9:10,1]
  error_cov <- sum_mcmc$statistics[11,1]
  random_var <- sum_mcmc$statistics[12:19,1]
  random_cov <- sum_mcmc$statistics[20:47,1]
  random_corr <- sum_mcmc$statistics[48:75,1]
  error_corr <- sum_mcmc$statistics[76,1]
  fixed_param1_CI <- cbind(sum_mcmc$quantiles[1:4,1], sum_mcmc$quantiles[1:4,5])
  fixed_param2_CI <- cbind(sum_mcmc$quantiles[5:8,1], sum_mcmc$quantiles[5:8,5])
  error_var_CI <- cbind(sum_mcmc$quantiles[9:10,1], sum_mcmc$quantiles[9:10,5])
  error_cov_CI <- cbind(sum_mcmc$quantiles[11,1], sum_mcmc$quantiles[11,5])
  random_var_CI <- cbind(sum_mcmc$quantiles[12:19,1], sum_mcmc$quantiles[12:19,5])
  random_cov_CI <- cbind(sum_mcmc$quantiles[20:47,1], sum_mcmc$quantiles[20:47,5])
  random_corr_CI <- cbind(sum_mcmc$quantiles[48:75,1], sum_mcmc$quantiles[48:75,5])
  error_corr_CI <- cbind(sum_mcmc$quantiles[76,1], sum_mcmc$quantiles[76,5])
  
  Results <- list("Gelman.msrf" = mean_psrf,
                  "DIC" = dic,
                  "y1.mean" = fixed_param1,
                  "y2.mean" = fixed_param2,
                  "error.var" = error_var,
                  "error.cov" = error_cov,
                  "error.corr" = error_corr,
                  "random.var" = random_var,
                  "random.cov" = random_cov,
                  "random.corr" = random_corr,
                  "y1.mean.CI" = fixed_param1_CI,
                  "y2.mean.CI" = fixed_param2_CI,
                  "error.var.CI" = error_var_CI,
                  "error.cov.CI" = error_cov_CI,
                  "error.corr.CI" = error_corr_CI,
                  "random.var.CI" = random_var_CI,
                  "random.cov.CI" = random_cov_CI,
                  "random.corr.CI" = random_corr_CI)
  if(SaveChains==TRUE){Results$Chains=fit_biv}
  class(Results) <- 'BPLMEM'
  return(Results)
}


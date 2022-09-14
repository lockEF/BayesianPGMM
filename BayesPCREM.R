##### Bayesian Piecewise Crossed Random Effects Model with an unknown Changepoint (BayesPCREM) ##### 

##### Packages needed to run function #####
# library(MASS)
# library(matrixStats)
# library(rjags)
# library(coda)
##### 

BayesPCREM <- function(data, ind_id_var, cross_id_var, time_var, y_var, 
                       iters_adapt=500, iters_BurnIn=20000, iters_sampling=30000, Thin=15, 
                       SaveChains=FALSE, TracePlots=FALSE){
  
  ###########################################################################################################################
  ##### Model Function Information #####
  ## INPUT: 
  ##    data = data in long form:
  ##           ind_id_var cross_id_var time_var y_var
  ##                    x            x        x     x
  ##                    x            x        x     x
  ##    ind_id_var = name of column that contains ids for individuals with repeated measures in a longitudinal dataset (e.g., students)
  ##    cross_id_var = name of column that contains ids for the crossed factor (e.g., teachers)
  ##    time_var = name of column that contains the time variable
  ##    y_var = name of column that contains the outcome variable
  ## ESTIMATION:
  ##    iters_adapt = 500, number of iterations for adaptation of jags model
  ##    iters_BurnIn = 20000, size of burn-in 
  ##    iters_sampling = 30000, number of iterations to monitor from posterior
  ##    Thin = 15, sample every Thin'th iteration .
  ## OUTPUT:
  ##    SaveChains = FALSE, logical expression to determine whether you want the 'Results' to store the MCMC chains
  ##    TracePlots = FALSE, logical expression to determine whether you want the 'Results' to store the mcmc.list object
  ##                 which can be used to create traceplots using 'plot()' or 'traceplot()' (`coda` package)
  ##    psrf = convergence criteria
  ##    DIC = for model comparison
  ##    parameter estimates (including CIs)
  ##    run times
  ###########################################################################################################################
  
  ##### load module to compute DIC #####
  load.module("dic")
  
  ##### Define values the models will require #####
  y <- data[,paste0(y_var)] # vector of all outcome observations
  t <- data[,paste0(time_var)] # vector of all timepoints observations
  ind_id <- as.numeric(factor(data[,paste0(ind_id_var)])) # vector of individual ids
  cross_id <- as.numeric(factor(data[,paste0(cross_id_var)])) # vector of crossed factor ids 
  
  ##### Error messages for input #####
  if(!is.data.frame(data)) stop('Expecting data.frame for data')
  if(sum(length(y)==length(t))<1) stop('Columns for y_var and t_var must have the same dimensions')
  if(sum(is.na(t))>0) stop('Columng for t_var cannot have NA values (but y_var can)')
  if(min(t)!=0) warning('Prior assumes first time point measured at t=0')
  
  ##### Define number of individuals, crossed factor clusters, and observations ##### 
  n <- length(unique(data[,paste0(ind_id_var)]))  # number of individuals (e.g., students) 
  K <- length(unique(data[,paste0(cross_id_var)])) # number of clusters in the crossed factor (e.g., teachers)
  N <- length(y) # number of total observations
  
  ##### Changepoint (CP) info ##### 
  minCP <- unique(sort(t))[2] # minimum timepoint for CP (second smallest timepoint value)
  maxCP <- unique(sort(t, decreasing = TRUE))[2] # maximum timepoint for CP (second largest timepoint value)
  minT <- min(unique(t)) # minimum timepoint for prior specification
  maxT <- max(unique(t)) # maximum timepoint for prior specification
  
  ##### Define mean, hyper mean, precision for random coefficients: ##### 
  mean = c(mean(subset(data, data[,time_var] == 0)[,y_var], na.rm = TRUE), 0, 0) # [mean at first timepoint (\beta_0), 0 (\beta_1), 0 (\beta_2)]
  PrecParamInt = 1/var(subset(data, data[,time_var] == 0)[,y_var], na.rm = TRUE) # precision for intercept
  PrecParam = (1/(sd(y, na.rm = TRUE)/(sd(t, na.rm = TRUE))))^2 # precision for betas
  prec = c(PrecParamInt, rep(PrecParam,2)) # precision vector for all beta params
  
  ##### Parameters for which I want to have the posterior distribution reported ##### 
  par_recov <- c("ForConv", "deviance", "pD")
  
  #########################
  ##### RUN TIME INFO #####
  #########################
  start_total <- Sys.time()
  
  ###################################
  ##### GENERATE INITIAL VALUES #####
  ###################################
  start_initialvals <- Sys.time()
  
  ##### Piecewise FEM #####
  fixed_pw_long <- "model{
  for (i in 1:N) {
      y[i] ~ dnorm(mu_y[i], tau_error)
      mu_y[i] <- b[1] + b[2]*t[i] + b[3]*(max(0, t[i]-b[4])) 
  }
  
  ## Level 1 precision and error variance ##
  sigma2_error <- 1/tau_error
  tau_error ~ dgamma(0.001,0.001)

  ## Priors for fixed-effects ##
  b[1] ~ dnorm(mean[1], prec[1])
  b[2] ~ dnorm(mean[2], prec[2])
  b[3] ~ dnorm(mean[3], prec[3])
  b[4] ~ dunif(minCP,maxCP)
  }"
  
  ##### Estimate Piecewise FEM -- to get initial values #####
  cat('Computing initial values...\n')
  fixed_spec <- textConnection(fixed_pw_long)
  dataList <- list('mean' = mean, 'prec' = prec, 't' = t, 'y' = y, 'N' = N, 'minCP' = minCP, 'maxCP' = maxCP)
  fixed.jags <- jags.model(fixed_spec, data = dataList, n.chains=3, n.adapt=500, quiet=TRUE)
  #update model and save samples
  update(fixed.jags, 2000) 
  Save.fixed <- jags.samples(fixed.jags, c('b','sigma2_error'), 2000, thin=4, quiet=TRUE) 
  
  ##### Compile Info #####
  Inits = list()
  Inits[[1]] = list(); Inits[[2]] = list(); Inits[[3]] = list()
  Inits[[1]]$beta0 = Inits[[2]]$beta0 = Inits[[3]]$beta0 = c()
  for(i in 1:3){ # 3 chains
    Inits[[i]]$beta0 = rowMeans(Save.fixed$b[,,i])
  }
  Inits[[1]]$taub = Inits[[2]]$taub = Inits[[3]]$taub = c()
  for(i in 1:3){
    Inits[[i]]$taub[1] <- prec[1]
    Inits[[i]]$taub[2:3] <- prec[2:3]
    Inits[[i]]$taub[4] <- (maxT-minT)/(4)
  }
  
  end_initialvals <- Sys.time()
  ###################################
  
  #########################
  ##### PIECEWISE CREM ####
  #########################
  # b0 + b1*t + b2*(max(0, t-gam))
  pw_crem <- "model{
  ##### Level-1 Model #####
  for (j in 1:N) {
    y[j] ~ dnorm(mu_y[j], tau_error)
    mu_y[j] <- betaik[j,1] + betaik[j,2]*t[j] + betaik[j,3]*(max(0, t[j]-betaik[j,4]))
    betaik[j,1] <- beta0[1] + b0[ind_id[j],1] + c0[cross_id[j],1]
    betaik[j,2] <- beta0[2] + b0[ind_id[j],2] + c0[cross_id[j],2]
    betaik[j,3] <- beta0[3] + b0[ind_id[j],3] + c0[cross_id[j],3]
    betaik[j,4] <- beta0[4] + b0[ind_id[j],4] + c0[cross_id[j],4]
  }

  ##### Loop over Individuals ##### 
  for(i in 1:n){
    for(p in 1:4){
      b0[i,p] ~ dnorm(0, taub[p])
    }
  }
  
  ##### Loop over Crossed Factor ##### 
  for(k in 1:K){
    for(p in 1:4){
      c0[k,p] ~ dnorm(0, tauc[p])
    }
  }
  
  ##### Priors #####
  ## Priors for Fixed Effects ##
  for(p in 1:3){
   beta0[p] ~ dnorm(0, 0.00001)
  }
  beta0[4] ~ dunif(minCP,maxCP)

  ## Priors for Residual Variance Components ##
  tau_error ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_error
  
  ## Priors for Random Effects ##
  for(p in 1:4){
  ## Priors for Individual Effects ##
    taub[p] ~ dgamma(0.001, 0.001)
    sigma2_b[p] <- 1/taub[p]
  
  ## Priors for Crossed Factor Effects ##
    tauc[p] ~ dgamma(0.001, 0.001)
    sigma2_c[p] <- 1/tauc[p]
  }

  ##### Collect important parameters to assess convergence #####
  ForConv[1:4] <- beta0
  ForConv[5] <- sigma2_error
  ForConv[6:9] <- sigma2_b
  ForConv[10:13] <- sigma2_c
  }"
  
  #########################################
  ##### ESTIMATE PIECEWISE CREM MODEL #####
  #########################################
  
  #############################
  ##### Define JAGS model ##### 
  start_mcmc <- Sys.time()
  cat('Calibrating MCMC...\n')
  ## Model definition for JAGS (with random effects): 
  pw_crem_spec<-textConnection(pw_crem)
  dataList <- list('t' = t, 'y' = y,
                   'N' = N, 'n' = n, 'K' = K,
                   'ind_id' = ind_id, 'cross_id' = cross_id,
                   'minCP' = minCP, 'maxCP' = maxCP)
  mod.jags <- jags.model(pw_crem_spec,
                         data = dataList,
                         inits = Inits, 
                         n.chains = 3,
                         n.adapt = iters_adapt, 
                         quiet = TRUE)
  end_mcmc <- Sys.time()
  #############################

  ###################
  ##### Burn-in ##### 
  start_burnin <- Sys.time()
  
  cat('Burn in of jags model...\n')
  update(mod.jags, iters_BurnIn)
  
  end_burnin <- Sys.time()
  ###################
  
  ###################
  ##### Samples ##### 
  start_samples <- Sys.time()
  
  cat('Collecting samples...\n')
  fit_pw_crem <- jags.samples(mod.jags,
                              variable.names = par_recov, # character vector giving names of variables to be monitored
                              n.iter = iters_sampling, # number of iterations to monitor
                              thin = Thin) # thinning interval for monitors
  
  end_samples <- Sys.time()
  ###################
  
  ## Run time total
  end_total <- Sys.time()
  #########################
  
  ###############################
  ##### ELEMENTS TO RECOVER #####
  param.names <- c("beta00", "beta10", "beta20", "betag0",
                   "ve", 
                   "var.b00", "var.b10", "var.b20", "var.bg0",
                   "var.c00", "var.c10", "var.c20", "var.cg0")
  ## Convergence
  mcmcList <- as.mcmc.list(fit_pw_crem$ForConv) # used for traceplots
  for(i in 1:3){colnames(mcmcList[[i]]) <- param.names}
  gelman_msrf <- gelman.diag(mcmcList) # individual parameter psrf, and multivariate psrf
  parameter_psrf <- data.frame(
    point_est = gelman_msrf$psrf[,1],
    upper_ci = gelman_msrf$psrf[,2]
  )
  multivariate_psrf <- gelman_msrf$mpsrf
  mean_psrf <- mean(parameter_psrf$point_est) # mean psrf across parameters
  ## Model Fit
  deviance <- mean(fit_pw_crem$deviance) # expected deviance
  pD <- mean(fit_pw_crem$pD) # effective number of parameters (expected deviance - fitted deviance)
  dic <- deviance + pD
  ## Parameter Estimates
  sum_mcmc <- summary(mcmcList) # parameter estimates
  parameter_estimates <- data.frame(
    Mean = sum_mcmc$statistics[,1],
    SD = sum_mcmc$statistics[,2],
    CI.Lower = sum_mcmc$quantiles[,1],
    CI.Upper = sum_mcmc$quantiles[,5]
  )
  ## Run Time
  run_time_initialvals <- end_initialvals - start_initialvals
  run_time_mcmc <- end_mcmc - start_mcmc
  run_time_burnin <- end_burnin - start_burnin
  run_time_samples <- end_samples - start_samples
  run_time_total <- end_total - start_total
  
  run_time <- data.frame(
    Operation = c("Obtaining Initial Values", "MCMC", "Burn-In", "Obtaining Samples", "Total"),
    Run.Time = c(format(run_time_initialvals), format(run_time_mcmc), format(run_time_burnin), format(run_time_samples), format(run_time_total))
  )
  ###############################
  
  ##########################
  ##### RESULTS OUTPUT #####
  ##########################
  Results <- list("parameter.psrf" = parameter_psrf,
                  "multivariate.psrf" = multivariate_psrf,
                  "mean.psrf" = mean_psrf,
                  "DIC" = dic,
                  "parameter.estimates" = parameter_estimates,
                  "run.times" = run_time)
  # optional results output
  if(SaveChains==TRUE){Results$Chains=fit_pw_crem}
  if(TracePlots==TRUE){Results$mcmcList=mcmcList} 
  class(Results) <- 'PCREM'
  return(Results)
}

BayesPGM = function(X, Y, max_cp=2, cp_prior='binomial', scale_prior='uniform', binom_prob=0.5, iters_BurnIn=20000, iters_sampling=30000, Thin=15, SaveChains=FALSE){
  ##X[i,j] gives timepoint for subject i at measurement j
  ##Y[i,j] gives value for subject i at measurement j
  
  ##Error messages for input
  if(!is.matrix(Y)) stop('Expecting matrix for input Y')
  if(!is.matrix(X)) stop('Expecting matrix for input X')
  if(sum(dim(X)==dim(Y))<2) stop('X and Y must have the same dimensions')
  if(sum(is.na(X)>0)) stop('X cannot have NA values (but Y can)')
  if(cp_prior!='binomial'&&cp_prior!='uniform') stop('Input for cp_prior must be \'binomial\' or \'uniform\'')
  if(scale_prior!='uniform'&&scale_prior!='hc') stop('Input for scale_prior must be \'uniform\' or \'hc\'')
  if(min(X)!=0) warning('Prior assumes first time point measured at x=0')
  
  #load module to compute DIC  
  load.module("dic")
  
  ##Define mean, hyper mean, precision for random coefficients:
  mean = c(mean(Y[,1],na.rm=TRUE),0,rep(0,max_cp))
  PrecParamInt = 1/var(Y[,1],na.rm=TRUE)
  PrecParam = (1/(max_cp*sd(Y,na.rm=TRUE)/(sd(X))))^2 ##Scale precision for beta 
  prec = matrix(c(PrecParamInt,rep(0,max_cp+2),rep(c(PrecParam,rep(0,max_cp+2)),max_cp+1)),nrow=max_cp+2,ncol=max_cp+2, byrow=TRUE)
  nsubj = dim(Y)[1]
  ntime = dim(Y)[2]
  maxT = max(X)
  minT = min(X)
  minCPmean = unique(sort(X))[2]
  maxCPmean = unique(sort(X, decreasing=TRUE))[2]


  ##Specify JAGS model, without random effects (for initialization) 
  if(cp_prior=='binomial') pw.spec<-textConnection(model.binomial.fixed.oneclass)
  if(cp_prior=='uniform') pw.spec<-textConnection(model.uniform.fixed.oneclass)
  dataList <- list('mean' = mean, 'prec' = prec, 'x' = X, 'y' = Y,
                   'nsubj' = nsubj, 'ntime' = ntime, 'minCPmean' = minCPmean,
                   'maxCPmean' = maxCPmean,'max_cp'=max_cp)
  if(cp_prior=='binomial') dataList$binom_prob <- binom_prob
  if(cp_prior=='uniform') dataList$P_I <- P_I
  cat('Computing initial values...\n')
  Model.jags <- jags.model(pw.spec, data = dataList, n.chains=3, n.adapt=500, quiet=TRUE)
  #update model and save samples
  system.time(update(Model.jags,2000))
  Save.fixed <- jags.samples(Model.jags, c('b', 'cp', 'K','I','sigma2y','C'),
                           2000,thin =4)

  ###Permute changepoint labels if necessary, so they are ordered
  for(i in 1:500){ for(j in 1:3){ for(k in 1:1){
    if(Save.fixed$K[k,i,j]>0){
      ind = c(1:max_cp)[as.logical(Save.fixed$I[k,,i,j])]
      Ord = order(Save.fixed$cp[k,ind,i,j])
      Save.fixed$I[k,,i,j] <- c(Save.fixed$I[k,ind[Ord],i,j],Save.fixed$I[k,-ind,i,j])
      Save.fixed$cp[k,,i,j] = c(Save.fixed$cp[k,ind[Ord],i,j],Save.fixed$cp[k,-ind,i,j])
      Save.fixed$b[k,1:max_cp+2,i,j] = c(Save.fixed$b[k,ind[Ord]+2,i,j],Save.fixed$b[k,c(3:(max_cp+2))[-ind],i,j])
    }
  }}}

  ###Generate initial values::
  Inits = list()
  Inits[[1]] = list(); Inits[[2]] = list(); Inits[[3]] = list()
  Inits[[1]]$mub = Inits[[2]]$mub = Inits[[3]]$mub = matrix(nrow=1,ncol=2+max_cp)
  Inits[[1]]$cp.mean = Inits[[2]]$cp.mean = Inits[[3]]$cp.mean=matrix(nrow=1,ncol=max_cp)
  for(i in 1:3){ for(j in 1:1){
    Inits[[i]]$mub[j,] = rowMeans(Save.fixed$b[j,,Save.fixed$K[j,,i]==getmode(Save.fixed$K[j,,i]),i])
    if(length(Inits[[i]]$cp.mean)<2){
        Inits[[i]]$cp.mean[j, ] = mean(Save.fixed$cp[j,,Save.fixed$K[j, , i] == getmode(Save.fixed$K[j,, i]), i])
      }
      if(length(Inits[[i]]$cp.mean)>1){
      Inits[[i]]$cp.mean[j, ] = rowMeans(Save.fixed$cp[j, 
                                                       , Save.fixed$K[j, , i] == getmode(Save.fixed$K[j, 
                                                                                                      , i]), i])
      }
  }}
  Inits[[1]]$sdb = Inits[[2]]$sdb = Inits[[3]]$sdb = matrix(nrow=1,ncol=max_cp+2)
  Inits[[1]]$cp.sd = Inits[[2]]$cp.sd = Inits[[3]]$cp.sd=matrix(nrow=1,ncol=max_cp)
  for(i in 1:3){ for(j in 1:1){
    for(k in 1:max_cp){ Inits[[i]]$cp.sd[j,k] <- (maxT-minT)/(4*max_cp)}
    Inits[[i]]$sdb[j,1] <- sqrt(1/PrecParamInt)
    for(k in 2:(max_cp+2)) Inits[[i]]$sdb[j,k] <- sqrt(2/PrecParam)
  }}
  
  if(cp_prior=='binomial'){
    Inits[[1]]$I = Inits[[2]]$I = Inits[[3]]$I=matrix(nrow=1,ncol=max_cp)
    for(i in 1:3){ for(j in 1:1){ for(k in 1:max_cp){
      Inits[[i]]$I[j,k] = getmode(Save.fixed$I[j,k,,i])  
    }}}
  } 
  
  if(cp_prior=='uniform'){
    Inits[[1]]$Temp = Inits[[2]]$Temp = Inits[[3]]$Temp =matrix(nrow=1,ncol=max_cp)
    for(i in 1:3){ for(j in 1:1){ for(k in 1:max_cp){
      Inits[[i]]$Temp[j,k] = getmode(Save.fixed$I[j,k,,i])   
    }}}
  } 

  ###Model definition for JAGS (with random effects):
  if(cp_prior=='binomial'&&scale_prior=='uniform') pw.spec<-textConnection(model.binomial.scaleUnif.oneclass)
  if(cp_prior=='uniform'&&scale_prior=='uniform') pw.spec<-textConnection(model.uniform.scaleUnif.oneclass)
  if(cp_prior=='binomial'&&scale_prior=='hc') pw.spec<-textConnection(model.binomial.scaleHC.oneclass)
  if(cp_prior=='uniform'&&scale_prior=='hc') pw.spec<-textConnection(model.uniform.scaleHC.oneclass)
  
  cat('Calibrating MCMC...\n')
  dataList = list('mean' = mean, 'prec' = prec, 'x' = X, 'y' = Y,
                  'nsubj' = nsubj, 'ntime' = ntime, 'maxT' = maxT, 'minT' = minT,
                  'minCPmean' = minCPmean, 'maxCPmean' = maxCPmean,
                  'NumForConv'=5+4*max_cp,'max_cp' = max_cp)
  if(cp_prior=='binomial') dataList$binom_prob<-binom_prob
  if(cp_prior=='uniform') dataList$P_I<-P_I
  Model.jags <- jags.model(pw.spec, data = dataList,inits = Inits,n.chains=3,
                           n.adapt=1000,quiet=TRUE)

  ##run model for  burn in
  cat("Running burn-in...\n")
  system.time(update(Model.jags,iters_BurnIn))

  ##Sample every Thin'th iteration for iters_sampling more iterations
  cat("Collecting samples...\n")
  Save <- jags.samples(Model.jags,
                     c('b', 'mub','cp.mean','muy','cp.sd','K','sigma2y','cp','sdb',
                       'C','ForConv','I','pD','deviance'), iters_sampling,thin = Thin)

  for(i in 1:(length(Save$sigma2y)/3)){ for(j in 1:3){ for(k in 1:1){
    if(Save$K[k,i,j]>0){
      ind = c(1:max_cp)[as.logical(Save$I[k,,i,j])]
      Ord = order(Save$cp.mean[k,ind,i,j])
      Save$cp.mean[k,,i,j] = c(Save$cp.mean[k,ind[Ord],i,j],Save$cp.mean[k,-ind,i,j])
      Save$cp.sd[k,,i,j] = c(Save$cp.sd[k,ind[Ord],i,j],Save$cp.sd[k,-ind,i,j])
      Save$mub[k,3:(max_cp+2),i,j] = c(Save$mub[k,ind[Ord]+2,i,j],Save$mub[k,c(3:(max_cp+2))[-ind],i,j])
      Save$sdb[k,3:(max_cp+2),i,j] = c(Save$sdb[k,ind[Ord]+2,i,j],Save$sdb[k,c(3:(max_cp+2))[-ind],i,j])
    }##add ordering for individual components...
  }}}


###Collect important parameters to assess convergence

  Save$ForConv[1,,] <- Save$sigma2y
  for(k in 1:1){
    Save$ForConv[2:(1+max_cp)+(4+4*max_cp)*(k-1),,] <- Save$cp.mean[k,,,]
    Save$ForConv[(2+max_cp):(1+2*max_cp)+(5+4*max_cp)*(k-1),,] <- Save$cp.sd[k,,,]
    Save$ForConv[(2+2*max_cp):(3+3*max_cp)+(5+4*max_cp)*(k-1),,] <- Save$mub[k,,,]
    Save$ForConv[(4+3*max_cp):(5+4*max_cp)+(5+4*max_cp)*(k-1),,] <- Save$sdb[k,,,]
  }
  mcmcList <- as.mcmc.list(Save$ForConv)
  Gelman.msrf = try(gelman.diag(mcmcList), silent=TRUE)
  RowNames = c('error.var')
  for(k in 1:1){
    RowNames = c(RowNames,paste0('C',k,'_cp',1:max_cp,'.mean'),
               paste0('C',k,'_cp',1:max_cp,'.sd'),  paste0('C',k,'_b',1:(max_cp+2),'.mean'),
               paste0('C',k,'_b',1:(max_cp+2),'.sd'))
  }
  if(class(Gelman.msrf)!="try-error"){row.names(Gelman.msrf$psrf)=RowNames}

  C <- list()
  for(j in 1:1){
    K_prob = matrix(nrow=(max_cp+1),ncol=1)
    for(k in 0:max_cp) {K_prob[k+1]= sum(Save$K[j,,]==k)/length(Save$K[j,,])}
    colnames(K_prob) = 'Probability'
    rownames(K_prob) = c(paste0('K=',0:max_cp))
    K=list()
    for(k in 0:max_cp){
      K[[k+1]] <- list()
      if(K_prob[k+1]>0.01){
        b.array = array(Save$b,dim=c(dim(Save$b)[1],dim(Save$b)[2],dim(Save$b)[3]*dim(Save$b)[4]))
        b.array = b.array[,,Save$K[j,,]==k]
        b = apply(b.array,c(1,2),'mean')[,1:(2+k)]
        mub.array = array(Save$mub[j,,,],dim=c(dim(Save$mub)[2],dim(Save$mub)[3]*dim(Save$mub)[4]))
        mub.array = mub.array[1:(2+k),Save$K[j,,]==k]
        b.mean = apply(mub.array,c(1),'mean')
        b.mean.CI = apply(mub.array,c(1),'quantile',probs=c(0.025,0.975))
        sdb.array = array(Save$sdb[j,,,],dim=c(dim(Save$sdb[j,,,])[1],dim(Save$sdb[j,,,])[2]*dim(Save$sdb[j,,,])[3]))
        sdb.array = sdb.array[1:(2+k),Save$K[j,,]==k]
        b.sd = apply(sdb.array,c(1),'mean')
        b.sd.CI = apply(sdb.array,c(1),'quantile',probs=c(0.025,0.975))
        cp=cp.mean=cp.mean.CI=cp.sd=cp.sd.CI=NULL
        if(k>0){
          cp.array = array(Save$cp,dim=c(dim(Save$cp)[1],dim(Save$cp)[2],dim(Save$cp)[3]*dim(Save$cp)[4]))
          cp.array = cp.array[,,Save$K[j,,]==k]
          cp = apply(cp.array,c(1,2),'mean')[,1:k]
          cpmean.array = array(Save$cp.mean[j,,,],dim=c(dim(Save$cp.mean)[2],dim(Save$cp.mean)[3]*dim(Save$cp.mean)[4]))
          cpmean.array = cpmean.array[1:k,Save$K[j,,]==k,drop=FALSE]
          cp.mean = apply(cpmean.array,c(1),'mean')
          cp.mean.CI = apply(cpmean.array,c(1),'quantile',probs=c(0.025,0.975))
          cpsd.array = array(Save$cp.sd[j,,,],dim=c(dim(Save$cp.sd)[2],dim(Save$cp.sd)[3]*dim(Save$cp.sd)[4]))
          cpsd.array = cpsd.array[1:k,Save$K[j,,]==k,drop=FALSE]
          cp.sd = apply(cpsd.array,c(1),'mean')
          cp.sd.CI = apply(cpsd.array,c(1),'quantile',probs=c(0.025,0.975))
        }
        K[[k+1]] = list('b' = b, 'b.mean' = b.mean, 'b.mean.CI' = b.mean.CI, 'b.sd' = b.sd, 
                      'b.sd.CI' = b.sd.CI, 'cp' = cp, 'cp.mean' = cp.mean, 'cp.mean.CI'=cp.mean.CI,'cp.sd'=cp.sd,'cp.sd.CI'=cp.sd.CI)
      }
      C[[j]] <- list('K_prob'=K_prob,'K'=K)
    }}
  ##return mean values
  DIC <- mean(Save$deviance)+mean(Save$pD)
  Results<-list('Gelman.msrf' = Gelman.msrf,
                'DIC'=DIC,
                'y.mean' = summary(Save$muy, FUN = 'mean')[[1]],
                'K_prob'=K_prob,'K'=K,'error.sd'=mean(sqrt(Save$sigma2y)),'error.sd.CI'= quantile(sqrt(Save$sigma2y),probs=c(0.025,0.975)))
  if(SaveChains==TRUE){Results$Chains=Save}
  class(Results)=PGM
  return(Results)
}




################################Jags models
model.binomial.fixed.oneclass <- "model {
  for( i in 1 : nsubj ) {
C[i] <- 1
for( j in 1 : ntime ) {
y[i , j] ~ dnorm(b[C[i],1] + b[C[i],2]*x[i,j]+ sum(I[C[i],1:max_cp]*b[C[i],3:(2+max_cp)]*X.cp[C[i],1:max_cp,i,j]),tauy)
} }
##ncomponents components
for(k in 1 : 1){
for(l in 1 : max_cp){
I[k,l] ~dbern(binom_prob) 
cp[k,l] ~ dunif(minCPmean,maxCPmean)
for(j in 1 : ntime){ for(i in 1 : nsubj){
X.cp[k,l,i,j] <- max(0,x[i , j]-cp[k,l])
}}}
K[k] <- sum(I[k,1:max_cp])
# prior precision for y
#prior distribution of the fixed parameters
b[k,1:(max_cp+2)]~dmnorm(mean[1:(max_cp+2)],prec[1:(max_cp+2),1:(max_cp+2)])
###Order changepoints for saving
##    Ord[k,1:max_cp] <- ifelse(K[k]<max_cp,c(rank(cp[k,1:K[k]]), (K[k]+1):max_cp),rank(cp[k,1:max_cp]))
##    Ord[k,1:max_cp] <- rank(cp[k,1:max_cp])
##    cp.ord[k,1:max_cp] <- cp[k,Ord[k,1:max_cp]] 
##    b.cps[k,1:max_cp] <- b[k,3:(max_cp+2)]
##    b.ord[k,3:(max_cp+2)] <- b.cps[k,Ord[k,1:max_cp]]
}
sigma2y <- 1/tauy
tauy ~ dgamma(0.001,0.001)
}
"

model.binomial.scaleUnif.oneclass <- "model {
for( i in 1 : nsubj ) {
C[i] <-1
for( j in 1 : ntime ) {
y[i , j] ~ dnorm(muy[i,j],tauy)
muy[i , j] <- b[i,1] + b[i,2]*x[i,j]+ sum(I[C[i],1:max_cp]*b[i,3:(max_cp+2)]*X.cp[C[i],1:max_cp,i,j])
} }
# distribution for the random-effects parameters
for( i in 1 : nsubj ) {
for(j in 1:(max_cp+2)){
b[i,j]~dnorm(mub[C[i],j], taub[C[i],j])
}
for(l in 1:max_cp){
cp[i,l]~dnorm(cp.mean[C[i],l],cp.prec[C[i],l]) T(minT,maxT)
}}
for(k in 1 : 1){
for(l in 1 : max_cp){
I[k,l] ~dbern(binom_prob) 
cp.mean[k,l] ~ dunif(minCPmean,maxCPmean)
cp.prec[k,l] <- 1/(cp.sd[k,l]^2)
cp.sd[k,l] ~ dunif(0,(maxT-minT)/4)
for(j in 1 : ntime){ for(i in 1 : nsubj){
X.cp[k,l,i,j] <- max(0,x[i , j]-cp[i,l])}}
}
K[k] <- sum(I[k,1:max_cp])
#prior distribution of the fixed parameters
mub[k,1:(max_cp+2)]~dmnorm(mean[1:(max_cp+2)],prec[1:(max_cp+2),1:(max_cp+2)])
for(j in 1:(max_cp+2)){
sdb[k,j] ~ dunif(0,(2/prec[j,j])^0.5)
taub[k,j] <- 1/(sdb[k,j]^2) 
}}
ForConv[1:NumForConv] <- c(rep(0,NumForConv)) ##to be filled in later
tauy ~ dgamma(0.001,0.001)
sigma2y <- 1/tauy
}
"

model.binomial.scaleHC.oneclass <- "model {
for( i in 1 : nsubj ) {
C[i] <-1
for( j in 1 : ntime ) {
y[i , j] ~ dnorm(muy[i,j],tauy)
muy[i , j] <- b[i,1] + b[i,2]*x[i,j]+ sum(I[C[i],1:max_cp]*b[i,3:(max_cp+2)]*X.cp[C[i],1:max_cp,i,j])
} }
# distribution for the random-effects parameters
for( i in 1 : nsubj ) {
for(j in 1:(max_cp+2)){
b[i,j]~dnorm(mub[C[i],j], taub[C[i],j])
}
for(l in 1:max_cp){
cp[i,l]~dnorm(cp.mean[C[i],l],cp.prec[C[i],l]) T(minT,maxT)
}}
for(k in 1 : 1){
for(l in 1 : max_cp){
I[k,l] ~dbern(binom_prob) 
cp.mean[k,l] ~ dunif(minCPmean,maxCPmean)
cp.prec[k,l] <- 1/(cp.sd[k,l]^2)
cp.sd[k,l] ~ dt(0, 1/(((maxT-minT)/4)/tan(0.45*3.1416))^2, 1)T(0,)
for(j in 1 : ntime){ for(i in 1 : nsubj){
X.cp[k,l,i,j] <- max(0,x[i , j]-cp[i,l])}}
}
K[k] <- sum(I[k,1:max_cp])
#prior distribution of the fixed parameters
mub[k,1:(max_cp+2)]~dmnorm(mean[1:(max_cp+2)],prec[1:(max_cp+2),1:(max_cp+2)])
for(j in 1:(max_cp+2)){
sdb[k,j] ~ dt(0, 1/(((2/prec[j,j])^0.5)/tan(0.45*3.1416))^2, 1)T(0,)
taub[k,j] <- 1/(sdb[k,j]^2) 
}}
ForConv[1:NumForConv] <- c(rep(0,NumForConv)) ##to be filled in later
tauy ~ dgamma(0.001,0.001)
sigma2y <- 1/tauy
}
"

model.uniform.fixed.oneclass <- "model {
  for( i in 1 : nsubj ) {
C[i] <- 1
for( j in 1 : ntime ) {
y[i , j] ~ dnorm(b[C[i],1] + b[C[i],2]*x[i,j]+ sum(I[C[i],1:max_cp]*b[C[i],3:(2+max_cp)]*X.cp[C[i],1:max_cp,i,j]),tauy)
} }
##ncomponents components
for(k in 1 : 1){
for(l in 1 : max_cp){
Temp[k,l] ~dbern(P_I[l])
I[k,l] <- prod(Temp[k,1:l])
cp[k,l] ~ dunif(minCPmean,maxCPmean)
for(j in 1 : ntime){ for(i in 1 : nsubj){
X.cp[k,l,i,j] <- max(0,x[i , j]-cp[k,l])
}}}
K[k] <- sum(I[k,1:max_cp])
# prior precision for y
#prior distribution of the fixed parameters
b[k,1:(max_cp+2)]~dmnorm(mean[1:(max_cp+2)],prec[1:(max_cp+2),1:(max_cp+2)])
}
sigma2y <- 1/tauy
tauy ~ dgamma(0.001,0.001)
}
"

model.uniform.scaleUnif.oneclass <- "model {
for( i in 1 : nsubj ) {
C[i] <-1
for( j in 1 : ntime ) {
y[i , j] ~ dnorm(muy[i,j],tauy)
muy[i , j] <- b[i,1] + b[i,2]*x[i,j]+ sum(I[C[i],1:max_cp]*b[i,3:(max_cp+2)]*X.cp[C[i],1:max_cp,i,j])
} }
# distribution for the random-effects parameters
for( i in 1 : nsubj ) {
for(j in 1:(max_cp+2)){
b[i,j]~dnorm(mub[C[i],j], taub[C[i],j])
}
for(l in 1:max_cp){
cp[i,l]~dnorm(cp.mean[C[i],l],cp.prec[C[i],l]) T(minT,maxT)
}}

for(k in 1 : 1){
for(l in 1 : max_cp){
Temp[k,l] ~dbern(P_I[l])
I[k,l] <- prod(Temp[k,1:l])
cp.mean[k,l] ~ dunif(minCPmean,maxCPmean)
cp.prec[k,l] <- 1/(cp.sd[k,l]^2)
cp.sd[k,l] ~ dunif(0,(maxT-minT)/4)
for(j in 1 : ntime){ for(i in 1 : nsubj){
X.cp[k,l,i,j] <- max(0,x[i , j]-cp[i,l])}}
}
K[k] <- sum(I[k,1:max_cp])
#prior distribution of the fixed parameters
mub[k,1:(max_cp+2)]~dmnorm(mean[1:(max_cp+2)],prec[1:(max_cp+2),1:(max_cp+2)])
for(j in 1:(max_cp+2)){
sdb[k,j] ~ dunif(0,(2/prec[j,j])^0.5)
taub[k,j] <- 1/(sdb[k,j]^2) 
}}

ForConv[1:NumForConv] <- c(rep(0,NumForConv)) ##to be filled in later
tauy ~ dgamma(0.001,0.001)
sigma2y <- 1/tauy
}
"

model.uniform.scaleHC.oneclass <- "model {
for( i in 1 : nsubj ) {
C[i] <-1
for( j in 1 : ntime ) {
y[i , j] ~ dnorm(muy[i,j],tauy)
muy[i , j] <- b[i,1] + b[i,2]*x[i,j]+ sum(I[C[i],1:max_cp]*b[i,3:(max_cp+2)]*X.cp[C[i],1:max_cp,i,j])
} }
# distribution for the random-effects parameters
for( i in 1 : nsubj ) {
for(j in 1:(max_cp+2)){
b[i,j]~dnorm(mub[C[i],j], taub[C[i],j])
}
for(l in 1:max_cp){
cp[i,l]~dnorm(cp.mean[C[i],l],cp.prec[C[i],l]) T(minT,maxT)
}}

for(k in 1 : 1){
for(l in 1 : max_cp){
Temp[k,l] ~dbern(P_I[l])
I[k,l] <- prod(Temp[k,1:l])
cp.mean[k,l] ~ dunif(minCPmean,maxCPmean)
cp.prec[k,l] <- 1/(cp.sd[k,l]^2)
cp.sd[k,l] ~ dt(0, 1/(((maxT-minT)/4)/tan(0.45*3.1416))^2, 1)T(0,)
for(j in 1 : ntime){ for(i in 1 : nsubj){
X.cp[k,l,i,j] <- max(0,x[i , j]-cp[i,l])}}
}
K[k] <- sum(I[k,1:max_cp])
#prior distribution of the fixed parameters
mub[k,1:(max_cp+2)]~dmnorm(mean[1:(max_cp+2)],prec[1:(max_cp+2),1:(max_cp+2)])
for(j in 1:(max_cp+2)){
sdb[k,j] ~ dt(0, 1/(((2/prec[j,j])^0.5)/tan(0.45*3.1416))^2, 1)T(0,)
taub[k,j] <- 1/(sdb[k,j]^2) 
}}

ForConv[1:NumForConv] <- c(rep(0,NumForConv)) ##to be filled in later
tauy ~ dgamma(0.001,0.001)
sigma2y <- 1/tauy
}
"

################################Function to get mode (for initializing number of changepoints)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



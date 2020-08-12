plotPGMM <- function(X, Y, Results=NA, xlab='X', ylab='Y', Colors=NULL, MeanColors=NULL, ...){
  
  N= dim(Y)[1]
  M = dim(X)[2]
  xvec = seq(min(X),max(X),length.out=100)
if(!is.list(Results)){
    plot(X[1,!is.na(Y[1,]) ], Y[1, !is.na(Y[1,])], type = "l", ylim = c(min(Y,na.rm=TRUE), max(Y,na.rm=TRUE)),xlim= c(min(X,na.rm=TRUE), max(X,na.rm=TRUE)), 
         xlab = xlab, ylab = ylab, ...)
    for (i in 2:N) lines(X[i,!is.na(Y[i,])], Y[i,!is.na(Y[i,])], type = "l")
  return('Spaghetti plot with no PGMM results')
}
  
if(is.list(Results)){
  UniqClust <- unique(Results$Clust)
  n_clust = length(UniqClust)
  Clus=list()
  for(p in 1:n_clust){
    k<-UniqClust[p]
    Clus[[k]] <- c(1:N)[Results$Clust==k]}
  Clus_mean = list()
  max_cp <- length(Results$C[[1]][[2]])-1
  for(p in 1:n_clust){
    k<-UniqClust[p]
  for(i in 1:(max_cp+1)){
  if(which.max(Results$C[[k]]$K_prob)==i){
    if(i==1) I=rep(0,max_cp)
    if(i>1) I=c(rep(1,i-1),rep(0,max_cp-i+1))
    int = Results$C[[k]]$K[[i]]$b.mean[1]
    b1 = Results$C[[k]]$K[[i]]$b.mean[2]
    b.cp = rep(0,max_cp)
    cp = rep(0,max_cp)
    if(i>1){ 
      b.cp[1:(i-1)] = Results$C[[k]]$K[[i]]$b.mean[3:(1+i)] 
      cp[1:(i-1)] = Results$C[[k]]$K[[i]]$cp.mean[1:(i-1)]
    }
    ##Compute fixed-effect model for cluster i:
    Clus_mean[[k]] = rep(0,100)
    for(j in 1:100){ 
      temp <- int  + b1*xvec[j]
      if(i>1){
        for(l in 1:(i-1)) temp = temp+b.cp[l]*(max(0,xvec[j]-cp[l]))
      }
    Clus_mean[[k]][j] <- temp
  }
  }  
}}
  if(is.null(Colors)) Colors = c('blue','red','green','gold','gray')
  if(is.null(MeanColors)) MeanColors = c('darkblue','darkred','darkgreen','gold4','darkgray')
    plot(0, 0, xlim = c(min(X), max(X)), ylim = c(min(Y,na.rm=TRUE), 
                                                  max(Y,na.rm=TRUE)), ylab = ylab, xlab = xlab, ...)
  for(p in 1:n_clust){
    k <- UniqClust[p]
  for(i in Clus[[k]]){
    points(X[i,!is.na(Y[i,])], Y[i,!is.na(Y[i,])], type = "l", col=Colors[k])
  }
  points(xvec,Clus_mean[[k]], type='l', col=MeanColors[k],lwd=4)
}
}}



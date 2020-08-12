plotPGM <- function(X, Y, Results=NA, xlab='X', ylab='Y', ...){
  
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
  max_cp <- length(Results$K)-1
  for(i in 1:(max_cp+1)){
  if(which.max(Results$K_prob)==i){
    if(i==1) I=rep(0,max_cp)
    if(i>1) I=c(rep(1,i-1),rep(0,max_cp-i+1))
    int = Results$K[[i]]$b.mean[1]
    b1 = Results$K[[i]]$b.mean[2]
    b.cp = rep(0,max_cp)
    cp = rep(0,max_cp)
    if(i>1){ 
      b.cp[1:(i-1)] = Results$K[[i]]$b.mean[3:(1+i)] 
      cp[1:(i-1)] = Results$K[[i]]$cp.mean[1:(i-1)]
    }
    ##Compute fixed-effect model for cluster i:
    Clus_mean = rep(0,M)
    for(j in 1:100){ 
      temp <- int  + b1*xvec[j]
      if(i>1){
        for(l in 1:(i-1)) temp = temp+b.cp[l]*(max(0,xvec[j]-cp[l]))
      }
      Clus_mean[j] <- temp
  }
  }  
}
    plot(0, 0, xlim = c(min(X), max(X)), ylim = c(min(Y,na.rm=TRUE), 
                                                  max(Y,na.rm=TRUE)), ylab = ylab, xlab = xlab, ...)
    for (i in 1:dim(X)[1]) {
      points(X[i,!is.na(Y[i,])], Y[i,!is.na(Y[i,])], type = "l", col = "darkgray")
    }
  points(xvec,Clus_mean, type='l',lwd=4)
}
}

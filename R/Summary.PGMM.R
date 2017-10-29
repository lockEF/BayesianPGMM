summary.PGMM <- function (object, ...) {
  n_clust <- length(object$ClusProb)
  Changepoints <- c()
  for(i in 1:n_clust) Changepoints[i] <- which.max(object$C[[i]]$K_prob)-1
  maxCP <- max(Changepoints)
  CoeffMat <- matrix(rep(NA,4*(maxCP+1)*n_clust), nrow=4+4*maxCP,ncol=n_clust)
  for(i in 1:n_clust){
    CoeffMat[1,i] <- object$C[[i]]$K[[Changepoints[i]+1]]$b.mean[1] 
    CoeffMat[2,i] <- object$C[[i]]$K[[Changepoints[i]+1]]$b.sd[1] 
    CoeffMat[3,i] <- object$C[[i]]$K[[Changepoints[i]+1]]$b.mean[2] 
    CoeffMat[4,i] <- object$C[[i]]$K[[Changepoints[i]+1]]$b.sd[2]
    if(Changepoints[i]>0){
      for(k in 1:Changepoints[i]){
      CoeffMat[4*k+1,i] <- object$C[[i]]$K[[Changepoints[i]+1]]$cp.mean[k]
      CoeffMat[4*k+2,i] <- object$C[[i]]$K[[Changepoints[i]+1]]$cp.sd[k]
      CoeffMat[4*k+3,i] <- object$C[[i]]$K[[Changepoints[i]+1]]$b.mean[k+2]
      CoeffMat[4*k+4,i] <- object$C[[i]]$K[[Changepoints[i]+1]]$b.sd[k+2]
      }}
  }
  ClusData <- rbind(t(object$ClusProb),Changepoints,CoeffMat)
  colnames(ClusData) <- paste("Class",c(1:n_clust))
  RowNames <- rep('',6+4*maxCP)
  RowNames[1] = "Probability"
  RowNames[2] = "Changepoints"
  RowNames[3] = "Mean intercept"
  RowNames[4] = "Intercept SD"
  RowNames[5] = "Mean slope"
  RowNames[6] = "Slope SD"
  for(k in 1:maxCP){
    RowNames[2+4*k+1] = paste("Mean changepoint",k)
    RowNames[2+4*k+2] = paste("Changepoint",k,"SD")
    RowNames[2+4*k+3] = paste("Mean slope change",k)
    RowNames[2+4*k+4] = paste("Slope change SD", k)
  }  
    
  rownames(ClusData) <- RowNames
  print(ClusData,digits=3)
  cat("\n")
  cat("Gelman's msrf:", object$Gelman.msrf$mpsrf, "\n")
  cat("DIC:", object$DIC)
  return(invisible(ClusData))
}
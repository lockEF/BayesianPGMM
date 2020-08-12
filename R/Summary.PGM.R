summary.PGM <- function (object, ...) {
  n_cp <- length(object$K_prob)
  Changepoints <- c()
  for(i in 1:(n_cp)) Changepoints[i] <- i-1
  maxCP <- n_cp-1
  CoeffMat <- matrix(rep(NA,4*(maxCP+1)*n_cp), nrow=4+4*maxCP,ncol=n_cp)
  for(i in 1:n_cp){
    if(object$K_prob[i]>0){
    CoeffMat[1,i] <- object$K[[i]]$b.mean[1] 
    CoeffMat[2,i] <- object$K[[i]]$b.sd[1] 
    CoeffMat[3,i] <- object$K[[i]]$b.mean[2] 
    CoeffMat[4,i] <- object$K[[i]]$b.sd[2]
    if(Changepoints[i]>0){
      for(k in 1:Changepoints[i]){
        CoeffMat[4*k+1,i] <- object$K[[i]]$cp.mean[k]
        CoeffMat[4*k+2,i] <- object$K[[i]]$cp.sd[k]
        CoeffMat[4*k+3,i] <- object$K[[i]]$b.mean[k+2]
        CoeffMat[4*k+4,i] <- object$K[[i]]$b.sd[k+2]
      }}
    }
  }
  ClusData <- rbind(t(object$K_prob),CoeffMat)
  colnames(ClusData) <- paste(c(0:maxCP), "CPs")
  RowNames <- rep('',5+4*maxCP)
  RowNames[1] = "Probability"
  RowNames[2] = "Mean intercept"
  RowNames[3] = "Intercept SD"
  RowNames[4] = "Mean slope"
  RowNames[5] = "Slope SD"
  for(k in 1:maxCP){
    RowNames[1+4*k+1] = paste("Mean changepoint",k)
    RowNames[1+4*k+2] = paste("Changepoint",k,"SD")
    RowNames[1+4*k+3] = paste("Mean slope change",k)
    RowNames[1+4*k+4] = paste("Slope change SD", k)
  }  
  
  rownames(ClusData) <- RowNames
  print(ClusData,digits=3)
  cat("\n")
  cat("Gelman's msrf:", object$Gelman.msrf$mpsrf, "\n")
  cat("DIC:", object$DIC)
  return(invisible(ClusData))
}
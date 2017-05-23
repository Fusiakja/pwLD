`estCredibleInterval` <-
function( hapCounts, nSim=1000, Dir=rep(1,9), alpha=.05, LD="Q", tol=.Machine$double.eps^.6, digits=12, HSweight =4)
{
  ci <- matrix(0, nrow=length(LD), ncol=2)
  rownames(ci) <- LD
  
  colnames(ci) <- c( paste(alpha/2*100, "%", sep=""), paste( (1- (alpha/2))*100, "%", sep="" ) )
  
  res <- .C("credibleInterval", genoCounts=as.integer(t(hapCounts)),  nSim=as.integer(nSim), Dir=as.double(Dir), LD=as.character(LD), LDnumb=as.integer(length(LD)), tol=as.double(tol), digits=as.integer(digits), LDdist=as.double(numeric(length(LD) * nSim) ), HSweight=as.double(HSweight), alpha=as.double(alpha), cilow=as.double(rep(0,length(LD))), ciup=as.double(rep(0,length(LD))),PACKAGE="pwLD")
  
  # the distribution of the LD measures 
  LDdist <- matrix(res$LDdist, nrow=length(LD), byrow=T)
  rownames(LDdist) <- LD
  
  # determine the shortest interval containing (1-alpha)*100 percent of the distribution
  #for(i in LD)
    #ci[i, ] <- shortest.cred.int(LDdist[i,],  alpha)
  #ci[i, ] <- quantile(LDdist[i, ], c(alpha/2, 1- (alpha/2)) )
  
  
  # the distribution of the LD measures 
  LDdist <- matrix(res$LDdist, nrow=length(LD), byrow=T)
  rownames(LDdist) <- LD
  
  # confidence intervals
  for(i in 1:length(LD))
  {
    ci[i, ] <- c(res$cilow[i], res$ciup[i])
  }
  output <- vector("list")
  output[[1]] <- LDdist
  output[[2]] <- ci

  
  return(output)

}


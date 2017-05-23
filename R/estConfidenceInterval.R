`estConfidenceInterval` <-
  function(hapCounts, hapFreq, N, nSim= 1000, alpha=.05, LD="Dprime", tol=.Machine$double.eps^.6, digits=12, HSweight=4, strategy = "bootstrap", paradigm="freq", Dir=c(1,1,1,1), intervall=c(0,1))
  {	
   
    
    alpha=alpha
    ci <- matrix(0, nrow=length(LD), ncol=2)
    rownames(ci) <- LD
    
    
    colnames(ci) <- c( paste(alpha/2*100, "%", sep=""), paste( (1- (alpha/2))*100, "%", sep="" ) )
    
    if(strategy=="bootstrap")
    {  
      res <- .C("confidenceGenoInterval", genoCounts=as.double(t(hapCounts)), genoFreq=as.double(t(hapFreq)), N=as.integer(N), paradigm=as.character(paradigm),nSim=as.integer(nSim), LD=as.character(LD), LDnumb=as.integer(length(LD)), tol=as.double(tol), digits=as.integer(digits), LDdist=as.double(numeric(length(LD) * nSim) ), HSweight=as.double(HSweight), alpha=as.double(alpha), strategy=as.character(strategy), cilow=as.double(rep(0,length(LD))), ciup=as.double(rep(0,length(LD))), tables=matrix(rep(0,nSim*9), nrow = 9, ncol = nSim,byrow = T),Dir=as.double(Dir),vars=as.double(rep(0,length(LD))),intervall=as.integer(intervall),PACKAGE="pwLD")
    }  
    
    if(strategy=="jackknife")
    {  
      res <- .C("confidenceGenoInterval", genoCounts=as.double(t(hapCounts)), genoFreq=as.double(t(hapFreq)), N=as.integer(N), paradigm=as.character(paradigm),nSim=as.integer(nSim), LD=as.character(LD), LDnumb=as.integer(length(LD)), tol=as.double(tol), digits=as.integer(digits), LDdist=as.double(numeric(length(LD) * N) ), HSweight=as.double(HSweight), alpha=as.double(alpha), strategy=as.character(strategy), cilow=as.double(rep(0,length(LD))), ciup=as.double(rep(0,length(LD))), tables=matrix(rep(0,9), nrow = 9, ncol = 9,byrow = T),Dir=as.double(Dir),vars=as.double(rep(0,length(LD))),intervall=as.integer(intervall),PACKAGE="pwLD")
    }  
    # the distribution of the LD measures 
    LDdist <- matrix(res$LDdist, nrow=length(LD), byrow=T)
    rownames(LDdist) <- LD
    
    # confidence intervals
    for(i in 1:length(LD))
    {
      ci[i, ] <- c(res$cilow[i], res$ciup[i])
    }
    output <- vector("list")
    Erg <- list()
    Erg[[1]] <- res$tables
    Erg[[2]] <- LDdist
    Erg[[3]] <- res$vars
    output[[1]] <- Erg
    output[[2]] <- ci
    
    return(output)	
    
  }


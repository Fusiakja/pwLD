`estHapFreq` <-
  function(hapCounts, Neff, paradigm="freq" ,dirich = rep(1,4), mc, seed=F )
  {

    res=.C("estimateFrequencies", genoFreq=as.double(t(hapCounts)), haploFreq=as.double(matrix(0, nrow=3, ncol=3)), Neff=as.integer(Neff), Dir=as.double(dirich), mc=as.integer(mc), paradigm=as.character(paradigm), PACKAGE="pwLD" )
    
    
    output <- list()
    output[[1]] <- matrix(res$haploFreq, nrow=3, byrow=T, dimnames=list(c("0", "1", "Sum"), c("0", "1", "Sum")))
    
    #output[[1]] <- matrix(c(res$hapoo, res$hapo1, res$hapoo+res$hapo1, res$hap1o, res$hap11, res$hap1o+res$hap11, res$hapoo+res$hap1o, res$hap1o+res$hap11, res$hapoo+res$hapo1+res$hap1o+res$hap11), ncol=3, nrow=3, byrow=T)#, dimnames=list(c("0", "1", "Sum"), c("0", "1", "Sum")))
    # if (paradigm == "freq") 
    # {
    #  hapFreq <-  hapCounts/N 
    # }
    # if (paradigm == "bayes") 
    # {
    #   hapFreq[1,1] <- (dirich[1]+hapCounts[1,1])/(4*dirich[1]+N)
    #   hapFreq[1,2] <- (dirich[2]+hapCounts[1,2])/(4*dirich[2]+N)
    #   hapFreq[2,1] <- (dirich[3]+hapCounts[2,1])/(4*dirich[3]+N)
    #   hapFreq[2,2] <- (dirich[4]+hapCounts[2,2])/(4*dirich[4]+N)
    #   
    #   
    #   hapFreq[3,1]=hapFreq[1,1]+hapFreq[2,1]
    #   hapFreq[3,2]=hapFreq[1,2]+hapFreq[2,2]
    #   hapFreq[1,3]=hapFreq[1,1]+hapFreq[1,2]
    #   hapFreq[2,3]=hapFreq[2,1]+hapFreq[2,2]
    #   hapFreq[3,3]=hapFreq[1,3]+hapFreq[2,3]
    # }
    # 
    # if (paradigm == "fullbayes") 
    # {
    #   d = c(dirich[1]+hapCounts[1,1], dirich[2]+hapCounts[1,2], dirich[3]+hapCounts[2,1], dirich[4]+hapCounts[2,2])
    #   if (seed == T) 
    #   {
    #     set.seed(42)
    #   }
    #   fb=colMeans(rdirichlet(mc,d))
    #   hapFreq[1,1]=fb[1]
    #   hapFreq[1,2]=fb[2]
    #   hapFreq[2,1]=fb[3]
    #   hapFreq[2,2]=fb[4]
    #   
    #   hapFreq[3,1]=hapFreq[1,1]+hapFreq[2,1]
    #   hapFreq[3,2]=hapFreq[1,2]+hapFreq[2,2]
    #   hapFreq[1,3]=hapFreq[1,1]+hapFreq[1,2]
    #   hapFreq[2,3]=hapFreq[2,1]+hapFreq[2,2]
    #   hapFreq[3,3]=hapFreq[1,3]+hapFreq[2,3]
    # }

    
    
    return(output[1])
  }
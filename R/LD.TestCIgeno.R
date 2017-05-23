`LD.TestCI_geno` <-
  function(hapt ,snps, code=c(0,1,2,3), LD=c("Dprime", "Q", "r", "Y", "HS"), HSweight=4 , CI=T, strategy=c("jackknife","bootstrap", "zapata"), alpha=.05, n.sim=5000, returnLDdist=F, paradigm=c("freq", "bayes", "fullbayes"), dirich=rep(1,9) , all.solutions=F, tol=.Machine$double.eps^.6, digits=12, seed=F, mc=100000, counts=matrix(0, nrow=4, ncol=4, byrow=T), Ns=100, gen_hap=c("geno","hap"), intervall=c(0,1))
  {
    paradigm <- match.arg(paradigm)
    
    if(class(snps) == "data.frame") snps <- as.matrix(snps)
    if(dim(snps)[1] < 2) stop("\n'snps' must be a 2 x N matrix!! \n\n")
    if(dim(snps)[1] > 2) stop("\nMore than two SNPs detected! Please refer to function 'LD.all'!\n\n")
    if(is.null(rownames(snps))) rownames(snps) <- c("snp1","snp2")
    #genotypes <- genotypes.3x3(snps, code=code, method=paradigm, dirich=dirich)
    
    # estimate genotypic frequencies	
    genotypes <- estGeno(snps, code=code, paradigm=paradigm, dirich=dirich, mc=mc)
    genotypes.freqs <- hapt
    #set.seed(1)
    gent <- rmultinom(1, Ns/2,c(hapt[1,1],hapt[1,2],hapt[1,3],hapt[2,1],hapt[2,2],hapt[2,3],hapt[3,1],hapt[3,2],hapt[3,3]))
    genotypes.counts <- matrix(c(gent[1],gent[2],gent[3],sum(gent[1],gent[2],gent[3]),gent[4],gent[5],gent[6],sum(gent[4],gent[5],gent[6]),gent[7],gent[8],gent[9],sum(gent[7],gent[8],gent[9]),sum(gent[1],gent[4],gent[7]),sum(gent[2],gent[5],gent[8]), sum(gent[3],gent[6],gent[9]),sum(gent)), ncol = 4,nrow = 4,byrow = TRUE)
    #rownames(genocounts) <- c("00", "01", "11", "Sum")
    #colnames(genocounts) <- c("00", "01", "11", "Sum")
    
    genotypes.freqs <- genotypes.counts/(Ns/2)
    
    
    #genocounts <- hapt#matrix(c(rep(0,16)), nrow = 4, ncol = 4)
    #genotypes.counts <- genotypes$count# hapt*Ns/2
    #rownames(genocounts) <- c("00", "01", "11", "Sum")
    #colnames(genocounts) <- c("00", "01", "11", "Sum")
    
    
    # if(!(((genotypes.freqs[4,1] > 0)&&(genotypes.freqs[4,1]<genotypes.freqs[4,4]))&&((genotypes.freqs[1,4]>0)&&(genotypes.freqs[1,4]<genotypes.freqs[4,4]))))
    # {
    #   #warning("Allele frequencies of 0 or 1 !!")
    #   output <- vector("list", length(LD)+1)
    #   output[[1]] <- genotypes.freqs
    #   for(i in 1:length(LD))
    #     output[[i+1]] <- NA
    #   names(output) <- c("table", LD)
    # 
    #   return(output)
    # }
    #genotypes.counts <- genotypes$count
    # # # sample size without missings
    Neff <- 0.5*Ns 
    # # 
    # # # estimate allele frequencies
    # # af<- estAF(genotypes.counts)
    # 
    # 
    # 
    # # estimate pairwise haplotype frequencies
    res <- estHap( genotypes.freqs, tol=tol, digits=digits, Neff = Neff)
    hapCounts <- res$table*2*Neff
    # hapCounts <- matrix(rep(0,9), nrow = 3, ncol = 3)
    # hapCounts[1,1] <- round(res$table[1,1]*2*Neff, digits = 0)#hapt
    # hapCounts[2,1] <- round(res$table[2,1]*2*Neff, digits = 0)
    # hapCounts[1,2] <- round(res$table[1,2]*2*Neff, digits = 0)
    # hapCounts[2,2] <- round(res$table[2,2]*2*Neff, digits = 0)
    # hapCounts[3,1] <- hapCounts[1,1]+hapCounts[2,1]
    # hapCounts[1,3] <- hapCounts[1,1]+hapCounts[1,2]
    # hapCounts[3,2] <- hapCounts[1,2]+hapCounts[2,2]
    # hapCounts[2,3] <- hapCounts[2,1]+hapCounts[2,2]
    # hapCounts[3,3] <- sum(hapCounts[1,1],hapCounts[1,2], hapCounts[2,1],hapCounts[2,2])
    # 
    # hapCounts <- hapt
    allSol <- res$solutions
    # #hapCounts <- hapt
    
    if(seed==TRUE)
    {
      set.seed(1)
    }
    
    hapFreq <- estHapFreq(hapCounts = hapCounts, paradigm = paradigm, Neff=hapCounts[3,3], dirich=dirich, mc=mc, seed=seed)
    hapFreq <- hapFreq[[1]]
    
    # #if any allele frequency is zero or one, all LD measures are set to 'NA'
    if(!(((hapFreq[3,1] > 0)&&(hapFreq[3,1]<1))&&((hapFreq[1,3]>0)&&(hapFreq[1,3]<1))))
    {
      #warning("Allele frequencies of 0 or 1 !!")
      output <- vector("list", length(LD)+1)
      output[[1]] <- hapFreq
      for(i in 1:length(LD))
        output[[i+1]] <- NA
      names(output) <- c("table", LD)
      
      return(output)
    }
    
    
    # estimate LD measures
    LDm <- estLD(hapFreq, LD, HSweight = HSweight)
    
    
    if(CI)
    {
      if(seed==TRUE)
      {
        set.seed(1)
      }
      #if(paradigm == "freq")
      if(gen_hap=="hap")
      {      
        ci.res <- estConfidenceInterval(hapCounts = hapCounts, hapFreq = hapFreq, N=Neff, paradigm = paradigm,nSim= n.sim, LD=LD,  tol=tol , alpha=alpha, digits=digits, strategy=strategy, Dir=dirich, intervall=intervall)
      }
      else
      {
        ci.res <- kestConfidenceInterval(hapCounts = genotypes.counts, hapFreq = genotypes.freqs, N=Neff, paradigm = paradigm,nSim= n.sim, LD=LD,  tol=tol , alpha=alpha, digits=digits, strategy=strategy, Dir=dirich, intervall=intervall)
      }
      #if(paradigm == "bayes" || paradigm == "fullbayes")
      #ci.res <- estCredibleInterval(hapCounts, nSim=n.sim, LD=LD, Dir=dirich, tol=tol, alpha=alpha, digits=digits)
    }
    
    #generate the output
    output <- list()
    output[[1]] <- hapFreq					# table of haplotype frequencies
    output[[2]] <- estLD(hapFreq*Neff, "chi2")			# chi square statisitcs, deviation from complete linkage equilibrium
    output[[3]] <- Neff						# sample size, without missings
    output[[4]] <- genotypes.counts[3,3]
    output[[5]] <- genotypes$SNP
    names(output) <- c("table", "chi2", "N", "doubleHet", "SNP")
    
    if(all.solutions)
    {
      output[[5]] <- allSol
      names(output)[length(output)] <- "Poo"
    }
    
    for(i in 1:length(LD))
    {
      output<-append(output, LDm[i])
      #names(output)[ length(output) ] <- LD[i]
    }
    names(output)[(length(output)- length(LD) +1 ) : length(output)] <- LD
    if(CI)									# interval estimates
    {
      if(returnLDdist)						# the interval estimates AND the distribution resulting from 'n.sim'
      {								# simulations are returned
        output <- append(output, ci.res[1])
        output <- append(output, ci.res[2])
        names(output)[(length(output) - 1):length(output)] <- c("LDdist", "CI" )
      }
      else								# only the interval estimates are returned
      {
        output <- append(output, ci.res[2])
        names(output)[length(output)] <- "CI"
      }
    }
    
    #########################################
    return(output)
  }


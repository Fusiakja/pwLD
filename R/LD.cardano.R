`LD.cardano` <-
function( snps, code=c(0,1,2,3), LD=c("Dprime", "Q", "r", "Y", "HS"), HSweight=4 , CI=T, strategy=c("jackknife", "bootstrap"), alpha=.05, n.sim=1000, returnLDdist=F, paradigm=c("freq", "bayes", "fullbayes"), dirich=rep(1,9) , all.solutions=F, tol=.Machine$double.eps^.6, digits=12, seed=F, mc=1000, intervall=c(0,1))
{
	paradigm <- match.arg(paradigm)
					
	if(class(snps) == "data.frame") snps <- as.matrix(snps)
	if(dim(snps)[1] < 2) stop("\n'snps' must be a 2 x N matrix!! \n\n")
	if(dim(snps)[1] > 2) stop("\nMore than two SNPs detected! Please refer to function 'LD.all'!\n\n")
	if(is.null(rownames(snps))) rownames(snps) <- c("snp1","snp2")
	#genotypes <- genotypes.3x3(snps, code=code, method=paradigm, dirich=dirich)
	

	# estimate genotypic frequencies	
	genotypes <- estGeno(snps, code=code, paradigm=paradigm, dirich=dirich)
	genotypes.freqs <- genotypes$count/genotypes$N
	genotypes.counts <- genotypes$count	
	# sample size without missings
	Neff <- genotypes$N 
	
	# estimate allele frequencies
	af<- estAF(genotypes.counts/genotypes.counts[4,4])
	
	# estimate pairwise haplotype frequencies
	res <- estHap( genotypes.counts/genotypes.counts[4,4], tol=tol, digits=digits, genotypes.counts[4,4])
	hapCounts <- res$table*2*Neff
	allSol <- res$solutions
	
	if(seed==TRUE)
	{
	  set.seed(1)
	}
	hapFreq <- estHapFreq(hapCounts = hapCounts, paradigm = paradigm, Neff=2*Neff, dirich=dirich, mc=mc, seed=seed)
  hapFreq <- hapFreq[[1]]
  
	#if any allele frequency is zero or one, all LD measures are set to 'NA'
	if(any(af == 0) | any(af == 1))
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

	# interval estimates
	if(CI)
	{
	  if(seed==TRUE)
	  {
	    set.seed(1)
	  }
	  ci.res <- estConfidenceInterval(hapCounts = genotypes.counts, hapFreq = genotypes.freqs, N=Neff, paradigm = paradigm,nSim= n.sim, LD=LD,  tol=tol , alpha=alpha, digits=digits, strategy=strategy, Dir=dirich, intervall=intervall, mc=mc)
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


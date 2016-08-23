`LD.cardano` <-
function( snps, code=c(0,1,2,3), LD=c("Dprime", "Q", "r2"),  CI=T, alpha=.05, n.sim=5000, returnLDdist=F, paradigm=c("freq", "bayes"), dirich=rep(1,9) , all.solutions=F, tol=.Machine$double.eps^.6, digits=12)
{
	paradigm <- match.arg(paradigm)
					
	if(class(snps) == "data.frame") snps <- as.matrix(snps)
	if(dim(snps)[1] < 2) stop("\n'snps' must be a 2 x N matrix!! \n\n")
	if(dim(snps)[1] > 2) stop("\nMore than two SNPs detected! Please refer to function 'LD.all'!\n\n")
	if(is.null(rownames(snps))) rownames(snps) <- c("snp1","snp2")
	#genotypes <- genotypes.3x3(snps, code=code, method=paradigm, dirich=dirich)
	

	# estimate genotypic frequencies	
	genotypes <- estGeno(snps, code=code, paradigm=paradigm, dirich=dirich)
	genotypes.freqs <- genotypes$freq
	genotypes.counts <- genotypes$count	
	# sample size without missings
	Neff <- genotypes$N 
	
	# estimate allele frequencies
	af<- estAF(genotypes.freqs)
	
	# estimate pairwise haplotype frequencies
	res <- estHap( genotypes.freqs, tol=tol, digits=digits)
	hapFreq <- res$table
	allSol <- res$solutions


	# if any allele frequency is zero or one, all LD measures are set to 'NA'
	if(any(af == 0) | any(af == 1))
	{
		warning("Allele frequencies of 0 or 1 !!")
		output <- vector("list", length(LD)+1)
		output[[1]] <- hapFreq
		for(i in 1:length(LD))
			output[[i+1]] <- NA
		names(output) <- c("table", LD)

		return(output)
	}

	
	# estimate LD measures
	LDm <- estLD(hapFreq, LD)
	
	# interval estimates
	if(CI)
	{
		if(paradigm == "freq")
			ci.res <- estConfidenceInterval(genotypes.freqs, N=Neff, nSim= n.sim, LD=LD,  tol=tol , alpha=alpha, digits=digits) 

		if(paradigm == "bayes")
			ci.res <- estCredibleInterval(genotypes.counts, nSim=n.sim, LD=LD, Dir=dirich, tol=tol, alpha=alpha, digits=digits)
	}


	#########################################
	# generate the output
	output <- list()
	output[[1]] <- hapFreq					# table of haplotype frequencies
	output[[2]] <- estLD(hapFreq*2*Neff, "chi2")			# chi square statisitcs, deviation from complete linkage equilibrium
	output[[3]] <- Neff						# sample size, without missings
	output[[4]] <- genotypes.counts[as.character(code[2]), as.character(code[2])]
	output[[5]] <- genotypes$SNP
	names(output) <- c("table", "chi2", "N", "doubleHet", "SNP")

	if(all.solutions) 
	{
		output[[5]] <- allSol
		names(output)[length(output)] <- "Poo"
	}

	for(i in LD)
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


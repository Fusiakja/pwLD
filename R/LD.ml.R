`LD.ml` <-
function(snps, code = c(0, 1, 2, 3), LD = c("Dprime"), paradigm = c("freq","bayes"), dirich = rep(1,9), tol = .Machine$double.eps, always.opt = T)
{
	paradigm <- match.arg(paradigm)


	
	if(class(snps) == "data.frame") snps <- as.matrix(snps)
	if(dim(snps)[1] < 2) stop("\n'snps' must be a 2 x N matrix!! \n\n")
	if(dim(snps)[1] > 2) stop("\nMore than two SNPs detected! Please refer to function 'LD.all'!\n\n")
	if(is.null(rownames(snps))) rownames(snps) <- c("snp1","snp2")
		
	
	# estimate genotypic frequencies	
	genotypes <- estGeno(snps, code=code, paradigm=paradigm, dirich=dirich)
	genotypes.freqs <- genotypes[[ 1 ]] 
	genotypes.counts <- genotypes[[ 2 ]]

	# sample size without missings
	N <- genotypes[[ 3 ]] 
	
	# estimate allele frequencies
	af<- estAF(genotypes.freqs)		
	po. <- af[1]
	p.o <- af[2]

	if(any(af == 0) | any(af == 1)) stop("Allele frequencies 0 or 1!!")

	if(paradigm=="freq")
	{
		
		k <- genotypes.counts[as.character(code[2]), as.character(code[2])]

		if((k == 0) & (!always.opt))
		{
			 Poo <- (genotypes.freqs[as.character(code[1]),as.character(code[1])] + 0.5*genotypes.freqs[as.character(code[1]),as.character(code[2])] + 0.5*genotypes.freqs[as.character(code[2]),as.character(code[1])])
		
			likelihood <- NA
			Lmax <- NA
		}
	
		else
		{	# determine the likelihood

			likelihood <- determine.likelihood(genotypes.counts, af)$L
			likelihood.fct <- function(Poo) eval(parse(text=likelihood))


			# valid range of Poo
			Poo.low <- max(0, sum(af) - 1)
			Poo.up <- min(af)

			# OPTIMIZATION
			if(Poo.low != Poo.up)
			{	
				res <- optimise(likelihood.fct, c(Poo.low,Poo.up), maximum=T, tol=tol)
				Poo <- res$maximum
				Lmax <- res$objective
			}
			else 
			{	
				Poo=Poo.low
				Lmax <- likelihood.fct(Poo)
			}		
		}
	}

	#
	if(paradigm == "bayes")
	{	

		# number of ambiguous cases
		k <- genotypes.freqs[as.character(code[2]), as.character(code[2])]*N
	
		# determine the likelihood
		likelihood <- determine.likelihood(genotypes.freqs*N, af)$L
		likelihood.fct <- function(Poo) eval(parse(text=likelihood))


		# valid range of Poo
		Poo.low <- max(0, sum(af) - 1)
		Poo.up <- min(af)

		# OPTIMIZATION
		if(Poo.low != Poo.up)
		{	
			res <- optimise(likelihood.fct, c(Poo.low,Poo.up), maximum=T, tol=tol)
			Poo <- res$maximum
			Lmax <- res$objective
		}
		else 
		{	
			Poo=Poo.low
			Lmax <- likelihood.fct(Poo)
		}	

	}	
	###################################################
	# OUTPUT
	output <- list()
	output[1] <- Poo
	for(i in LD) output<-append(output, linkage(Poo, af[1], af[2], method=i))
	names(output) <- c("Poo", LD)

 return(output)
 
}


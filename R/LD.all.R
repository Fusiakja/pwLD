`LD.all` <-
function(data,  code=c(0,1,2,3), LD=c("D", "Dprime", "Q", "r2", "OR", "MI", "chi2"), MAF=0, paradigm=c("freq", "bayes"), dirich=rep(1,9),  verbose=T, tol=.Machine$double.eps^.6, digits=12)
{
	LD <- match.arg(LD, several.ok=T)

	time <- proc.time()[3]

	paradigm <- match.arg(paradigm)
	
	# check if the data is in the right format
	data <- validate.DS(data, code=code)	

	if(! ((paradigm == "bayes" ) && (MAF == 0)) )
	{
		# remove SNPs with MAF beolw the threshold 
		allele.freqs <- apply(data, 1, allele.freq, code)
    		data <- data[ which(allele.freqs > MAF & allele.freqs < (1-MAF) ), ]
	}
	if(verbose) cat("\n",dim(data)[1],"SNPs with MAF >", MAF, "in ",dim(data)[2], "samples left.\n\n Estimating LD...")
	
	# number of SNPs
	n <- dim(data)[1]

	# number of entries of the triangular LD matrix
	Nentries <- (n*(n-1)/2)

	# the externat C code
	# the output is stored as a vector representing the triangular matrix of the LD matrix
	res <- .C("LDall", data=as.character(as.vector(data)), nrow=as.integer(dim(data)[1]), ncol=as.integer(dim(data)[2]), LD=as.character(LD), LDnumb=as.integer(length(LD)), code=as.character(code), paradigm=as.character(paradigm), Dir=as.double(dirich), MAF=as.double(MAF), tol=as.double(tol), as.integer(digits), LDmatPtr=as.double(rep(-5, length(LD)*(n*(n-1)/2))), PACKAGE="pwLD"  )

	
	output <- vector("list", length(LD))
	names(output) <- LD
	
	# generate the output
	LDnumb = 1
	for(i in 1:length(LD))
	{	# convert the triangular matrix into a symmetric LD matrix
		output[[i]] <- vec2sm( res$LDmatPtr[ ((LDnumb-1)*Nentries+1)  : (((LDnumb-1)*Nentries+1) + (n*(n-1)/2) - 1 ) ] )
		dimnames(output[[i]]) <- list(rownames(data), rownames(data))
		LDnumb = LDnumb+1
	}


	if(verbose) cat(paste("Done\n\n computed in ",round(proc.time()[3]-time, 1)," seconds (",round(((proc.time()[3]-time)/60),1)," minutes).\n\n", sep=""))

	return(output)
}


`LD.all` <-
function(data,  code=c(0,1,2,3), LD=c("D", "Dprime", "Q", "r", "OR", "MI", "chi2", "Y", "HS"), MAF=0, paradigm=c("freq", "bayes", "fullbayes"), strategy = c("bootstrap", "jackknife","zapata"),dirich=rep(1,4),  verbose=T, tol=.Machine$double.eps^.6, digits=12, CI=F, HSweight=4, alpha=0.2, nSim=1000, seed=F, inervall=c(0,1))
{
	LD <- match.arg(LD, several.ok=T)

	time <- proc.time()[3]

	paradigm <- match.arg(paradigm)
  
	sure="y"
	if(CI==T&&strategy=="bootstrap")
	{
	  #cat ("This may take a while. Are you sure? Press 'y'.")
	  #sure <- readline()
	}
	if(!((sure == "Y") || (sure=="y")) )
	{
	 #  return()
	}
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
	ifci <- 0
	
	if(CI==TRUE)
	{ifci <- 1}

	if(seed==TRUE)
	{
	  set.seed(1)
	}
	res <- .C("LDall", data=as.character(as.vector(data)), nrow=as.integer(dim(data)[1]), ncol=as.integer(dim(data)[2]), LD=as.character(LD), LDnumb=as.integer(length(LD)), code=as.character(code), paradigm=as.character(paradigm), Dir=as.double(dirich), MAF=as.double(MAF), tol=as.double(tol), digits=as.integer(digits), LDmatPtr=as.double(rep(-5, length(LD)*(n*(n-1)/2))), HSweight=as.double(HSweight), ci=as.integer(ifci), mc=as.integer(1000), strategy=as.character(strategy), alpha=as.double(alpha), cilow=as.double(rep(-5, length(LD)*(n*(n-1)/2))), ciup=as.double(rep(-5, length(LD)*(n*(n-1)/2))), nSim=as.integer(nSim), LDdist=as.double(rep(-5, nSim)), PACKAGE="pwLD"  )

	if(CI==TRUE)
	{
	output <- vector("list", length(LD))
	names(output) <- LD
	outputlow <- vector("list", length(LD))
	names(outputlow) <- LD
	outputup <- vector("list", length(LD))
	names(outputup) <- LD
	
	
	LDnumb = 1
	for(i in 1:length(LD))
	{	# convert the triangular matrix into a symmetric LD matrix
	  output[[i]] <- vec2sm( res$LDmatPtr[ ((LDnumb-1)*Nentries+1)  : (((LDnumb-1)*Nentries+1) + (n*(n-1)/2) - 1 ) ] )
	  dimnames(output[[i]]) <- list(rownames(data), rownames(data))
	  LDnumb = LDnumb+1
	}
	LDnumb = 1
	for(i in 1:length(LD))
	{	# convert the triangular matrix into a symmetric LD matrix
	  outputlow[[i]] <- vec2sm( res$cilow[ ((LDnumb-1)*Nentries+1)  : (((LDnumb-1)*Nentries+1) + (n*(n-1)/2) - 1 ) ] )
	  dimnames(outputlow[[i]]) <- list(rownames(data), rownames(data))
	  LDnumb = LDnumb+1
	}
	LDnumb = 1
	for(i in 1:length(LD))
	{	# convert the triangular matrix into a symmetric LD matrix
	  outputup[[i]] <- vec2sm( res$ciup[ ((LDnumb-1)*Nentries+1)  : (((LDnumb-1)*Nentries+1) + (n*(n-1)/2) - 1 ) ] )
	  dimnames(outputup[[i]]) <- list(rownames(data), rownames(data))
	  LDnumb = LDnumb+1
	}
	
	output <- c(output, outputlow, outputup)
	
	}
	else
	{
	  output <- vector("list", length(LD))
	  names(output) <- LD
	  
	  #generate the output
	  LDnumb = 1
	  for(i in 1:length(LD))
	  {	# convert the triangular matrix into a symmetric LD matrix
	    output[[i]] <- vec2sm( res$LDmatPtr[ ((LDnumb-1)*Nentries+1)  : (((LDnumb-1)*Nentries+1) + (n*(n-1)/2) - 1 ) ] )
	    dimnames(output[[i]]) <- list(rownames(data), rownames(data))
	    LDnumb = LDnumb+1
	  }
	  output <- c(output)
	  
	}
	
	if(verbose) cat(paste("Done\n\n computed in ",round(proc.time()[3]-time, 1)," seconds (",round(((proc.time()[3]-time)/60),1)," minutes).\n\n", sep=""))

	return(output)
}


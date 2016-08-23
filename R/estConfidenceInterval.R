`estConfidenceInterval` <-
function(genoFreq, N, nSim= 1000, alpha=.05, LD="Q", tol=.Machine$double.eps^.6, digits=12)
{	
	ci <- matrix(0, nrow=length(LD), ncol=2)
	rownames(ci) <- LD
	colnames(ci) <- c( paste(alpha/2*100, "%", sep=""), paste( (1- (alpha/2))*100, "%", sep="" ) )

	res <- .C("confidenceInterval", genoFreq=as.double(t(genoFreq)), N=as.integer(N), nSim=as.integer(nSim), LD=as.character(LD), LDnumb=as.integer(length(LD)), tol=as.double(tol), digits=as.integer(digits),LDdist=as.double(numeric(length(LD) * nSim) ), PACKAGE="pwLD")

	# the distribution of the LD measures 
	LDdist <- matrix(res$LDdist, nrow=length(LD), byrow=T)
	rownames(LDdist) <- LD

	# confidence intervals
	for(i in 1:length(LD))
		ci[i, ] <- quantile(LDdist[i, ], c(alpha/2, 1- (alpha/2)) )

	output <- vector("list", 2)
	output[[1]] <- LDdist
	output[[2]] <- ci

	return(output)	
}


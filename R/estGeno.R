`estGeno` <-
function(snps, code, paradigm="freq", dirich=rep(1,9), mc=1000)
{

	res=.C("estimateGenoFreq", genoA=as.character(snps[1,]), genoB=as.character(snps[2, ]), code=as.character(code), N=as.integer(length(snps[1, ])), paradigm=as.character(paradigm), Dir=as.double(dirich), genoFreq=matrix(0, nrow=4, ncol=4, byrow=T), genoCounts=matrix(0, nrow=4, ncol=4, byrow=T), Neff=as.integer(0), PACKAGE="pwLD" ) 

	#
	freqs <- t(res$genoFreq)
	counts <- t(res$genoCounts)
	dimnames(freqs) <- dimnames(counts) <- list(  c(code[1:3], "Sum"), c(code[1:3], "Sum"))

	out <- vector("list", 4)
	out[[1]] <- freqs
	out[[2]] <- counts
	out[[3]] <- res$Neff
	out[[4]] <- rownames(snps)
	names(out) <- c("freq", "count", "N", "SNP")	

	return(out)
}


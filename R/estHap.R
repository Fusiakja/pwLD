`estHap` <-
function( genoFreq, tol, digits=12)
{
	res=.C("estimateHaploFreq", genoFreq=as.double(t(genoFreq)), haploFreq=as.double(matrix(0, nrow=3, ncol=3)), tol=as.double(tol), digits=as.integer(digits), pooHat=as.double(rep(0,3)), PACKAGE="pwLD" )


	output <- list()
	output[[1]] <- matrix(res$haploFreq, nrow=3, byrow=T, dimnames=list(c("0", "1", "Sum"), c("0", "1", "Sum")))
	output[[2]] <- res$pooHat

	names(output) <- c("table", "solutions")
	return( output )
}


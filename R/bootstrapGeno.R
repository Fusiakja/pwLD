`bootstrapGeno` <-
function(genoFreq, N)
{
	res <- .C("bootstrapGenoFreq", genoFreq=as.double(t(genoFreq)), N=as.integer(N), genoFreqBS=matrix(0.0, nrow=4, ncol=4), PACKAGE="pwLD" )

	return(t(res$genoFreqBS))
}


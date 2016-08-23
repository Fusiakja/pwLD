`estAF` <-
function(genoFreq)
{
	res <- .C("alleleFreq", genoFreq=as.double(t(genoFreq)), af=as.double(c(0,0)), PACKAGE="pwLD")

	return(res$af)
}


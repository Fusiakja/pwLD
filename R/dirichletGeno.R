`dirichletGeno` <-
function(genoCounts, Dir)
{
	res <- .C("sampleDirichletGenoFreq", genoCounts=as.double(t(genoCounts)), Dir=as.double(Dir), genoFreqRS=matrix(0.0, nrow=4, ncol=4), PACKAGE="pwLD" )

	return(res$genoFreqRS)
}


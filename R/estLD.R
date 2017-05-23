`estLD` <-
function(haploFreq, LD=c("D", "Dprime", "r", "Q", "OR", "MI","Y","HS", "chi2" ), HSweight=4)
{
	LD = match.arg(LD, several.ok=T)

	res <- .C("estimateLD", tab=as.double(t(haploFreq)), what=as.character(LD), numb=as.integer(length(LD)), LD=as.double(rep(0, length(LD))), HSweight=as.double(HSweight), PACKAGE="pwLD"  )

	ld <- res$LD
	#names(ld) <- LD

	return( ld )
}


`estLD` <-
function(haploFreq, LD=c("D", "Dprime", "r2", "Q", "OR", "MI", "chi2" ))
{
	LD = match.arg(LD, several.ok=T)

	res <- .C("estimateLD", tab=as.double(t(haploFreq)), what=as.character(LD), numb=as.integer(length(LD)), LD=as.double(rep(0, length(LD))), PACKAGE="pwLD"  )

	ld <- res$LD
	names(ld) <- LD

	return( ld )
}


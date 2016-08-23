`linkage` <-
function(Poo,po.,p.o, method=c("Dprime"))
{
	method <- match.arg(method)

        # in order to prevent the function of chrashing, 'NA' is returned if the allel frequencies are zero or one
        if(any(c(po.,p.o) == 0) | any(c(po.,p.o) == 1)) return(NA)

	# conbtingeny table
	tab <- snp.contigency.table(Poo, po.,p.o)

 	P01 <- tab[1,2]
	P10 <- tab[2,1]
	P11 <- tab[2,2]
	
	# deviation from independence 
	D <- Poo - po.*p.o

	if(method=="Dprime")
	{
		if(D >= 0 ) D.max <- min( po.*(1-p.o) ,  (1-po.)*p.o  )	
		else D.max <- min( po.*p.o, (1-po.)*(1-p.o) )
		return(D/D.max)
	}
	#else stop("no valid method chosen!")
}


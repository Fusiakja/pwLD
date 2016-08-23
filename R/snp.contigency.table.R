`snp.contigency.table` <-
function(Poo, po., p.o)
{
	
	cont.table <- matrix(NA, 2, 2)
	dimnames(cont.table) <- list(c("0", "1"), c("0", "1"))

	###################################################################
	# - determine the entries of the contingency table
	# make sure that all entries are positive
	P01 <- max(0, po. - Poo)
	P10 <- max(0, p.o - Poo)
	P11 <- max(0, 1 - Poo - P01 - P10)

	# rescale the entries such that the p_ij sum to one
	# NOT NECESSARY since this is ensured by the calculation of Pij above
	
	sum.Pij <- (Poo + P01 + P10 + P11)

	cont.table[1,1]  <- Poo 
 	cont.table[1,2] <- P01 / sum.Pij
	cont.table[2,1] <- P10 / sum.Pij
	cont.table[2,2] <- P11 / sum.Pij

	return(cont.table)	

}


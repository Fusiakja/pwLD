`determine.likelihood` <-
function(genotypes.counts, marginals)
{
	# marginal distributions
	po. <- marginals[1]
	p.o <- marginals[2]
	
	code <- rownames(genotypes.counts)[1:3]

	###################################
	# count the number of different cases
	cases <- list()
	names.cases <- list()

	for(i in code)
	{
		for(j in code)
		{
			cases <- append(cases, genotypes.counts[i,j])
			names.cases <- append(names.cases, paste(i,j,sep="-"))
		}
	}
	names(cases) <- names.cases
	cases <- unlist(cases)

	# determine the log likelihood...
	likelihood <- "0"
	
	cases <- cases[which(cases > 0)]

	for(i in 1 : length(cases))
		likelihood <- paste(likelihood,"+", log.likelihood(names(cases)[i], code=code),"*", cases[i],sep="" )
		
	####################
	# output
	output <- vector("list",3)
	names(output) <- c("L", "po.", "p.o")
	output[1]<-likelihood
	output[2]<-po.
	output[3]<-p.o

	return(output)

}


`validate.DS` <-
function(data, code)
{
	# if 'data' is not a matrix
	if(is.null(dim(data))) stop("'data' is not a matrix!")
       
	# if 'data' has no unique rownames, generate ones
 	if(is.null(rownames(data)))
	{  # stop("No rownames of  matrix 'data' specified!\nPlease add unique identifiers for each SNP!\n  ")
		rownames(data) <- 1:dim(data)[1]
 	}
	if(mode(code) == "numeric")
	{
		# if 'data' is not of mode 'numeric'
		if(mode(data)!="numeric")
 		{
			# this conversation may take a lot of time...
			if(mode(data) == "list") data <- matrix(unlist(data), ncol=dim(data)[2], byrow=F,dimnames=dimnames(data))
		
			else mode(data) <- "numeric"
	 	}
	}

	 if(length(levels(as.factor(data))) > length(code) | sum(levels(as.factor(data)) %in% code) < length(levels(as.factor(data))))  stop("code used in 'data' not equal to parameter 'code'!")

return(data)
}


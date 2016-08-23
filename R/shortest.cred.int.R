`shortest.cred.int` <-
function(dist, alpha)
{
	# ordering
	dist.ord <- dist[order(dist)]

	int.length <- c()
	quant <- list()
	i <- 1
	for(q in seq(0,alpha, 1e-3))
	{
		quant[[i]] <- quantile(dist, c(q, (1-alpha)+q ))
	
		int.length[i] <- max(dist.ord[dist.ord <= quant[[i]][2]]) - min(dist.ord[dist.ord >= quant[[i]][1]])
	
		i <- i+1
	}

	# determine shortest interval
	return(quant[[which.min(int.length)]])

}


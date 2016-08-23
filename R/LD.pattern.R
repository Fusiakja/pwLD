`LD.pattern` <-
function( LD.mat , map=NULL, window=1.7e6,  minSNPs = 10, scale=c("Mb", "Kb", "bp"), plot=T)
  {
	unit <- match.arg(scale)
	if(unit == "Mb") scale <- 1e6
	else if(unit == "Kb") scale <- 1e3
	else scale <- 1

       if(is.null(map)) stop("No genomic position of SNPs specified!\n")

	# determine the LD measures used
	LD <- names(LD.mat)[1:(length(LD.mat))]

	# ids of all snps
        snp.ids <- rownames(LD.mat[[LD[1]]])

	# absolute genome positions in base pairs
	 BP <- map
	 names(BP) <- snp.ids
	
	BP <- BP[snp.ids]
 
        # number of LD measures to plot
        ld.numb <- length(LD)

        # list to store the mean LD values for each LD measure
        slided <- vector("list", ld.numb)
	for(ld in LD)
		slided[[ld]] <- c() #vector("numeric", (dim(data)[1] - window))
	names(slided) <- LD

	# vector to store the mean position in bp
	LD.BP <- c()
	numbSNP <- c()
	
	
	for(i in 1:(length(snp.ids)-minSNPs+1))
	{
		# deteremine all SNPs covered by window
		BP.start <- BP[ snp.ids[i] ]
		BP.end <- BP.start + window
		
		SNPs <- snp.ids[ which(BP >= BP.start & BP <= BP.end) ]

		

		if(length(SNPs) >= minSNPs) 
		{ 	 
			for(ld in LD)
			{	# extract the part of the matrix that contains 'SNPs'
				mat <- LD.mat[[ld]][ SNPs, SNPs]
				mat.tri <- mat[upper.tri(mat)]
				# mean LD value
				slided[[ld]][i]  <- mean(abs(mat.tri), na.rm=T)
			}
			# mean genomic position
			LD.BP[i] <- mean(BP[rownames(mat)], na.rm=T)

			# number of SNPs covered by 'window'
			numbSNP[i] <- length(SNPs)

		}	
	
        }

	# if 'window' is chosen to small such that no SNPs are covered, the function returns an error 
	if(sum(is.null(LD.BP))) stop("\nNo SNPs covered by 'window'! Increase window size or decrease parameter 'minSNPs'! \n\n")

	# if only some windows cover no SNPs, a warning is returned
	if(is.na(sum(LD.BP))) warning("\nSome windows contained no SNPs! \n\n")

	# make the plot 
	if(plot)
	{
		plot(LD.BP/scale, slided[[1]], ylim=c(0,1), xlab=paste("genome position (", unit," )"), ylab="LD", type="l", yaxp=c(0,1,10))
		if(ld.numb > 1)
       			for(i in 2: ld.numb)
              			lines(LD.BP/scale, slided[[i]], col=i, type="l")
        	legend("topright", legend=LD, fill = c(1:length(LD)))
		legend("top", legend=c(paste("#SNPs", length(snp.ids)), paste("window size:", (window/scale), unit) ) )
	}
	
       #######################################################
       # some informations
       info <- matrix(ncol=1, nrow=3, dimnames=list(c("# SNPs", "average distance between 2 SNPs", "window"), " "))
       info[1,1] <- length(snp.ids)
       info[2,1] <- paste(round((range(BP)[2]-range(BP)[1])/length(snp.ids)/scale, 2), unit)
       info[3,1] <- window


 	output <- vector("list", 4)
	names(output)<-c("LD.mean", "BP.mean", "#SNPs","info")


	output[[1]] <- slided
	output[[2]] <- LD.BP
	output[[3]] <- numbSNP
	output[[4]] <- info

 return(output)
}


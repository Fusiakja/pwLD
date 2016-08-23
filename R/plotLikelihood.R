`plotLikelihood` <-
function(snps,  code=c(0,1,2,3),   paradigm=c("freq", "bayes"),  dirich=rep(1,9), main="",  ylim=NULL,   legend=F, tol=.Machine$double.eps^.6, plot.tol=.Machine$double.eps^.25, digits=12)
{
	if(class(snps) == "data.frame") snps <- as.matrix(snps)
	if(is.null(dim(snps))) stop("\n'snps' must a 2 x N matrix! \n\n")
	if(dim(snps)[1] > 2) stop("\nMore than two SNPs detected!\n\n")

	paradigm <- match.arg(paradigm)

	# estimate genotypic frequencies	
	genotypes <- estGeno(snps, code=code, paradigm=paradigm, dirich=dirich)
	genotypes.freqs <- genotypes[[ 1 ]] 
	genotypes.counts <- genotypes[[ 2 ]]

	# sample size without missings
	N <- genotypes[[ 3 ]] 
	
	# estimate allele frequencies
	af<- estAF(genotypes.freqs)
	po. <- af[1]; p.o = af[2];

	if(any(af == 0) | any(af == 1)) stop("Allele frequencies of 0 or 1 !!")
	
	# determine the likelihood function 
        L <- function(Poo) eval(parse(text=determine.likelihood(genotypes.freqs*N, af )$L) )

	# MLE
	res.L <- LD.ml(snps, code=code,  paradigm=paradigm, dirich=dirich)

	# analytical fix-point solution, all solutions of the cubic polonomial are returned
	s <- LD.cardano(snps, all.solutions=T, code=code, paradigm=paradigm, dirich=dirich, tol=tol, digits=digits, CI=F)$Poo

	#####################################################
	# eliminate solutions that are not in the valid interval
	# vaild range of Poo
	low <- max(0, sum(af)-1)
 	up <- min(af)
		
	idx.low <- abs(low - s) < tol
  	idx.up <- abs(s - up) < tol

  	if(sum(idx.low) > 0) s[idx.low] <- low
  	if(sum(idx.up) > 0) s[idx.up] <- up

  	idx <- (s >= low) & (s <= up)
  	if(sum(idx) == 0) stop(paste("no solution in valid range\n"))

  	s <- s[idx]
  	#####################################################
	# valid interval	
	x <- seq(max(0+plot.tol, sum(af)-1+plot.tol), min(af)-plot.tol, (min(af)-max(0,sum(af)-1) )/1000 )
	# log-likelihood	
	L.x <- L(x)

	L.x[L.x == Inf] <- NA
	L.x[L.x == -Inf] <- NA

	#####################################################
	# scale of y-axis
	mi <- min(L.x, na.rm=T)
	ma <-max(L.x, na.rm=T)
	r <- (ma-mi)*0.25

	# plot the likelihood
	if(is.null(ylim)) plot(x,L.x, type="l", ylim=c(mi, ma+r),xlab=substitute(paste("valid range of p"[0*0], sep="")), ylab="ln( likelihood )", main=main, sub=paste(rownames(snps)[1], "vs.", rownames(snps)[2]))
	else plot(x,L.x, type="l", xlab=substitute(paste("valid range of p"[0*0], sep="")), ylab="ln( likelihood )", main=main, ylim=ylim, sub=paste(rownames(snps)[1], "vs.",rownames(snps)[2]))

	#####################################################
	# LEGEND: LD.cardano, ONE solution
	if(length(s) < 2) 
	{	poo <- round(s,3)
		L.poo <- round(L(s),3)
		
		text.ANA <- expression(1,1)
		text.ANA[[1]] <- substitute(paste(widehat(p)[0*0]," = ",poo, ", L(",widehat(p)[0*0],") = ", L.poo ,sep=""), list(poo=eval(poo), L.poo=eval(L.poo)))
		text.ANA[[2]] <- substitute(paste(widehat(D),"(",widehat(p)[0*0],") = ",Dprime, sep=""), list(Dprime=eval(linkage(s, po., p.o, method="Dprime"))))

		abline(v=s, col="red")
		abline(v=res.L$Poo, col="blue")
		
		if(legend)
			legend("topright", legend=text.ANA, title="function 'LD.cardano'", bg="white",col="red", lty=c(1,-1), y.intersp=1.5)
	}

	#####################################################
	# LEGEND: LD.cardano, non-unique solution
	else
	{ 	COL <- c("red", "darkred", "orange")
		poo <- round(s,3)
		L.poo <- round(L(s),3)

		if(length(s)==3) text.ANA <- expression(1,1,1,1)
		else text.ANA <- expression(1,1,1)

		for(i in 1:length(s))
			text.ANA[[i]]	<- substitute(paste(widehat(p)[0*0][i]," = ",poo,", L(",widehat(p)[0*0][i],") = ", L.poo, sep=""), list(i=eval(i),poo=eval(poo[i]), i= eval(i), L.poo=eval(L.poo[i])))	

		text.ANA[[length(s)+1]] <- substitute(paste(widehat(D),"'(",widehat(p)[0*0][Poomax],") = ",Dprime, sep=""), list(Poomax=eval(which.max(L(s))), Dprime=eval(linkage(s[which.max(L(s))], po., p.o, method="Dprime"))))
		
		abline(v=s, col=COL)
		
		# result ML
		abline(v=res.L$Poo, col="blue")

		if(legend)
			legend("topright", legend=text.ANA, title="function 'LD.cardano'", bg="white",col=COL[1:length(s)], lty=c(rep(1,length(s)),-1),y.intersp=1.5)

	}


	####################################################
	# LEGEND: LD.ml
	ML.poo <- round(res.L$Poo, 3)
	ML.L.poo <- round(L(res.L$Poo),3)
	
	text.ML <- expression(1, 1)
	text.ML[[1]] <-substitute(paste (widehat(p)[0*0]," = ", ML.poo, ", L(",widehat(p)[0*0],") = ",ML.L.poo, sep=""), list(ML.poo=eval(ML.poo), ML.L.poo = eval(ML.L.poo)))
	text.ML[[2]] <- substitute(paste(widehat(D),"'(",widehat(p)[0*0],") = ",Dprime, sep=""), list(Dprime=eval(linkage(res.L$Poo, po., p.o,method="Dprime"))))

	if(legend)
		legend("topleft", legend= text.ML,title="numerical likelihood optimisation", bg="white", col="blue", lty=c(1, -1, -1), y.intersp=1.5)


	#####################################################	
	# LEGEND:  data properties
	text.data <- expression(1,1,1)

	poo.int <- c(round(max(po.+p.o-1, 0),4), round(min(po.,p.o),4))
	text.data[[1]] <- substitute( paste("valid interval:  ",widehat(p)[0*0] %in%" [",pmin,",  ",pmax,"]"), list(pmin=eval(poo.int[1]), pmax=eval(poo.int[2])) )
	text.data[[2]] <-  paste("N = ", dim(snps)[2], sep="")
	text.data[[3]] <- paste("# double heterozygotes: ",round(genotypes.counts[as.character(code[2]), as.character(code[2])],3), sep="")
	
	if(legend)
		legend("bottom", legend=text.data, bg="white", title="data properties" ,y.intersp=1.5)
	
}


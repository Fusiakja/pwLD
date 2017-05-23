`LD.HaploBlocks` <-
function(data,  code=c(0,1,2,3), LD=c("D", "Dprime", "Q", "r2", "OR", "MI", "chi2"), MAF=0, paradigm=c("freq", "bayes"), strategy = c("bootstrap", "zapata"),dirich=rep(1,9),  verbose=T, tol=.Machine$double.eps^.6, digits=12, CI=F, HSweight=4, alpha=0.2, nSim=1000, seed=F, map=NULL)
{
  LD <- match.arg(LD, several.ok=T)
  
  time <- proc.time()[3]
  
  paradigm <- match.arg(paradigm)
  
  # check if the data is in the right format
  data <- validate.DS(data, code=code)	
  
  if(! ((paradigm == "bayes" ) && (MAF == 0)) )
  {
    # remove SNPs with MAF beolw the threshold 
    allele.freqs <- apply(data, 1, allele.freq, code)
    data <- data[ which(allele.freqs > MAF & allele.freqs < (1-MAF) ), ]
  }
  if(verbose) cat("\n",dim(data)[1],"SNPs with MAF >", MAF, "in ",dim(data)[2], "samples left.\n\n Estimating LD...")
  
  # number of SNPs
  n <- dim(data)[1]
  
  # number of entries of the triangular LD matrix
  Nentries <- (n*(n-1)/2)
  
  # the externat C code
  # the output is stored as a vector representing the triangular matrix of the LD matrix
  res <- .C("MIG", data=as.character(as.vector(data)), nrow=as.integer(dim(data)[1]), ncol=as.integer(dim(data)[2]), LD=as.character(LD), LDnumb=as.integer(length(LD)), code=as.character(code), paradigm=as.character(paradigm), Dir=as.double(dirich), MAF=as.double(MAF), tol=as.double(tol), digits=as.integer(digits), LDmatPtr=as.double(rep(-5, length(LD)*(n*(n-1)/2))), HSweight=as.double(HSweight), ci=as.double(0), mc=as.integer(1000), strategy=as.character(strategy), alpha=as.double(alpha), cilow=as.double(rep(-5, length(LD)*(n*(n-1)/2))), ciup=as.double(rep(-5, length(LD)*(n*(n-1)/2))), nSim=as.integer(nSim), LDdist=as.double(rep(-5, nSim)), MIG1=as.double(rep(-5, length(LD)*(n*(n-1)/2))), MIG2=as.double(rep(-5, length(LD)*(n*(n-1)/2))),PACKAGE="pwLD"  )
  
  MIG_res <- list()
  for (li in 1:length(res$MIG1)) 
  {
    MIG_res[[li]] <- rownames(data)[as.integer(res$MIG1[li]): as.integer(res$MIG2[li])]  
  }
  
  MIG_map <- list()
  for (as in 1:n) 
  {
    MIG_map[[as]] <- numeric()
    for (bs in 1:length(MIG_res[[as]])) 
    {
      MIG_map[[as]] <-  c(MIG_map[[as]], map[MIG_res[[as]][bs], "pos"])
    } 
  }

  
  plot_line <- 0
  plot_line_max <- 0
  ploT_lines <- numeric()
  for (i in 1:n) 
  {
    if (!is.na(MIG_map[[i]][1])) 
    {
      if (i>1 && i<= length(MIG_map)) {
        if (length(intersect(MIG_map[[i]],MIG_map[[i-1]]))>0) 
        {
          plot_line <- plot_line+1
          ploT_lines <- c(ploT_lines, plot_line)
          if (plot_line_max < plot_line) {
            plot_line_max <- plot_line
          }
        }
        else
        {
          plot_line <- 0
          ploT_lines <- c(ploT_lines, plot_line)
        }
      }
    }
  }

  #par(new=T)
  #dev.off()
  
  #plot.new()
  for (i in 1:n) {
    if (i==1) {
      plot(MIG_map[[i]], rep(ploT_lines[i],length(MIG_map[[i]])), ylab = "", ylim = c(1,plot_line_max),type = "l", xlim = c(as.integer(map[rownames(data)[1], "pos"]),as.integer(map[rownames(data)[length(rownames(data))], "pos"])), xlab = "genome position", yaxt="n") 
      #mtext("Genome position",1, line = 3)
      }
    if (!is.na(MIG_map[[i]][1])){
      if (i>1 && i<= length(MIG_map)) {
        lines(MIG_map[[i]], rep(ploT_lines[i],length(MIG_map[[i]]))) 
      }
    }
  }
  if(verbose) cat(paste("Done\n\n computed in ",round(proc.time()[3]-time, 1)," seconds (",round(((proc.time()[3]-time)/60),1)," minutes).\n\n", sep=""))
  
  return(MIG_map)
}
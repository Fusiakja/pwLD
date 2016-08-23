`allele.freq` <-
function(SNP, code)
  {
    SNP <- SNP[ which(SNP != code[4]) ]

    N <- length(SNP)	

    freq <-  (2*sum(SNP == code[1]) + sum(SNP == code[2]) ) /(2*N)
    
    return( freq )
  }


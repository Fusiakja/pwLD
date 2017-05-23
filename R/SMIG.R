`SMIG` <-
  function(datas, HSweight=4)
  {
    cis <- numeric()
    for (i in 1:dim(datas)[1]) {
      
      LDs <- LD.cardano(datas[c(i,i+1),], paradigm = "freq", CI=T, LD="Dprime", strategy="zapata")
      cis <- rbind(cis, LDs$CI[1:2])
      
    }
    
    
    return( cis )
  }

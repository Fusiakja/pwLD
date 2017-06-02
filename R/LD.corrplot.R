`LD.corrplot` <- function(valuematrix, values =c("values", "upperCI", "lowerCI"), LD=c("D", "Dprime", "Q", "r", "OR", "MI", "chi2", "Y", "HS"), snp.zoom=FALSE, snp.target, snp.region, corr.labels=FALSE)
{
  
  if(snp.zoom==TRUE)
  {
    snp <- match(snp.target, colnames(valuematrix[[1]]))
    
    if((snp - snp.region)<1)
    {
      min.snp <- 1
    }
    else
    {
      min.snp <- (snp - snp.region)
    }
    
    
    if((snp + snp.region)>ncol(valuematrix[[1]]))
    {
      max.snp <- ncol(valuematrix[[1]])
    }
    else
    {
      max.snp <- (snp + snp.region)
    }
    
    if(corr.labels==TRUE)
    {
      if(values == "values")
      {
        region <- valuematrix[[1]][min.snp:max.snp, min.snp:max.snp]
        corrplot(region, method = "color", type = "upper", na.label = ".",tl.cex = 0.5)
      }
      else if(values == "lowerCI")
      {
        region <- valuematrix[[2]][min.snp:max.snp, min.snp:max.snp]
        corrplot(region, method = "color", type = "upper",na.label = ".",tl.cex = 0.5)
      }
      else if(values == "upperCI")
      {
        region <- valuematrix[[3]][min.snp:max.snp, min.snp:max.snp]
        corrplot(region, method = "color", type = "upper",na.label = ".",tl.cex = 0.5)
      }
      else
      {
        stop("no correlationmatrix was found", call. = TRUE)
      }
    }
    else
    {
      if(values == "values")
      {
        region <- valuematrix[[1]][min.snp:max.snp, min.snp:max.snp]
        corrplot(region, method = "color", type = "upper", na.label = ".",tl.cex = 0.5, tl.pos="n")
      }
      else if(values == "lowerCI")
      {
        region <- valuematrix[[2]][min.snp:max.snp, min.snp:max.snp]
        corrplot(region, method = "color", type = "upper",na.label = ".",tl.cex = 0.5, tl.pos="n")
      }
      else if(values == "upperCI")
      {
        region <- valuematrix[[3]][min.snp:max.snp, min.snp:max.snp]
        corrplot(region, method = "color", type = "upper",na.label = ".",tl.cex = 0.5, tl.pos="n")
      }
    }
  }
  else
  {
    if(corr.labels==TRUE)
    {
      if(values == "values")
      {
        region <- valuematrix[[1]]
        corrplot(region, method = "color", type = "upper",na.label = ".",tl.cex = 0.5)
      }
      else if(values == "upperCI")
      {
        region <- valuematrix[[2]]
        corrplot(region, method = "color", type = "upper",na.label = ".",tl.cex = 0.5)
      }
      else if(values == "lowerCI")
      {
        region <- valuematrix[[3]]
        corrplot(region, method = "color", type = "upper",na.label = ".",tl.cex = 0.5)
      }
      else
      {
        stop("no correlationmatrix was found", call. = TRUE)
      }
    }
    else
    {
      if(values == "values")
      {
        region <- valuematrix[[1]]
        corrplot(region, method = "color", type = "upper",na.label = ".",tl.cex = 0.5, tl.pos="n")
      }
      else if(values == "upperCI")
      {
        region <- valuematrix[[2]]
        corrplot(region, method = "color", type = "upper",na.label = ".",tl.cex = 0.5, tl.pos="n")
      }
      else if(values == "lowerCI")
      {
        region <- valuematrix[[3]]
        corrplot(region, method = "color", type = "upper",na.label = ".",tl.cex = 0.5, tl.pos="n")
      }
      else
      {
        stop("no correlationmatrix was found", call. = TRUE)
      }
      
    }
    
  }
}
`LD.snpplot` <- function(valuematrix, annotation, values =c("values", "upperCI", "lowerCI"), LD=c("D", "Dprime", "Q", "r", "OR", "MI", "chi2", "Y", "HS"), snp.zoom=FALSE, snp.target, snp.region, unit=c("Mb","Kb", "bp"), legends=TRUE)
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
    
    y=valuematrix[[3]][min.snp:max.snp,match(snp.target, colnames(valuematrix[[1]]))]
    if(unit == "Mb") scale <- 1e6
    else if(unit == "Kb") scale <- 1e3
    else scale <- 1
    x <- (annotation[names(y),"pos"]/scale)
    pos.target <- (annotation[snp.target,"pos"]/scale)
    rbPal <- colorRampPalette(c('red','blue'))
    color <- rbPal(10)[as.numeric(cut(abs(y),breaks = 10))]
    color2 <- rbPal(3)[as.numeric(cut(abs(y),breaks = 3))]
    plot(x = x,y = y, xlab = paste("Position on chromosome in", unit, sep = " "), ylab = LD, main = paste("LD Region for ", snp.target, sep = ""), pch=17,cex=(abs(y)+0.4), col=color, ylim = c(-1,1))
    abline(v=pos.target,lty=2)
    if(legends==TRUE)
    {    
      legend("bottom",legend = c(-1:1), col = c('blue', 'red'), pch = 17, cex=0.7,horiz = TRUE, text.font=2, pt.cex = c(1.4,0.4,1.4))
    }  
  }
  else
  {
    y=valuematrix[[3]] [,match(snp.target, colnames(valuematrix[[1]]))]
    if(unit == "Mb") {scale <- 1e6}
    else if(unit == "Kb") {scale <- 1e3}
    else {scale <- 1}
    x <- annotation[colnames(valuematrix[[1]]),"pos"]/scale
    pos.target <- (annotation[snp.target,"pos"]/scale)
    rbPal <- colorRampPalette(c('red','blue'))
    color <- rbPal(10)[as.numeric(cut(abs(y),breaks = 10))]
    plot(x = x,y = y, xlab = paste("Position on chromosome in", unit, sep = " "), ylab = LD, main = paste("LD Region for ", snp.target, sep = ""), pch=17,cex=(abs(y)+0.4), col=color,ylim = c(-1,1))
    abline(v=pos.target,lty=2)
    if(legends==TRUE)
    {  
      legend("bottom",legend = c(-1:1), col = c('blue', 'red'), pch = 17, cex=0.7,horiz = TRUE, text.font=2, pt.cex = c(1.4,0.4,1.4))  
    }  
  }
}

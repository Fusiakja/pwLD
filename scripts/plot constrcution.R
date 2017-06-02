
#### plot of coverage against minimal haplotype frequency
tabelle <- Jack_se_geno
measures = c("r","D'","Q","Y","HS")
for(measure in measures){
  main_titel="Jackknife leave one out method"
  main_titel= paste(main_titel, measure,sep = ": ")
  pathtitle=paste("Jackknife_loo_geno",measure,sep = "_")
{j=11}
if(measure=="D'")
{j=13}
if(measure=="Q")
{j=15}
if(measure=="Y")
{j=17}
if(measure=="HS")
{j=19}
patyh <- "C:/Users/Jakub F/Desktop/Masterarbeit/Ergebnisse/Plots of Coverage/new plots/"
mypath <- file.path(paste(patyh, pathtitle, ".pdf",sep = ""))
line_colors <- c("black","red","darkgreen","blue","skyblue4")
pdf(file=mypath)
for(k in 1:5)
{
smallest <- c()
for(i in (((k-1)*42)+1):((k*42)))
{
  smallest <- rbind(smallest,c(min(tabelle[i,6:9]),tabelle[i,j]))
}
small_ord <- smallest[order(smallest[,1]),]
if(k==1)
{
  ticks<-c(0,25,50,75,95,100)
  plot(small_ord[,1],small_ord[,2],type = "l", lwd=2.5, ylim = c(0,100), col=line_colors[k] ,xlab = "", ylab = "", main = "", yaxt="n", xaxt="n", cex=3)
  title(main=main_titel, cex.main=1.75)
  axis(1,lwd=1.5,cex.axis=1.2)
  mtext("min. Haplotype Frequency", side=1, line=2.5, cex=1.5)
  axis(2,at=ticks,labels=ticks, lwd=1.5,cex.axis=1.2)
  mtext("% of Values in  95% CI", side=2, line=2.5, cex=1.5)
  abline(h=95,lty = 2)
  }
else
{
  lines(small_ord[,1],small_ord[,2],type = "l", lwd=2.5, ylim = c(0,100), col=line_colors[k])
}
}
legend("bottomright", legend = c("N=100","N=500","N=1000","N=5000","N=10000"), col = line_colors,lwd=c(1.5,1.5))
dev.off()
}


## calcuation of runtime
was<- c("bootstrap", "jackknife")
perfom <- c()
samples <- c(100,500,1000,5000,10000)
for(j in was)
{
  for(l in 0:1)
  {
    for(k in samples)
    {
    start_time <- Sys.time()
    for(i in 1:1000)
    {
      a <- LD.TestCI_geno(hapt = matrize1,snps = HapMap_geno[c(199,203),],LD = c("Dprime", "r","Q", "Y", "HS"),paradigm = "bayes",strategy = j, alpha = 0.05, CI = T, n.sim = 100, returnLDdist = T, seed = F, Ns = k, dirich = c(0.5,0.5,0.5,0.5), gen_hap = "geno", intervall = l, HSweight = 4)
    }
    end_time <- Sys.time()
    difer <- (end_time-  start_time)/1000
    
    perfom <- rbind(perfom, c(j, l, difer))
    }
  }
}

ticks<-samples
line_colors <- c("black","red","darkgreen","blue","skyblue4")
plot(samples,perfom[1:5,3], ylim = c(min( as.numeric(perfom[,3])), as.numeric(max(perfom[,3]))), xaxt="n",yaxt="n", type = "l", col=line_colors[1], xlab = "", ylab = "", main = "", lwd=2.5)
lines(samples,perfom[6:10,3], col=line_colors[2], lwd=2.5, lty=2)
lines(samples,perfom[11:15,3], col=line_colors[3],lwd=2.5, lty=3)
lines(samples,perfom[16:20,3], col=line_colors[4], lwd=2.5, lty=4)
title(main="Computational performance of CI methods", cex.main=1.5)
mtext("sample size N", side=1, line=2.5, cex=1.5)
axis(2,lwd=1.5,cex.axis=1.2)
mtext("Duration in sec", side=2, line=2.5, cex=1.5)
axis(1,at=ticks,labels=ticks, lwd=1.5,cex.axis=1.2)
legend("left",legend = c("Bootstrap q","Bootstrap se","Jackknife pv","Jackknife loo"), col = line_colors,lwd=c(1.5,1.5), ncol = 1,cex = .8,pt.cex = 14,bg="transparent",bty = "n")

#dev.off()




##### table of coverage investigation
tabe <- Boot_se_geno
min_liste <- list()
for(k in 1:5)
{
ergebnis <- c()
Counting <- c("mean coverage","min coverage","max coverage","percentage < 90%")
mittel_cov <- c(mean(tabe[(((k-1)*42)+1):((k*42)),11]),mean(tabe[(((k-1)*42)+1):((k*42)),13]),mean(tabe[(((k-1)*42)+1):((k*42)),15]),mean(tabe[(((k-1)*42)+1):((k*42)),17]),mean(tabe[(((k-1)*42)+1):((k*42)),19]))
min_cov <-  c(min(tabe[(((k-1)*42)+1):((k*42)),11]),min(tabe[(((k-1)*42)+1):((k*42)),13]),min(tabe[(((k-1)*42)+1):((k*42)),15]),min(tabe[(((k-1)*42)+1):((k*42)),17]),min(tabe[(((k-1)*42)+1):((k*42)),19]))
max_cov <- c(max(tabe[(((k-1)*42)+1):((k*42)),11]),max(tabe[(((k-1)*42)+1):((k*42)),13]),max(tabe[(((k-1)*42)+1):((k*42)),15]),max(tabe[(((k-1)*42)+1):((k*42)),17]),max(tabe[(((k-1)*42)+1):((k*42)),19]))
smaller_than <- c(sum(tabe[(((k-1)*42)+1):((k*42)),11] <90)/42*100,sum(tabe[(((k-1)*42)+1):((k*42)),13] <90)/42*100,sum(tabe[(((k-1)*42)+1):((k*42)),15]<90)/42*100,sum(tabe[(((k-1)*42)+1):((k*42)),17]<90)/42*100,sum(tabe[(((k-1)*42)+1):((k*42)),19]<90)/42*100)
ergebnis <- rbind(mittel_cov, min_cov, max_cov,smaller_than)
rownames(ergebnis) <- c(Counting[1],Counting[2],Counting[3],Counting[4])
colnames(ergebnis) <- c(colnames(tabe)[11],colnames(tabe)[13],colnames(tabe)[15],colnames(tabe)[17],colnames(tabe)[19])
min_liste[[k]] <- ergebnis
}




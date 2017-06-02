# construct tables from Scholz aet al
scholztablef <- function(OR, shape){
  x <- log(sqrt(OR))
  if(shape == 1)
  {
    y <- 0
    z <- 0
  }
  if(shape == 2)
  {
    y <- x
    z <- -x
  }

  if(shape == 3)
  {
    y <- 10
    z <- -y
  }

  if(shape == 4)
  {
    y <- 10
    z <- -x
  }

  if(shape == 5)
  {
    y <- 10
    z <- y
  }

  freqs <- (1/(exp(x+y+z)+exp(x)+exp(y)+exp(z)))*c(exp(x+y+z), exp(y), exp(z), exp(x))
  return(freqs)
}
# construct tables with defined marginals
ratiosfunc <- function(All_1, All_2, OR){
g <- function(poo) poo^2*(1-OR)+poo*((1-All_1)-All_2+OR*(All_1+All_2))-(OR*All_1*All_2)
br <- c(0,min(All_1,All_2))
if( OR > 1 || OR < 1)
{
  poo_1 <- (-((1-All_1)-All_2+OR*(All_1+All_2))+sqrt((((1-All_1)-All_2+OR*(All_1+All_2))^2)-(4*(1-OR)*(-(OR*All_1*All_2)))))/(2*(1-OR))
  poo_2<- (-((1-All_1)-All_2+OR*(All_1+All_2))-sqrt((((1-All_1)-All_2+OR*(All_1+All_2))^2)-(4*(1-OR)*(-(OR*All_1*All_2)))))/(2*(1-OR))

  if(poo_2>0 && poo_2 < min(All_1,All_2))
  {poo <- poo_2}
  else if(poo_1>0 && poo_1 < min(All_1,All_2))
  {poo <- poo_1}
  else
  {poo <- 0}
}
else
{
  OR=1.0000000000001
  poo_1 <- (-((1-All_1)-All_2+OR*(All_1+All_2))+sqrt((((1-All_1)-All_2+OR*(All_1+All_2))^2)-(4*(1-OR)*(-(OR*All_1*All_2)))))/(2*(1-OR))
  poo_2<- (-((1-All_1)-All_2+OR*(All_1+All_2))-sqrt((((1-All_1)-All_2+OR*(All_1+All_2))^2)-(4*(1-OR)*(-(OR*All_1*All_2)))))/(2*(1-OR))

  if(poo_2>0 && poo_2 < min(All_1,All_2))
  {poo <- poo_2}
  else if(poo_1>0 && poo_1 < min(All_1,All_2))
  {poo <- poo_1}
  else
  {poo <- 0}
}

#poo <- uniroot(g,br, tol = 6e-12)$root

po1 <- All_1-poo
p1o <- All_2-poo
p11 <- 1-po1-p1o-poo

z = matrix(0,2,2)
z[1,] = c(poo,po1)
z[2,] = c(p1o,p11)

return(z)
}


# number of Haplotype counts
haploCounts <- c(100,500,1000,5000,10000)
# odds ratios investigated
odds_ratio <- c(1,2,5,10,20,50,100)
# initialize
res <- c()
bres <- c()
hs <- c()
all_sums <- c()
prozent <- 0

# Which marginal frequencies
allelfrqs <- matrix(c(0.4,0.6,0.6,0.6,0.2,0.1),ncol = 2,nrow = 3, byrow = TRUE)
# Which method: jackkinfe, bootstrap 
for(jack_boot in c("jackknife"))
{
  for(count in haploCounts)
  {
    for(ratio in odds_ratio)
    {
      for(bia in 1:2)
      {
        for(dia in 1:3)
        { 
          # calc haplotype freuqncy table
          if(bia ==1)
          {
            hap <- scholztablef(OR = ratio,shape = dia)
            matrize <- t(matrix(c(hap[1],hap[3],hap[1]+hap[3],hap[2],hap[4],hap[2]+hap[4],hap[1]+hap[2],hap[3]+hap[4], sum(hap)), nrow = 3, byrow = T))
          }
          else
          {
            hap <- ratiosfunc(All_1 =allelfrqs[dia,1], All_2 =  allelfrqs[dia,2],OR = ratio)
            matrize <- t(matrix(c(hap[1],hap[3],hap[1]+hap[3],hap[2],hap[4],hap[2]+hap[4],hap[1]+hap[2],hap[3]+hap[4], sum(hap)), nrow = 3, byrow = T))
          }
          # set row and colnames
          rownames(matrize) <- c("0", "1","Sum")
          colnames(matrize) <- c("0", "1","Sum")
          # calc OR
          oar <- matrize[1,1]*matrize[2,2]/(matrize[1,2]*matrize[2,1])
          (sqrt(oar)-1)/(sqrt(oar)+1)
          # calc original LD value
          should_be <- estLD(matrize, LD=c("Dprime", "r","Q", "Y", "HS"), HSweight = 4)
          
          bvars <- c(0,0,0,0,0)
          boots <- list()
          vars <- c()
          prozentDp <- 0
          prozentHS <- 0
          prozentQ <- 0
          prozentr2 <- 0
          prozentY <- 0
          bprozentDp <- 0
          bprozentHS <- 0
          bprozentQ <- 0
          bprozentr2 <- 0
          bprozentY <- 0
          not_in_ci_Dp <- 0
          not_in_ci_hs <- 0
          not_in_ci_q <- 0
          not_in_ci_r <- 0
          not_in_ci_y <- 0
          not_nas <- 0
          nas <- 0
          bnot_in_ci_Dp <- 0
          bnot_in_ci_hs <- 0
          bnot_in_ci_q <- 0
          bnot_in_ci_r <- 0
          bnot_in_ci_y <- 0
          bnot_nas <- 0
          bnas <- 0
          # random seed
          set.seed(1)
          # draw samples
          zap <- rmultinom(1000,count,prob = hap)
          for(j in 1:dim(zap)[2])
          {
            ############## From haplotype Level #######################################
            # sap <- zap[,j]
            # matrize1 <- t(matrix(c(sap[1],sap[3],sap[1]+sap[3],sap[2],sap[4],sap[2]+sap[4],sap[1]+sap[2],sap[3]+sap[4], sum(sap)), nrow = 3, byrow = T))
            # a <- LD.TestCI(hapt = matrize1,snps = HapMap_geno[c(199,203),],LD = c("Dprime", "r","Q", "Y", "HS"),paradigm = "bayes",strategy = "jackknife", alpha = 0.05, CI = T, n.sim = 1000, returnLDdist = T, seed = T, Ns = count, dirich = c(0.5,0.5,0.5,0.5), intervall=0,  HSweight = 4)
            ############# From genotype level ########################################
            sap <- zap[,j]/count
            matrize <-t(matrix(c(sap[1],sap[3],sap[1]+sap[3],sap[2],sap[4],sap[2]+sap[4],sap[1]+sap[2],sap[3]+sap[4], sum(sap)), nrow = 3, byrow = T))

            matrize1 <- t(matrix(c(hap[1]^2,2*hap[3]*hap[1],hap[3]^2,
                                   2*hap[2]*hap[1],(2*hap[3]*hap[2])+(2*hap[1]*hap[4]),2*hap[3]*hap[4],
                                   hap[2]^2,2*hap[2]*hap[4], hap[4]^2), nrow = 3, byrow = T))

            matrize1 <- matrix(c(matrize1[1,1],matrize1[1,2],matrize1[1,3],sum(matrize1[1,])
                                 ,matrize1[2,1],matrize1[2,2],matrize1[2,3],sum(matrize1[2,])
                                 ,matrize1[3,1],matrize1[3,2],matrize1[3,3],sum(matrize1[3,])
                                 ,sum(matrize1[,1]),sum(matrize1[,2]), sum(matrize1[,3]), sum(matrize1)), ncol = 4,nrow = 4, byrow = TRUE)
            a <- LD.TestCI_geno(hapt = matrize1,snps = HapMap_geno[c(199,203),],LD = c("Dprime", "r","Q", "Y", "HS"),paradigm = "bayes",strategy = "jackknife", alpha = 0.05, CI = T, n.sim = 1000, returnLDdist = T, seed = F, Ns = count, dirich = c(0.5,0.5,0.5,0.5), gen_hap = "geno", intervall = 0, HSweight = 4)
            
            # check if LD value in CI
            if(!anyNA(a$Q))
            {
              if(!((anyNA(a$CI))||(any(is.nan(a$CI)))))
              {
                

                if(should_be[1] >= a$CI[1]  & should_be[1] <= a$CI[6])
                {
                  prozentDp <- prozentDp+1
                }
                else
                {
                  not_in_ci_Dp <- not_in_ci_Dp+1
                }
                if(should_be[2] >= a$CI[2]  & should_be[2] <= a$CI[7])
                {
                  prozentr2 <- prozentr2+1
                }
                else
                {
                  not_in_ci_r <- not_in_ci_r+1
                }
                if(should_be[3] >= a$CI[3]  & should_be[3] <= a$CI[8])
                {
                  prozentQ <- prozentQ+1
                }
                else
                {
                  not_in_ci_q <- not_in_ci_q+1
                }
                if(should_be[4] >= a$CI[4]  & should_be[4] <= a$CI[9])
                {
                  prozentY <- prozentY+1
                }
                else
                {
                  not_in_ci_y <- not_in_ci_y+1
                }
                if(should_be[5] >= a$CI[5]  & should_be[5] <= a$CI[10])
                {
                  prozentHS <- prozentHS+1
                }
                else
                {
                  not_in_ci_hs <- not_in_ci_hs+1
                }
                not_nas <- not_nas+1
              }
              else
              {
                nas <- nas+1 
              }  
            }
          }
          if(not_nas>0)
          {
            # calc procentage
            prozentHS <- (prozentHS/not_nas)*100
            prozentDp <- (prozentDp/not_nas)*100
            prozentr2 <- (prozentr2/not_nas)*100
            prozentQ <- (prozentQ/not_nas)*100
            prozentY <- (prozentY/not_nas)*100
          }
          else
          {
            prozentHS <- 0
            prozentDp <- 0
            prozentr2 <- 0
            prozentQ <- 0
            prozentY <- 0
          }
          # constr table
          res <- rbind(res, c( oar,hap[1]+hap[2],hap[3]+hap[4],hap[1]+hap[3], hap[2]+hap[4], hap,should_be[2],prozentr2, should_be[1],prozentDp,should_be[3],prozentQ,should_be[4],prozentY,should_be[5],prozentHS))

        }
      }
    }
  }
}

#rownames(res) <- rep(rep(rep(odds_ratio,c(rep(3, length(haploCounts)))), length(samples)),3)
# set colnames of result table
colnames(res) <- c("OR","po.","p1.","p.o","p.1","f00", "f01", "f10", "f11","true r", "r","true D'","D'","true Q", "Q","true Y", "Y","true HS", "HS")
#colnames(bres) <- c("Odds Ratio","Ao+","A1+","A+o","A+1","f00", "f01", "f10", "f11","true r", "mean var r","r","true D'", "mean var D'","D'","true Q", "mean var Q", "Q","true Y", "mean var Y","Y","true HS", "mean var HS","HS")

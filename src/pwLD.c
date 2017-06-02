/* Karsten Krug, 2008-06-12 
 
 C code for package  pwLD
 
 */

# include<R.h>
# include<Rmath.h>
# include <Rdefines.h>
# include <Rinternals.h>
# include "my_solve_poly.h"
# include "pwLD.h"
# include "LDmeasures.h"
# include "distribution_statistics.h"



/*############################################################
 given	 the pairwise genotypic frequencies the function estimates
 the haplotype frequencies by the fix-point method
 
 genoFreq		- 4x4 matrix of genotypic frequencies with
 marginals
 haploFreq		- 3x3 matrix of haplotype frequencies with
 marginals
 tol			- epsilon, tolerance value in order to check
 equality of double values
 - used to determine which solutions of the cubic
 polynomial are within
 the valid interval, defined by the given
 allele frequencies
 ( is a pointer for reasons of the C interface
 of R, using function .C() )
 */
void estimateHaploFreq( double genoFreq[4][4], double haploFreq[3][3], double
                          *tol, int *digits, double pooHat[3] )
{
  double k=0, p00=0, p01=0, p10=0, p11=0, af[2]={0,0};
  double c1=0, c2=0, c3=0, c4=0, a=0, b=0, c=0, x[3]={0,0,0}, pooHatML=0, *xptr=0, L[3]={0,0,0};	
  int i=0, j=0, intervalIdx[3]={0,0,0}, validSolNumb=0, validSolNumbIdx[3]={0,0,0};
  
  for(i=0; i < 3; i++) pooHat[i] = NA_REAL; 
  
  /* 		frequency double heterozygote cases, k 
   since the cubic polynom is defined on p_ij 's of unambiguously inferred 
   haplotype frequencies (n_ij/2N ), whereas the parameter k is a 
   genotypic frequency  N_i^j / N, this frequencies has to be deivided
   by 2 */
  k = genoFreq[1][1]/2;
  
  /* estimate allele frequencies from genoytpe frequencies */
  alleleFreq(genoFreq, af);
  
  
  /* if there are no double heterozygotes */
  if(k == 0) 
  {
    pooHatML = genoFreq[0][0] + 0.5* genoFreq[0][1] + 0.5*genoFreq[1][0];
    pooHat[0] = pooHatML;
    
  }	
  
  else
  {	
    /* non-amiguous haplotype frequencies */
    p00 = genoFreq[0][0] + 0.5* genoFreq[0][1] + 0.5*genoFreq[1][0]; 
    p01 = genoFreq[0][2] + 0.5* genoFreq[0][1] + 0.5*genoFreq[1][2];
    p10 = genoFreq[2][0] + 0.5* genoFreq[1][0] + 0.5*genoFreq[2][1];
    p11 = genoFreq[2][2] + 0.5* genoFreq[1][2] + 0.5*genoFreq[2][1];
    
    /* the cubic polynomial */
    c1 = -2*pow(k, 2);
    c2 = (3*pow(k, 2) - k*p00 - k*p11 + k*p10 + k*p01); 
    c3 = (k*p00 + k*p11 - p00*p11 - p01*p10 - k*p01 - k*p10 - pow(k,2));
    c4 = (p00*p11);	
    
    // 		a = fprec(c2/c1, *digits);
    // 		b = fprec(c3/c1, *digits);
    // 		c = fprec(c4/c1, *digits);
    
    a= c2/c1; b= c3/c1; c=c4/c1;
    
    /* solve the polynomial */
    xptr = x;
    
    /* cubic function */
    if( c != 0)
    {	
      /* solve the cubic polynomial */
      gsl_poly_solve_cubic( a, b, c, xptr, xptr+1, xptr+2, *tol, *digits);
      
      
      //Rprintf("\na=%f, b=%f, c=%f\nx0=%f, x1=%f, x2=%f\n", a, b, c, x[0], x[1], x[2]);
    }
    /* quadratic function */
    else 
    {
      x[2] = 0;
      gsl_poly_solve_quadratic( 1, a, b, xptr,xptr+1);	
    }
    
    /* calculate the solutions of haplotype frequency p00*/
    for(i=0; i < 3; i++)
      pooHat[i] = (p00 + k* x[i]);
    
    /* determine which solutions are within the valid interval */
    whichSolutionInterval( pooHat,  intervalIdx,  af,  *tol );
    //Rprintf("\nvalid interval: [%f, %f]\nintervallIdx: %d, %d, %d\npooHat: %f, %f, %f\n", fmax2(0, af[0]+af[1]-1), fmin2(af[0],af[1]), intervalIdx[0], intervalIdx[1], intervalIdx[2], pooHat[0], pooHat[1], pooHat[2]);
    
    /* how many and which solutions are within the interval ? */
    for(i=0;i<3;i++) if(intervalIdx[i] == 1)
    {
      validSolNumb++;
      validSolNumbIdx[i] = 1;
    }
    /* if there is an unique solution*/
    if( validSolNumb == 1){ for(i=0;i<3;i++) if(validSolNumbIdx[i]== 1 ) pooHatML = pooHat[i];}
    else
    {	double ML=0;
      
      /* calculate the likelihood score of each solution */
      whichSolutionML(genoFreq, pooHat, 3,  validSolNumbIdx, af, L);
      
      ML = L[0];
      pooHatML = pooHat[0];
      
      /* maximum likelihood */
      for(i=1;i<3; i++)
      {
        if( fmax2(ML, L[i]) > ML ) 
        {
          pooHatML = pooHat[i];
          ML = L[i];
        }
      }
    }
  }
  
  /* contingency table of haplotype frequencies */
  haploFreq[0][0] = pooHatML;
  haploFreq[0][1] = fmax2(0, af[0] - pooHatML);
  haploFreq[0][2] = af[0];
  haploFreq[1][0] = fmax2(0, af[1] - pooHatML);
  haploFreq[1][1] = fmax2(0, 1 - haploFreq[0][0] - haploFreq[0][1] - haploFreq[1][0]);
  haploFreq[1][2] = fmax2(0,1-af[0]);
  haploFreq[2][0] = af[1];
  haploFreq[2][1] = fmax2(0,1- af[1]);
  haploFreq[2][2] = 0;
  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      haploFreq[2][2]+= haploFreq[i][j];
  /* rescale the table */
  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      haploFreq[i][j] /= haploFreq[2][2];
  
}


/* #################################################
 nine different likelihood terms with 
 respect to nine genotypic frequencies 
 */
double term00(double Poo) { return 2*log(Poo); }
double term01(double Poo, double pox){ return (log(2)+log(Poo)+log(pox-Poo));}
double term02(double Poo, double pox){ return 2*log(pox-Poo); }
double term10(double Poo, double pxo){ return (log(2)+log(Poo)+log(pxo-Poo)); }
double term11(double Poo, double pox, double pxo){ return log(2*Poo*(1-pox-pxo+Poo) + 2*(pox-Poo)*(pxo-Poo)); }
double term12(double Poo, double pox, double pxo){ return (log(2)+log(pox-Poo)+log(1-pox-pxo+Poo));}
double term20(double Poo, double pxo){ return (2*log(pxo-Poo)); }
double term21(double Poo, double pox, double pxo){ return (log(2)+log(pxo-Poo)+log(1-pox-pxo+Poo));}
double term22(double Poo, double pox, double pxo){ return 2*log(1-pox-pxo+Poo);
}


/*####################################################
 determines the solution that is maximum likelihood
 
 
 poo				- points to a vector of solutions
 solNumb			- lenght of the vector
 validSolNumbIdx	- vector of the same length indicating which solutions, pooHat,
 are within the valid interval
 - solutions outside the interval have a likelihood score of -Inf
 af				- allele frequencies
 Lscore			- points to the vector that holds the Likelihood scores
 
 */
void whichSolutionML( double genoFreq[4][4], double *poo, int solNumb, int *validSolNumbIdx, double *af, double *Lscore)
{
  double L[solNumb], tmp=0;
  int i, j, s;
  
  /* initialize the likelihood scores of valid solution with zero,
   and with -Inf for  the  solution outside the valid interval*/
  for(i=0;i<solNumb;i++)  
  {
    if(validSolNumbIdx[i] == 1) L[i] = 0;
    else L[i] = R_NegInf;
  }
  
  /* evaluate the likelihood at the three solutions */
  for(i=0; i< 3; i++)
    for(j=0; j<3; j++)
    {	for(s=0; s < solNumb; s++)
    {
      if(validSolNumbIdx[s] == 1)
      {	
        if(i == 0 && j ==0)
        {
          tmp = genoFreq[0][0]*term00(poo[s]);
          if( !isnan( tmp ) ) L[s] += tmp;
          //Rprintf("%d-%d: poo%d: %f, L=%f\t ",i,j, s, tmp, L[s]);		
        }
        if(i == 0 && j ==1) 
        {
          tmp = genoFreq[0][1]*term01(poo[s], af[0]);
          if(!isnan(tmp)) L[s] += tmp;
          //Rprintf("%d-%d: poo%d: %f, L=%f\t ",i,j ,s, tmp, L[s]);
        }
        if(i == 0 && j ==2)
        {
          tmp = genoFreq[0][2]*term02(poo[s], af[0]);
          if(!isnan(tmp)) L[s] += tmp;
          //Rprintf("%d-%d: poo%d: %f, L=%f\t ",i,j, s, tmp, L[s]);
        }
        if(i == 1 && j ==0)
        {
          tmp = genoFreq[1][0]*term10(poo[s], af[1]);
          if(!isnan(tmp)) L[s] += tmp;	
          //Rprintf("%d-%d: poo%d: %f,L=%f\t ",i,j, s, tmp, L[s]);
        }
        if(i == 1 && j ==1)
        {
          
          tmp=genoFreq[1][1]*term11(poo[s], af[0], af[1]);
          if(!isnan(tmp)) L[s] += tmp;	
          //Rprintf("%d-%d: poo%d: %f,L=%f\t ",i,j, s, tmp, L[s]);
        }
        if(i == 1 && j ==2)
        {
          
          tmp=genoFreq[1][2]*term12(poo[s], af[0], af[1]);
          if(!isnan(tmp)) L[s] += tmp;
          //Rprintf("%d-%d: poo%d: %f, L=%f\t ",i,j, s, tmp, L[s]);
        }
        if(i == 2 && j ==0)
        {
          tmp=genoFreq[2][0]*term20(poo[s], af[1]);
          if(!isnan(tmp)) L[s] += tmp;
          //Rprintf("%d-%d: poo%d: %f,L=%f\t ",i,j, s, tmp, L[s]);
        }
        if(i == 2 && j ==1)
        {
          tmp=genoFreq[2][1]*term21(poo[s], af[0], af[1]);
          if(!isnan(tmp)) L[s] += tmp;
          //Rprintf("%d-%d: poo%d: %f,L=%f\t ",i,j, s, tmp, L[s]);
        }
        if(i == 2 && j ==2) 
        {
          tmp=genoFreq[2][2]*term22(poo[s], af[0], af[1]);
          if(!isnan(tmp)) L[s] += tmp;
          //Rprintf("%d-%d: poo%d: %f,L=%f\t ",i,j, s, tmp, L[s]);
        }
      }
    }
    
    }	
    for(i=0; i < solNumb; i++)
      Lscore[i] = L[i];
  
}


/*######################################
 allele frequencies
 
 genoFreq		- 4x4 matrix (with marginals) of genotypic frequencies
 af			- vector of length two store the allele frequencies
 */
void alleleFreq( double genoFreq[4][4], double af[2])
{
  af[0]=0; af[1]=0;
  
  af[0] += (genoFreq[0][3] + 0.5*genoFreq[1][3]);
  af[1] += (genoFreq[3][0] + 0.5*genoFreq[3][1]);
  
}


/*###########################
 check which solutions of the cubic
 are within the valid interval
 
 poo		- points to a vector of solutions
 idx		- points to a vector of the same length as poo
 - the vector is used to indicate which solution are within the valid interval
 if idx[i] == 1 -> poo[i] is a valid solutions
 if idx[i] == 0 -> poo[i] is outside the valid interval
 tol		- epsilon, tolerance value in order to check equality of double values
 - used to determine which solutions of the cubic polynomial are within
 the valid interval, defined by the given allele frequencies
 */
void whichSolutionInterval(double *poo, int *idx, double *af, double tol)
{
  double validInt[2];
  int i;
  
  /* determine the valid interval of paramter Poo */
  validInt[0] = fmax2(0, af[0] + af[1] - 1);
  validInt[1] = fmin2(af[0], af[1]);
  
  for(i=0; i < 3; i++)
  {	
    if(  ((poo[i] - validInt[0]) < 0)  && ( fabs(validInt[0] -poo[i]) < tol ) )
    {				
      poo[i] = validInt[0];		
      
    }
    else if( ((validInt[1] - poo[i]) < 0) && (fabs(validInt[1] -poo[i]) < tol)  )
    {	
      poo[i] = validInt[1];
      
    }
    if( (validInt[0] <= poo[i]) && (poo[i] <= validInt[1]) )
      idx[i]=1;
    
    else idx[i] =0;
    
  }
}





/*###################################################
 estimation of LD measures
 
 tab		- 3x3 matrix of haplotype frequencies with marginals
 what		- character vector specifying the LD measures to
 estimate
 numb		- number of LD measures (length of 'what')
 LD		- points to a vector of length 'numb' that stores the
 result 
 */
void estimateLD(double tab[3][3], char **what, int *numb,  double *LD, double *HSweight, double *exponent)
{
  int i=0;
  for(i=0; i<*numb; i++)
  {		
    if( strcmp(what[i], "Q") == 0  )
      LD[i] = Q(tab);
    if( strcmp(what[i], "Dprime") == 0 )
      LD[i] = Dprime(tab);
    if( strcmp(what[i], "r") == 0 )
      LD[i] = R(tab);
    if( strcmp(what[i], "OR") == 0 )
      LD[i] = oddsratio(tab);
    if( strcmp(what[i], "D") == 0 )
      LD[i] = D(tab);
    if( strcmp(what[i], "MI") == 0 )
      LD[i] = MI(tab);
    if( strcmp(what[i], "chi2") == 0 )
      LD[i] = ChiSquareHaplotypeFrequencies(tab);
    if( strcmp(what[i], "Y") == 0 )
      LD[i] = Y(tab);
    if( strcmp(what[i], "HS") == 0 )
      LD[i] = HS(tab, HSweight, exponent);
    
    
  }
}

/*###################################################
 
 estimation of genotypic frequencies
 
 - genoA		- points to a vetcor of genotypes
 - genoB		- points to a vetcor of genotypes
 - code		- points to a vector of length 4 representing the coding used for the genotypes
 - N			- sample size, the length of genoA / genoB
 - paradigm		- character indicating the paradigm to use in order to estimate frequencies
 - either "freq" or "bayes"
 - Dir			- vector of length 9 representing the shape parameter of the Dirichlet prioir
 - genoFreq
 - Neff		- points to a variable, that holds the effective sample size, i.e. without missing values, after calling this function 
 
 */
void estimateGenoFreq( char **genoA, char **genoB, char **code, int *N,  char **paradigm, double Dir[9], double genoFreq[4][4], double genoCounts[4][4], int *Neff )
{
  
  double DirSum=0;
  int i=0, j=0;
  
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
    {	genoFreq[i][j] = 0;
      genoCounts[i][j] = 0;
    }
    /* the sum of the 9 Dirichlet shape parameter */
    for(i=0; i<9; i++){ DirSum += Dir[i];}	
    
    /* genotype counts */
    for(i=0; i < *N; i++)
    {
      if( (strcmp( genoA[i], code[0]) == 0) &&  ( strcmp(genoB[i], code[0]) == 0 ))  genoCounts[0][0] += 1;  
      if( (strcmp( genoA[i], code[0]) == 0) &&  ( strcmp(genoB[i], code[1]) == 0 ))  genoCounts[0][1] += 1;
      if( (strcmp( genoA[i], code[0]) == 0) &&  ( strcmp(genoB[i], code[2]) == 0 ))  genoCounts[0][2] += 1;  
      if( (strcmp( genoA[i], code[1]) == 0) &&  ( strcmp(genoB[i], code[0]) == 0 ))  genoCounts[1][0] += 1; 
      if( (strcmp( genoA[i], code[1]) == 0) &&  ( strcmp(genoB[i], code[1]) == 0 ))  genoCounts[1][1] += 1;
      if( (strcmp( genoA[i], code[1]) == 0) &&  ( strcmp(genoB[i], code[2]) == 0 ))  genoCounts[1][2] += 1;  
      if( (strcmp( genoA[i], code[2]) == 0) &&  ( strcmp(genoB[i], code[0]) == 0 ))  genoCounts[2][0] += 1; 
      if( (strcmp( genoA[i], code[2]) == 0) &&  ( strcmp(genoB[i], code[1]) == 0 ))  genoCounts[2][1] += 1; 
      if( (strcmp( genoA[i], code[2]) == 0) &&  ( strcmp(genoB[i], code[2]) == 0 ))  genoCounts[2][2] += 1; 
    }
    
    
    /*genoCounts[0][0]=counts[0][0];
     genoCounts[0][1]=counts[0][1];
     genoCounts[0][2]=counts[0][2];
     genoCounts[1][0]=counts[1][0];
     genoCounts[1][1]=counts[1][1];
     genoCounts[1][2]=counts[1][2];
     genoCounts[2][0]=counts[2][0];
     genoCounts[2][1]=counts[2][1];
     genoCounts[2][2]=counts[2][2];*/
    
    
    /* marginal distributions */
    for(i=0;i<3;i++) 
    {
      genoCounts[0][3] += genoCounts[0][i];
      genoCounts[1][3] += genoCounts[1][i];
      genoCounts[2][3] += genoCounts[2][i];
      
      genoCounts[3][0] += genoCounts[i][0];
      genoCounts[3][1] += genoCounts[i][1];
      genoCounts[3][2] += genoCounts[i][2];
    }
    
    
    /* sample size without missings */
    genoCounts[3][3] = genoCounts[3][0] + genoCounts[3][1] + genoCounts[3][2];
  *Neff = genoCounts[3][3];
} 

void estimateFrequencies(double genoFreq[3][3], double haploFreq[3][3], int *Neff, double Dir[4], int *mc, char **paradigm)
{
  
  double N = *Neff;
  int monte_carlo = *mc;
  if ((strcmp( paradigm[0], "freq") == 0)) 
  {
    for(int i=0;i<4;i++)
    {
      for(int j=0;j<4;j++)
      {	
        haploFreq[i][j] = genoFreq[i][j]/N;
      }
    }
  }
  if ((strcmp( paradigm[0], "bayes") == 0)) 
  {
    haploFreq[0][0] = (Dir[0]+genoFreq[0][0])/(4*Dir[0]+N);
    haploFreq[0][1] = (Dir[1]+genoFreq[0][1])/(4*Dir[1]+N);
    haploFreq[1][0] = (Dir[2]+genoFreq[1][0])/(4*Dir[2]+N);
    haploFreq[1][1] = (Dir[3]+genoFreq[1][1])/(4*Dir[3]+N);
    
    haploFreq[2][0] = haploFreq[0][0]+haploFreq[1][0];
    haploFreq[2][1] = haploFreq[0][1]+haploFreq[1][1];
    haploFreq[0][2] = haploFreq[0][0]+haploFreq[0][1];
    haploFreq[1][2] = haploFreq[1][0]+haploFreq[1][1];
    haploFreq[2][2] = haploFreq[0][2]+haploFreq[1][2];
  }
  
  if ((strcmp( paradigm[0], "fullbayes") == 0)) 
  {
    double d[4] = {Dir[0]+genoFreq[0][0], Dir[1]+genoFreq[0][1], Dir[2]+genoFreq[1][0], Dir[3]+genoFreq[1][1]};
    double out[4] ={0,0,0,0};
    double mc_mean[4] = {0,0,0,0};
    for(int k=0; k< monte_carlo; k++)
    {
      out[0]=0;
      out[1]=0;
      out[2]=0;
      out[3]=0;
      rdirichlet(d,4,out);
      
      for (int g=0; g< 4; g++)
      {
        mc_mean[g]+=out[g];
      }
    }
    
    for (int g=0; g< 4; g++)
    {
      mc_mean[g]/= monte_carlo;
    }
    
    haploFreq[0][0] = mc_mean[0];
    haploFreq[0][1] = mc_mean[1];
    haploFreq[1][0] = mc_mean[2];
    haploFreq[1][1] = mc_mean[3];
    
    
    haploFreq[2][0] = haploFreq[0][0]+haploFreq[1][0];
    haploFreq[2][1] = haploFreq[0][1]+haploFreq[1][1];
    haploFreq[0][2] = haploFreq[0][0]+haploFreq[0][1];
    haploFreq[1][2] = haploFreq[1][0]+haploFreq[1][1];
    haploFreq[2][2] = haploFreq[0][2]+haploFreq[1][2];
  }
  
}



/*##################################
 
 estimation of pairwise haplotype frequencies and LD measures
 of all possible marker combinations contained in 'data'
 
 data
 nrow
 ncol
 LD
 LDnumb
 code
 paradigm
 Dir
 MAF
 tol	
 LDmatPtr
 strategy - which CI method
 ci - computed or not
 mc - how many MC iterations
 alpha - significance level
 clow, cup - lower and upper bounds of CI
 nSim - number of Bootstrap samples
 LDdist - out for internal poperperties (e. g. LDs for Bootstrap samples)
 vars - output for standard error of LDs of CIs
 intervall - which CI method
 
 */
void LDall(  char **data, int *nrow, int *ncol, char **LD, int *LDnumb,  char
               **code, char **paradigm, double *Dir, double *MAF, double *tol, int *digits, double *LDmatPtr, double *HSweight, int *ci, int *mc, char **strategy, double *alpha, double *cilow, double *ciup, int *nsim, double *LDdist, double *vars, int *intervall)
{
  int i=0, j=0, k=0, m=0, n=0, N=0, countEntries=0, Nentries=(*nrow)*(*nrow-1)/2;
  char *genoSnp1[*ncol], *genoSnp2[*ncol];
  double genoFreq[4][4], genoCounts[4][4], haploFreq[3][3], haploCounts[3][3], *LDptr=0, LDvalue[*LDnumb], pooHat[3]; 
  LDptr = LDvalue;
  double allelfreq[2];
  double tabels[*nsim][9];
  double Liste[*nsim];
  *LDdist=10;
  for(int s=0; s<*nsim; s++)
  {
    Liste[s]=0;
  }
  //Rprintf("%f ", *HSweight);
  
  /* loop over all possible marker combinations (the upper triangular matrix) */	
  for(i=0; i<(*nrow-1); i++)
  {	for(j=(i+1); j<(*nrow); j++)
  {
    /* extract two rows of the data matrix, representing the genotypes of two markers*/
    for(k=0;k<*ncol; k++) 
    {	genoSnp1[k] = data[i+k*(*nrow)];
      genoSnp2[k] = data[j+k*(*nrow)];
    } 
    
    /* initialize the matrices of genotype and haplotype frequencies */
    for(n=0; n<4;n++)
      for(m=0; m<4; m++)
      {	
        genoFreq[n][m] = NA_REAL;
        genoCounts[n][m] = NA_REAL;
        
        if((n < 3) && (m < 3))
        {	haploFreq[n][m] = NA_REAL;
          haploCounts[n][m] = NA_REAL;
        }
      }
      
      /* estimate genotypic frequencies*/
      estimateGenoFreq( genoSnp1, genoSnp2, code, ncol, paradigm,  Dir,  genoFreq, genoCounts, &N);
    
    //estimateGenoFreq( data[i], data[j], code, ncol, paradigm,  Dir,  genoFreq, Nptr );
    
    for(n=0;n<4;n++)
    {
      for(m=0;m<4;m++)
      {
        genoFreq[n][m]=genoCounts[n][m]/genoCounts[3][3];
      }
    }
    
    /* estimate allele frequencies from genotypic
     frequencies*/
    //alleleFreq(genoFreq, allelfreq);
    
    allelfreq[0]=genoFreq[0][3] + 0.5*genoFreq[1][3];
    allelfreq[1]=genoFreq[3][0] + 0.5*genoFreq[3][1];
    
    /* check wether the allele frequencies are above the MAF
     threshold */
    //if( (af[0] > *MAF) && (af[0] < (1-*MAF)) && (af[1] > *MAF) && (af[1] < (1-*MAF)) )	
    if( (allelfreq[0] >= *MAF) && (allelfreq[0] <= (1-*MAF)) && (allelfreq[1] >= *MAF) && (allelfreq[1] <= (1-*MAF)) )	
    {
      
      /* estimate haplotype frequencies */
      estimateHaploFreq( genoFreq, haploCounts, tol, digits, pooHat );
      
      for(n=0;n<3;n++)
      {
        for(m=0;m<3;m++)
        {
          haploCounts[n][m]= haploCounts[n][m]*2*N;
        }
      }
      
      int hap_N=2* N;
      estimateFrequencies(haploCounts, haploFreq, &hap_N, Dir, mc, paradigm);
      
      /*estimate LD*/
      double exponent=0;
      estimateLD(haploFreq, LD, LDnumb, LDptr, HSweight, &exponent);
      //Rprintf("%f ", haploFreq[0][0]);
      
      /* store the result */
      for(k=0; k < *LDnumb; k++)
      {
        //Rprintf("%f ", LDvalue[k]);
        LDmatPtr[ k*Nentries+ countEntries ] = LDvalue[k];
        //Rprintf("%f ", LDvalue[k]);
      }
      
      if(*ci==1)
      {
        for(int s=0; s<*nsim; s++)
        {
          LDdist[s]=0;
        }
        
        //Rprintf("%f ", allelfreq[0]);
        //alleleFreq(genoFreq, allelfreq);
        for(k=0; k < *LDnumb; k++)
        {
          int Num=1;
          double vars=1;
          int mc=1;
          confidenceGenoInterval(genoCounts,genoFreq, &N, paradigm, nsim, &LD[k], &Num, tol, digits, Liste, HSweight, alpha, strategy, &cilow[k*Nentries+ countEntries], &ciup[k*Nentries+ countEntries], tabels, Dir, &vars, intervall, &mc);
          //Rprintf("%d ", Num);
          
          //alleleFreq(genoFreq, allelfreq);
          //Rprintf("%f ", haploFreq[0][0]);
          if(ciup[k*Nentries+ countEntries]<=-1)
          {
            ciup[k*Nentries+ countEntries]=-1;
          }
          if(cilow[k*Nentries+ countEntries]<=-1)
          {
            cilow[k*Nentries+ countEntries]=-1;
          }
          if(ciup[k*Nentries+ countEntries]>=1)
          {
            ciup[k*Nentries+ countEntries]=1;
          }
          if(cilow[k*Nentries+ countEntries]>=1)
          {
            cilow[k*Nentries+ countEntries]=1;
          }
          
          for(int s=0; s<*nsim; s++)
          {
            //Rprintf("%f ", LDdist[s]);
            Liste[s]=0;
            //Rprintf("%f ", LDdist[s]);
            
          }
        }
      }
    }
    /* if any allele frequency does not suffice the MAF
     criteria, the respective LD value is set to 'NA' */
    else 
    {	
      for(k=0; k < *LDnumb; k++) 
      {	LDmatPtr[ k*Nentries + countEntries ] = NA_REAL;
        
        //LDmatPtr[k][countEntries] =NA_REAL;
        if(*ci==1)
        {
          if(ciup[k*Nentries+ countEntries]<=-1)
          {
            ciup[k*Nentries+ countEntries]=NA_REAL;
          }
          if(cilow[k*Nentries+ countEntries]<=-1)
          {
            cilow[k*Nentries+ countEntries]=NA_REAL;
          }
        }
      }
    }
    
    /* increment the number of entries of the upper
     triangular matrix */
    countEntries++;
  }	
  
  }
}

/*#########################################
 
 NOT USED IN THE CURRENT IMPLEMENTATION
 
 generates a bootstrap sample from genotypic frequencies by
 multinomial sampling from 'genoFreq' with sample size 'N'
 
 genoFreq		- 4x4 matrix of genotypic frequencies with marginals
 N			- sample size
 genoFreqBS	- 4x4 matrix to store the bootstrapped genotype frequencies
 */


void bootstrapHapFreq( double genoFreq[3][3], int *N, double genoFreqBS[3][3], char **paradigm )
{
  double genoFreqVector[4];
  int genoCountsBSVector[4], i=0, j=0, count=0;
  //double hapfreq[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  int Neff =2* *N;
  /*  represent genotypic frequencies  as vector */
  for(i=0; i<2; i++)
  {	for(j=0;j<2;j++)
  {
    genoFreqVector[count] = genoFreq[i][j]/Neff;
    genoCountsBSVector[count] = 0;
    count++;
    //Rprintf("%f ", genoFreq[i][j]);
  }	
  //Rprintf("\n");
  }
  
  /* draw a multinomial sample with probabilities  */
  GetRNGstate();
  rmultinom(Neff, genoFreqVector, 4, genoCountsBSVector);
  PutRNGstate();
  
  for(int l=0; l<4; l++)
  {
    if( genoCountsBSVector[l]==NAN)
    {bootstrapHapFreq(genoFreq, N, genoFreqBS, paradigm);} 
  }
  
  //for(i=0;i<9;i++) Rprintf("%d ", genoCountsBSVector[i]);
  //Rprintf("\n");
  
  /* 	- the drawn sample is stored in  'genoCountsBSVector' 
   - convert this vector into a 4x4 matrix of
   genotypic frequencies with marginals 
   */
  count=0;
  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
    {
      genoFreqBS[i][j] = (double)genoCountsBSVector[count];
      count++;
    }
    
  genoFreqBS[2][0] = genoFreqBS[0][0]+genoFreqBS[1][0];
  genoFreqBS[2][1] = genoFreqBS[0][1]+genoFreqBS[1][1];
  genoFreqBS[0][2] = genoFreqBS[0][0]+genoFreqBS[0][1];
  genoFreqBS[1][2] = genoFreqBS[1][0]+genoFreqBS[1][1];
  genoFreqBS[2][2] = genoFreqBS[0][2]+genoFreqBS[1][2];
  
  //estimateFrequencies(genoFreqBS, hapfreq, &Neff,  Dir, mc, paradigm);
  //Rprintf(" %f", hapfreq[0][0]);
  /* check wether the allele frequencies are valid */
  //alleleFreq(genoFreqBS, af);
  if( (genoFreqBS[2][0] > 0 && genoFreqBS[2][0] < genoFreqBS[2][2]) && ( genoFreqBS[0][2] > 0 && genoFreqBS[0][2] < genoFreqBS[2][2] )) 
  {
    return;
  }
  //if the allele frequencies are NOT valid, i.e. 0 or 1, the function is called recursively 
  else
  {
    bootstrapHapFreq(genoFreq, N, genoFreqBS, paradigm);
  }
  
}

/*##################################################### 
 
 NOt used in the currend implementation
 
 - estimation of frequentistic confidence intervals by bootstrapping
 - the function returns the distributions of LD measures over the simulations
 - the confidence intervals as quantiles of that distributions are calculated in R
 
 genoFreq		- 4x4 matrix of genotypic frequencies with marginals
 N			- sample size
 nSim			- number of bootstrap replicates to simulate
 LD			- character vector of LD measures
 LDnumb		- number of LD measures to estimate (length of vector 'LD')
 tol			- epsilon for robust comparison of float values
 LDdist		- vector to store the distributions of the LD measures over the simulations
 */
void confidenceInterval( double genoCounts[3][3],double genoFreq[3][3], int *N, char **paradigm,int *nSim, char **LD, int *LDnumb, double *tol, int *digits,  double *LDdisti, double *HSweight, double *alpha, char **strategy,  double *cilow, double *ciup, double tables[*nSim][4], double Dir[4], double vars[*LDnumb], int *intervall, int *mc)
{
  int i=0, j=0;
  double genoFreqBS[3][3], haploFreq[3][3]={{0,0,0},{0,0,0},{0,0,0}}; 
  double LDvaluei[*LDnumb];

  
  int sims=*nSim;
  /* initialize the matrix of bootstrapped genotype frequencies*/
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      genoFreqBS[i][j]=0;
  
  if((strcmp( strategy[0], "jackknife") == 0))
  {
    double exponent=0;
    double original_exponent=0;
    double original_Y=0;
    int y_num=1;
    
    double exponent_list[*LDnumb*(2*(*N))];
    char *org_Y="Y";
    for(int g=0; g<(*LDnumb*(2*(*N)));g++)
    {
      exponent_list[g]=0;
    }
    exponent_list[0]=exponent_list[0];
    double observed_LDs[*LDnumb];
    double Y_list[*LDnumb*(2*(*N))];
    double Y_2=0;
    double zero_weight=0;
    estimateLD(genoFreq, &org_Y, &y_num, &original_Y, &zero_weight, &original_exponent);
    estimateLD(genoFreq, LD, LDnumb, observed_LDs, HSweight, &original_exponent);
    int jacktimes=0;
    int element=0;
    int tab=0;
    int new_hap_N=(2* *N);
    for(int a=0; a<2; a++)
    {
      for(int b=0; b<2; b++)
      {
        
        for(int j=0;j<*LDnumb;j++)
        {LDvaluei[j] = 0;}
        
        jacktimes=(int) round(genoCounts[a][b]);
        if((int) round(genoCounts[a][b]) > 0)
        {
          
          genoCounts[a][b]=(int) round(genoCounts[a][b])-1;
          
          genoCounts[2][0] = genoCounts[0][0]+genoCounts[1][0];
          genoCounts[2][1] = genoCounts[0][1]+genoCounts[1][1];
          genoCounts[0][2] = genoCounts[0][0]+genoCounts[0][1];
          genoCounts[1][2] = genoCounts[1][0]+genoCounts[1][1];
          genoCounts[2][2] = genoCounts[0][2]+genoCounts[1][2];
          
          
          tables[tab][0]=genoCounts[0][0];
          tables[tab][1]=genoCounts[0][1];
          tables[tab][2]=genoCounts[1][0];
          tables[tab][3]=genoCounts[1][1];
          
          int hap_n_minus_one=new_hap_N-1;
          
          for(int j=0;j<*LDnumb;j++)
          {
            
            estimateFrequencies(genoCounts, haploFreq, &hap_n_minus_one, Dir, mc, paradigm);
            /*estimate LD*/
            estimateLD(haploFreq, &org_Y, &y_num, &Y_2, &zero_weight, &exponent);
            estimateLD(haploFreq, LD, LDnumb, LDvaluei, HSweight, &exponent);
          }
          
          for(int k=element; k<element+jacktimes; k++)
          {
            for(int f=0; f<*LDnumb; f++)
            {
              LDdisti[f*(2* *N)+k]=LDvaluei[f];
              //Rprintf(" %f", LDvaluei[f]);
              exponent_list[f*(2* *N)+k]=exponent;
              Y_list[f* (2* *N) + k]=Y_2;
            }
          }
          Y_list[0]=Y_list[0];
          genoCounts[a][b]=genoCounts[a][b]+1;
          genoCounts[2][0] = genoCounts[0][0]+genoCounts[1][0];
          genoCounts[2][1] = genoCounts[0][1]+genoCounts[1][1];
          genoCounts[0][2] = genoCounts[0][0]+genoCounts[0][1];
          genoCounts[1][2] = genoCounts[1][0]+genoCounts[1][1];
          genoCounts[2][2] = genoCounts[0][2]+genoCounts[1][2];
          element+=jacktimes;
          tab++;
        }
      }
    }
    
    for(int e=0;e<*LDnumb;e++)
    {
      if(*intervall == 0)
      { 
       
        /*double LDs[(2*(*N))];
        double expos_for_LDs[(2*(*N))];
        double Ys[(2*(*N))];
        for(int h=0; h<2*(*N); h++)
        {
          LDs[h]=0;
          expos_for_LDs[h]=0;
          Ys[h]=0;
        }
        int intList=0;
        //draw values from large List
        
        for(int g=e*(2*(*N)); g<e*(2*(*N))+(2*(*N)); g++)
        {
          LDs[intList]=LDdisti[g];
          expos_for_LDs[intList]=exponent_list[g];
          Ys[intList]=Y_list[g];
          intList++;
        }
        double var_ps=0;
        if((strcmp( LD[e], "HS") == 0))
        {
          intList=0;
          //draw values from large List
          double sum_ps=0;
          for(int a=0; a<2*(*N); a++)
          {
            sum_ps+= Ys[a];
          }
          
          sum_ps=sum_ps/(2*(*N));
          
          var_ps=0;
          
          for(int a=0; a<2*(*N); a++)
          {
            var_ps+= pow(Ys[a]-sum_ps,2);
          }
          
          double new_ns=2*(*N);
          //Rprintf(" %f", new_ns);
          double new_ns_min_one=new_ns-1;
          double ratio=(new_ns_min_one)/new_ns;
          //Rprintf(" %f", ratio);
          var_ps=pow(var_ps*ratio,0.5);
          var_ps=var_ps*original_exponent*pow(fabs(original_Y),(original_exponent-1));
        }
        else{
          //Mean of values
          double sum_ps=0;
          for(int a=0; a<2*(*N); a++)
          {
            sum_ps+=LDs[a];
          }
          sum_ps=sum_ps/(2*(*N));
          
          //calculate standard errror of values
          var_ps=0;        
          for(int a=0; a<2*(*N); a++)
          {
            var_ps+= pow(LDs[a]-sum_ps,2);
          }
          
          double new_ns=2*(*N);
          //Rprintf(" %f", new_ns);
          double new_ns_min_one=new_ns-1;
          double ratio=(new_ns_min_one)/new_ns;
          //Rprintf(" %f", ratio);
          var_ps=pow(var_ps*ratio,0.5);
        }
        
        //calculate cis
        double stud_t=0;
        double alphafull=1-(*alpha/2);
        double sim_minus_1=(2* *N)-1;
        int logger=0;
        int tail=1;
        stud_t= Rf_qt(alphafull, sim_minus_1,tail,logger);
        
        double Dupper = observed_LDs[e] + 1.964 * (var_ps);
        double Dlower = observed_LDs[e] - 1.964 * (var_ps);
        
        if(Dlower < -1){Dlower=-1;}
        if(Dupper > 1){Dupper=1;}	    
        
        cilow[e]=Dlower;
        ciup[e]=Dupper;*/
        
        
        /////////////////////////////7
        
        /*double LDs[(2*(*N))];
        double expos_for_LDs[(2*(*N))];
        double Ys[(2*(*N))];
        for(int h=0; h<2*(*N); h++)
        {
        LDs[h]=0;
        expos_for_LDs[h]=0;
        Ys[h]=0;
        }
        int intList=0;
        //draw values from large List
        for(int g=e*(2*(*N)); g<e*(2*(*N))+(2*(*N)); g++)
        {
        LDs[intList]=LDdisti[g];
        expos_for_LDs[intList]=exponent_list[g];
        LDdisti[g]=exponent_list[g];
        Ys[intList]=Y_list[g];
        intList++;
        }
        double var_ps=0;
        if((strcmp( LD[e], "HS") == 0))
        {
          intList=0;
          //draw values from large List
          
          for(int a=0; a<2*(*N); a++)
          {
            LDs[a]= Ys[a];
          }
        }
        //Mean of values
        double sum_ps=0;
        for(int a=0; a<2*(*N); a++)
        {
        sum_ps+=LDs[a];
        }
        sum_ps=sum_ps/(2*(*N));
        
        //calculate standard errror of values

        for(int a=0; a<2*(*N); a++)
        {
        var_ps+= pow(LDs[a]-sum_ps,2);
        }
        
        double new_ns=2*(*N);
        //Rprintf(" %f", new_ns);
        double new_ns_min_one=new_ns-1;
        double ratio=(new_ns_min_one)/new_ns;
        //Rprintf(" %f", ratio);
        var_ps=pow(var_ps*ratio,0.5);
        //calculate cis
        double stud_t=0;
        double alphafull=1-(*alpha/2);
        double sim_minus_1=(2* *N)-1;
        int logger=0;
        int tail=1;
        stud_t= Rf_qt(alphafull, sim_minus_1,tail,logger);
        
        double Dupper = observed_LDs[e] + 1.964 * (var_ps);
        double Dlower = observed_LDs[e] - 1.964 * (var_ps);
        
        if((strcmp( LD[e], "HS") == 0))
        {
        Dupper=sign(Dupper)*pow(fabs(Dupper), original_exponent);
        Dlower=sign(Dlower)*pow(fabs(Dlower), original_exponent);
        }
        
        if(Dlower < -1){Dlower=-1;}
        if(Dupper > 1){Dupper=1;}	    
        
        cilow[e]=Dlower;
        ciup[e]=Dupper;*/
        
        //////////////////      Alte Methode
        long double LDs[(2*(*N))];
        for(int h=0; h<(2*(*N)); h++)
        {
          LDs[h]=0;
        }
        int intList=0;
        for(int g=e*(2*(*N)); g<e*(2*(*N))+(2*(*N)); g++)
        {
          //LDs[intList]=atanh(LDdisti[g]);// - atanh(observed_LDs[e]));
          if((strcmp( paradigm[0], "freq") == 0))
          {
            LDs[intList] = ((2*(*N))*observed_LDs[e]) - ((2*(*N-1))*LDdisti[g]);
          }
          else
          {
            LDs[intList] = ((2*(*N))*atanh(observed_LDs[e])) - ((2*(*N-1))*atanh(LDdisti[g]));
          }
          
          //Rprintf ("%f ", LDdisti[g]);
          //Rprintf(" %f", observed_LDs[e]);
          
          intList++;
        }
        long double sum_ps=0;
        for(int a=0; a<(2*(*N)); a++)
        {
          sum_ps+=LDs[a];
        }
        sum_ps= sum_ps/(2*(*N));
        //Rprintf(" %f", sum_ps);
        
        //long double sum_ps_2=0;
        //long double sum_ps_3=0;
         
        //for(int a=0; a<*nSim; a++)
        //{
        //  if(LDs[a] >=0)
        //  {sum_ps_2+=LDs[a];}
        //  else
        //  {sum_ps_3+=fabs(LDs[a]);}
        //}
        
        //long double sum_ps_4=sum_ps_2-sum_ps_3;
        //calculate standard errror of values
        long double var_ps=0;
        //var_ps=pow(sum_ps_4,2)+pow((2*(*N))*sum_ps,2)-2*(2* *N)*sum_ps*sum_ps_4;
        
        
        for(int a=0; a<(2*(*N)); a++)
        {
          var_ps+=(pow((LDs[a]-sum_ps),2));
        }
        //Rprintf(" %f", var_ps);
        
        var_ps=pow(var_ps*1/((2*(*N)-1)*(2*(*N))),0.5);
        //double new_ns=2*(*N);
        //Rprintf(" %f", new_ns);
        //double new_ns_min_one=new_ns-1;
        //double ratio=(new_ns_min_one)/new_ns;
        //Rprintf(" %f", ratio);
        //var_ps=sqrt(var_ps*ratio);
        //Rprintf(" %f", var_ps);
        
        vars[e]=tanh(var_ps);
        //double Dupper =(sum_ps + 1.96 *sqrt((var_ps)/(*N)));//;
        //double Dlower =(sum_ps - 1.96 *sqrt((var_ps)/(*N)));//(new_hap_N));
        
        double stud_t=0;
        double alphafull=1-(*alpha/2);
        double sim_minus_1=(2* *N)-1;
        int logger=0;
        int tail=1;
        //Rf_qt(alphafull,sim_minus_1);
        stud_t= Rf_qt(alphafull, sim_minus_1,tail,logger);
        
      //  long double Dupper=tanh(atanh(observed_LDs[e])+stud_t*var_ps);
      //  long double Dlower=tanh(atanh(observed_LDs[e])-stud_t*var_ps);
        long double Dupper=0;
        long double Dlower=0;
        if((strcmp( paradigm[0], "freq") == 0))
        {
          Dupper=observed_LDs[e]+stud_t*var_ps;
          Dlower=observed_LDs[e]-stud_t*var_ps;        
        }
        else
        {
          Dupper=tanh(atanh(observed_LDs[e])+stud_t*var_ps);
          Dlower=tanh(atanh(observed_LDs[e])-stud_t*var_ps);        
        }
        
        if(Dlower < -1){Dlower=-1;}
        if(Dupper > 1){Dupper=1;}
        
        cilow[e]=Dlower;//D_prime - 1.644854*sqrt(Var_Zapata);
        ciup[e]=Dupper;
      }
      else
      {
        long double LDs[(2*(*N))];
        for(int h=0; h<(2*(*N)); h++)
        {
          LDs[h]=0;
        }
        int intList=0;
        for(int g=e*(2*(*N)); g<e*(2*(*N))+(2*(*N)); g++)
        {
          if((strcmp( paradigm[0], "freq") == 0))
          {
            LDs[intList]=LDdisti[g];
          }
          else
          {
            LDs[intList]=atanh(LDdisti[g]);
          }
          intList++;
          
        }
        long double sum_ps=0;
        for(int a=0; a<(2*(*N)); a++)
        {
          sum_ps+=LDs[a];
        }
        sum_ps= sum_ps/(2*(*N));


  
        double var_ps=0;
        for(int a=0; a<(2*(*N)); a++)
        {
          var_ps+=(pow((LDs[a]-sum_ps),2));
        }
        //Rprintf(" %f", var_ps);
        
        //var_ps=sqrt(var_ps*1/((2*(*N)-1)*(2*(*N))));
        double new_ns=2*(*N);
        //Rprintf(" %f", new_ns);
        double new_ns_min_one=new_ns-1;
        double ratio=(new_ns_min_one)/new_ns;
        //Rprintf(" %f", ratio);
        var_ps=pow(var_ps*ratio,0.5);
        //Rprintf(" %f", var_ps);
        
        vars[e]=tanh(var_ps);
        //double Dupper =(sum_ps + 1.96 *sqrt((var_ps)/(*N)));//;
        //double Dlower =(sum_ps - 1.96 *sqrt((var_ps)/(*N)));//(new_hap_N));
        
        double stud_t=0;
        double alphafull=1-(*alpha/2);
        double sim_minus_1=(2* *N)-1;
        int logger=0;
        int tail=1;
        //Rf_qt(alphafull,sim_minus_1);
        stud_t= Rf_qt(alphafull, sim_minus_1,tail,logger);
        
        //long double Dupper=tanh(atanh(observed_LDs[e])+stud_t*var_ps);
        //long double Dlower=tanh(atanh(observed_LDs[e])-stud_t*var_ps);
        long double Dupper=0;
        long double Dlower=0;
        if((strcmp( paradigm[0], "freq") == 0))
        {
          Dupper=observed_LDs[e]+stud_t*var_ps;
          Dlower=observed_LDs[e]-stud_t*var_ps;        
        }
        else
        {
          Dupper=tanh(atanh(observed_LDs[e])+stud_t*var_ps);
          Dlower=tanh(atanh(observed_LDs[e])-stud_t*var_ps);        
        }
        
        if(Dlower < -1){Dlower=-1;}
        if(Dupper > 1){Dupper=1;}
        
        cilow[e]=Dlower;//D_prime - 1.644854*sqrt(Var_Zapata);
        ciup[e]=Dupper;
        
      }
      
    }
    
  }
  if((strcmp( strategy[0], "bootstrap") == 0))
  {
    double exponent=0;
    double original_exponent=0;
    double original_Y=0;
    int y_num=1;
    
    double exponent_list[*LDnumb*(*nSim)];
    char *org_Y="Y";
    for(int g=0; g<(*LDnumb*(*nSim));g++)
    {
      exponent_list[g]=0;
    }
    exponent_list[0]=exponent_list[0];
    double observed_LDs[*LDnumb];
    double Y_list[*LDnumb*(*nSim)];
    double Y_2=0;
    double zero_weight=0;
    estimateLD(genoFreq, &org_Y, &y_num, &original_Y, &zero_weight, &original_exponent);
    estimateLD(genoFreq, LD, LDnumb, observed_LDs, HSweight, &original_exponent);

    int hap_N=2* *N;
    /* loop over the number of simulations */
    for(i=0; i<*nSim; i++ )
    {
      for(j=0;j<*LDnumb;j++)
      {LDvaluei[j] = 0;}
      
      /* generate a bootstraped 3x3 table of genotypic frequencies */
      bootstrapHapFreq(genoCounts, N, genoFreqBS, paradigm);
      
      tables[i][0]=genoFreqBS[0][0];
      tables[i][1]=genoFreqBS[0][1];
      tables[i][2]=genoFreqBS[1][0];
      tables[i][3]=genoFreqBS[1][1];
      /* estimate haplotype frequencies */
      //estimateHaploFreq(genoFreqBS, haploFreq, tol, digits, pooHat);
      
      estimateFrequencies(genoFreqBS, haploFreq, &hap_N, Dir, mc, paradigm);
      
      /* estimate LD */
      estimateLD(haploFreq, &org_Y, &y_num, &Y_2, &zero_weight, &exponent);
      estimateLD(haploFreq, LD, LDnumb, LDvaluei, HSweight, &exponent);
      /* store the result  in a vector */
      for(j=0;j<*LDnumb;j++)
      {
        LDdisti[ j* (*nSim) + i] = LDvaluei[j];
        exponent_list[j* (*nSim) + i]=exponent;
        Y_list[j* (*nSim) + i]=Y_2;
      }
      Y_list[0]=Y_list[0];
    }
    
    for(int e=0;e<*LDnumb;e++)
    {
      if(*intervall==0)
      {
        double LDs[*nSim];
        for(int h=0; h<*nSim; h++)
        {
          LDs[h]=0;
        }
        
        int intList=0;
        for(int g=e*(*nSim); g<e*(*nSim)+(*nSim); g++)
        {
          LDs[intList]=LDdisti[g];
          intList++;
        }
        sort(LDs, sims);
        double alphas=*alpha;
        double alpha_u = alphas/2;
        double tolle=quantile(sims, LDs, alpha_u);
        //Rprintf ("%f ", tolle);
        //digits=quantile(LDdist, *nSim);
        double alpha_l = 1-(alphas/2);
        double digi=quantile(sims, LDs, alpha_l);
        //Rprintf ("%f ", digi);
        // char strat=strategy;
        //LDdist=0;
        cilow[e]=tolle;
        ciup[e]=digi;
      }
      else
      {
        
        long double LDs[*nSim];
        for(int h=0; h<(*nSim); h++)
        {
          LDs[h]=0;
        }
        int intList=0;
        for(int g=e*(*nSim); g<e*(*nSim)+(*nSim); g++)
        {
          if((strcmp( paradigm[0], "freq") == 0))
          {
            LDs[intList]=LDdisti[g];
          }
          else
          {
            LDs[intList]=atanh(LDdisti[g]);
          }
          intList++;
        }
        long double sum_ps=0;
        for(int a=0; a<(*nSim); a++)
        {
          sum_ps+=LDs[a];
        }
        sum_ps= sum_ps/(*nSim);
        
        
        
        double var_ps=0;
        for(int a=0; a<(*nSim); a++)
        {
          var_ps+=pow((LDs[a]-sum_ps),2);
        }
        //Rprintf(" %f", var_ps);
        var_ps = pow(var_ps/(*nSim-1),0.5);
        
        vars[e]=tanh(var_ps);
        //double Dupper =(sum_ps + 1.96 *sqrt((var_ps)/(*N)));//;
        //double Dlower =(sum_ps - 1.96 *sqrt((var_ps)/(*N)));//(new_hap_N));
        
        double stud_t=0;
        double alphafull=1-(*alpha/2);
        double sim_minus_1=(2* *N)-1;
        int logger=0;
        int tail=1;
        //Rf_qt(alphafull,sim_minus_1);
        stud_t= Rf_qt(alphafull, sim_minus_1,tail,logger);
        long double Dupper=0;
        long double Dlower=0;
        if((strcmp( paradigm[0], "freq") == 0))
        {
          Dupper=observed_LDs[e]+stud_t*var_ps;
          Dlower=observed_LDs[e]-stud_t*var_ps;        
        }
        else
        {
          Dupper=tanh(atanh(observed_LDs[e])+stud_t*var_ps);
          Dlower=tanh(atanh(observed_LDs[e])-stud_t*var_ps);        
        }
        
  //      long double Dupper=tanh(atanh(observed_LDs[e])+stud_t*var_ps);
  //      long double Dlower=tanh(atanh(observed_LDs[e])-stud_t*var_ps);
        if(Dlower < -1){Dlower=-1;}
        if(Dupper > 1){Dupper=1;}
        
        cilow[e]=Dlower;//D_prime - 1.644854*sqrt(Var_Zapata);
        ciup[e]=Dupper;
        
        ///////////////////////////
        
        /*double LDs[(*nSim)];
        double expos_for_LDs[(*nSim)];
        double Ys[(*nSim)];
        for(int h=0; h<*nSim; h++)
        {
          LDs[h]=0;
          expos_for_LDs[h]=0;
          Ys[h]=0;
        }
        int intList=0;
        //draw values from large List
        
        for(int g=e*(*nSim); g<e*(*nSim)+(*nSim); g++)
        {
          LDs[intList] = LDdisti[g];
          expos_for_LDs[intList]=exponent_list[g];
          Ys[intList]=Y_list[g];
          intList++;
        }
        double var_ps=0;
        if((strcmp( LD[e], "HS") == 0))
        {
          intList=0;
          //draw values from large List
          double sum_ps=0;
          for(int a=0; a<*nSim; a++)
          {
            sum_ps+= Ys[a];
          }
          
          sum_ps=sum_ps/(*nSim);
          
          var_ps=0;
          
          for(int a=0; a<*nSim; a++)
          {
            var_ps+= pow(Ys[a]-sum_ps,2);
          }
          
          var_ps = pow(var_ps/(*nSim-1),0.5);
          var_ps=var_ps*original_exponent*sign(original_Y)* pow(fabs(original_Y),original_exponent-1);
        }
        else{
        //Mean of values
        double sum_ps=0;
        for(int a=0; a<*nSim; a++)
        {
          sum_ps+=LDs[a];
        }
        sum_ps=sum_ps/(*nSim);
        
        //calculate standard errror of values
        var_ps=0;        
        for(int a=0; a<*nSim; a++)
        {
          var_ps+= pow(LDs[a]-sum_ps,2);
        }
        
        var_ps = pow(var_ps/(*nSim-1),0.5);
        }
        
        //calculate cis
        double stud_t=0;
        double alphafull=1-(*alpha/2);
        double sim_minus_1=(2* *N)-1;
        int logger=0;
        int tail=1;
        stud_t= Rf_qt(alphafull, sim_minus_1,tail,logger);
        
        double Dupper = observed_LDs[e] + stud_t * (var_ps);
        double Dlower = observed_LDs[e] - stud_t * (var_ps);

        if(Dlower < -1){Dlower=-1;}
        if(Dupper > 1){Dupper=1;}	    

        cilow[e]=Dlower;
        ciup[e]=Dupper;*/
        
        /*double LDs[(*nSim)];
        double expos_for_LDs[(*nSim)];
        double Ys[(*nSim)];
        for(int h=0; h<(*nSim); h++)
        {
        LDs[h]=0;
        expos_for_LDs[h]=0;
        Ys[h]=0;
        }
        int intList=0;
        //draw values from large List
        
        for(int g=e*(*nSim); g<e*(*nSim)+(*nSim); g++)
        {
        LDs[intList]=LDdisti[g];
        expos_for_LDs[intList]=exponent_list[g];
        Ys[intList]=Y_list[g];
        intList++;
        }
        double var_ps=0;
        if((strcmp( LD[e], "HS") == 0))
        {
        intList=0;
        //draw values from large List
        double sum_ps=0;
        for(int a=0; a<(*nSim); a++)
        {
        sum_ps+= Ys[a];
        }
        
        sum_ps=sum_ps/((*nSim));
        
        var_ps=0;
        
        for(int a=0; a<(*nSim); a++)
        {
        var_ps+= pow(Ys[a]-sum_ps,2);
        }
        
        double new_ns=(*nSim);
        //Rprintf(" %f", new_ns);
        double new_ns_min_one=new_ns-1;
        double ratio=(new_ns_min_one)/new_ns;
        //Rprintf(" %f", ratio);
        var_ps=pow(var_ps*ratio,0.5);
        var_ps=var_ps*original_exponent*pow(fabs(original_Y),(original_exponent-1));
        }
        else{
        //Mean of values
        double sum_ps=0;
        for(int a=0; a<(*nSim); a++)
        {
        sum_ps+=LDs[a];
        }
        sum_ps=sum_ps/(*nSim);
        
        //calculate standard errror of values
        var_ps=0;        
        for(int a=0; a<(*nSim); a++)
        {
        var_ps+= pow(LDs[a]-sum_ps,2);
        }
        
        double new_ns=(*nSim);
        //Rprintf(" %f", new_ns);
        double new_ns_min_one=new_ns-1;
        double ratio=(new_ns_min_one)/new_ns;
        //Rprintf(" %f", ratio);
        var_ps=pow(var_ps*ratio,0.5);
        }
        
        //calculate cis
        double stud_t=0;
        double alphafull=1-(*alpha/2);
        double sim_minus_1=(2*(*N))-1;
        int logger=0;
        int tail=1;
        stud_t= Rf_qt(alphafull, sim_minus_1,tail,logger);
        
        double Dupper = observed_LDs[e] + 1.964 * (var_ps);
        double Dlower = observed_LDs[e] - 1.964 * (var_ps);
        
        if(Dlower < -1){Dlower=-1;}
        if(Dupper > 1){Dupper=1;}	    
        
        cilow[e]=Dlower;
        ciup[e]=Dupper;*/
        
        /*double LDs[(*nSim)];
        double expos_for_LDs[(*nSim)];
        double Ys[(*nSim)];
        for(int h=0; h<*nSim; h++)
        {
          LDs[h]=0;
          expos_for_LDs[h]=0;
          Ys[h]=0;
        }
        int intList=0;
        //draw values from large List
        
        for(int g=e*(*nSim); g<e*(*nSim)+(*nSim); g++)
        {
          LDs[intList] = LDdisti[g];
          expos_for_LDs[intList]=exponent_list[g];
          Ys[intList]=Y_list[g];
          intList++;
        }
        if((strcmp( LD[e], "HS") == 0))
        {
          intList=0;
          //draw values from large List
  
          for(int a=0; a<*nSim; a++)
          {
            LDs[a]= Ys[a];
          }
        }
        //Mean of values
        long double sum_ps=0;
        for(int a=0; a<*nSim; a++)
        {
          sum_ps+=LDs[a];
        }
        sum_ps=sum_ps/(*nSim);
        
        //calculate standard errror of values
        long double var_ps=0;
        
        for(int a=0; a<*nSim; a++)
        {
          var_ps+= pow(LDs[a]-sum_ps,2);
        }
        
        var_ps = pow(var_ps/(*nSim-1),0.5);

        
        //calculate cis
        double stud_t=0;
        double alphafull=1-(*alpha/2);
        double sim_minus_1=(2* *N)-1;
        int logger=0;
        int tail=1;
        stud_t= Rf_qt(alphafull, sim_minus_1,tail,logger);
        
        double Dupper = observed_LDs[e] + 1.964 * (var_ps);
        double Dlower = observed_LDs[e] - 1.964 * (var_ps);
        
        if((strcmp( LD[e], "HS") == 0))
        {
        Dupper=sign(Dupper)*pow(fabs(Dupper), original_exponent);
        Dlower=sign(Dlower)*pow(fabs(Dlower), original_exponent);
        }
        
        if(Dlower < -1){Dlower=-1;}
        if(Dupper > 1){Dupper=1;}	    

        cilow[e]=Dlower;
        ciup[e]=Dupper;*/
      }
    }
  }
 
  return;
}



/*############################################
 
 generates a new set of genotypic frequencies assuming 
 Dirichlet prior
 
 */
void sampleDirichletGenoFreq( int genoCounts[3][3], double *Dir, double genoFreqRS[3][3] )
{
  double shapeDirichlet[4]={0,0,0,0}, dirichletGenoFreq[4]={0,0,0,0};
  int i=0, j=0, count=0; // N=genoCounts[3][3];
  
  
  /*  determine the shape parameters of the Dirichlet distribution to sample from, i.e. n_ij + Dir_ij */
  for(i=0; i<2; i++)
  {	for(j=0;j<2;j++)
  {
    shapeDirichlet[count] = genoCounts[i][j] + Dir[count];
    count++;
    //Rprintf("%f ", genoCounts[i][j]);
  }	
  //Rprintf("\n");
  }
  
  /* sample from Dirichlet( n_ij + Dir_ij )  */
  rdirichlet( shapeDirichlet, 4, dirichletGenoFreq);
  
  
  //for(i=0;i<9;i++) Rprintf("%f ", dirichletGenoFreq[i]);
  //Rprintf("\n");
  
  /* 	- the drawn sample is stored in  'dirichletGenoFreq' 
   - convert this vector into a 4x4 matrix of
   genotypic frequencies with marginals 
   */
  count=0;
  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
    {
      genoFreqRS[i][j] = dirichletGenoFreq[count];
      count++;
    }
    
    /* initialize the marginals of the table with zero */	
    /*for(i=0; i<2; i++) 
     { 	
     genoFreqRS[i][2]=0;
     genoFreqRS[2][i]=0;
     }*/
    
  genoFreqRS[2][0] = genoFreqRS[0][0]+genoFreqRS[1][0];
  genoFreqRS[2][1] = genoFreqRS[0][1]+genoFreqRS[1][1];
  genoFreqRS[0][2] = genoFreqRS[0][0]+genoFreqRS[0][1];
  genoFreqRS[1][2] = genoFreqRS[1][0]+genoFreqRS[1][1];
  genoFreqRS[2][2] = genoFreqRS[0][2]+genoFreqRS[1][2];
  
  /* estimate marginal distributions */
  /*for(i=0;i<2;i++)
   {
   genoFreqRS[0][2] +=genoFreqRS[0][i];
   genoFreqRS[1][2] +=genoFreqRS[1][i];
   
   genoFreqRS[2][0] +=genoFreqRS[i][0];
   genoFreqRS[2][1] +=genoFreqRS[i][1];
   }
   
   genoFreqRS[2][2] = genoFreqRS[2][1] + genoFreqRS[2][0];*/
  
  /* check wether the allele frequencies are valid */
  //alleleFreq(genoFreqRS, af);
  //if( (af[0] > 0 && af[0] < 1) && ( af[1] > 0 && af[1] < 1 )) return;
  /* if the allele frequencies are NOT valid, i.e. 0 or 1, the function is called recursively */
  //else
  //sampleDirichletGenoFreq(genoCounts, Dir, genoFreqRS);
  
  return;
}
/*###########################################
 
 estimation of Bayesian Credible Intervals 
 
 */

void credibleInterval(int genoCounts[3][3], int *nSim, double Dir[4], char **LD, int *LDnumb, double *tol, int *digits, double *LDdist, double *HSweight, double *alpha, double *cilow, double *ciup)
{
  int i=0, j=0;
  double genoFreqRS[3][3], LDvalue[*LDnumb];
  int sims = *nSim;
  /* initialize the matrix of simulated genotype frequencies*/
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      genoFreqRS[i][j]=0;
  
  for(i=0; i< *nSim;i++)
  {
    /* simulate genotypic frequencies assuming a Dirichlet prior with shape parameter n_ij + Dir_ij */
    sampleDirichletGenoFreq(genoCounts, Dir, genoFreqRS);
    
    
    /* estimate haplotype frequencies */
    //estimateFrequencies(genoFreqRS, haploFreq, &Neff, Dir, mc, &paradigm);
    
    /* estimate LD */
    double exponent=0;
    estimateLD(genoFreqRS, LD, LDnumb, LDvalue, HSweight, &exponent);
    
    /* store the result  in a vector */
    for(j=0;j<*LDnumb;j++)
      LDdist[ j* (*nSim) + i] = LDvalue[j];
  }
  
  
  for(j=0;j<*LDnumb;j++)
  {
    double LDs[*nSim];
    int intList=0;
    for(int g=j; g<*LDnumb* *nSim; g=g+*LDnumb)
    {
      LDs[intList]=LDdist[g];
      //Rprintf ("%d ", g);
      intList++;
    }
    
    sort(LDs, sims);
    
    double alphas=*alpha;
    double tolle=0;//quantile(sims, LDdist, alpha_u);
    double digi=0;//quantile(sims, LDdist, alpha_l);
    double lengths[(int)(alphas/0.001)];
    double ups[(int)(alphas/0.001)];
    double lows[(int)(alphas/0.001)];
    for(int l=0; l<alphas/0.001; l++)
    {
      lengths[l]=-100;
      ups[l]=-100;
      lows[l]=-100;
    }
    int counts=0;
    double max =0;
    double min =100;
    for(double k=0; k<=alphas; k=k+0.001)
    {
      tolle=quantile(sims, LDs, k);
      digi=quantile(sims, LDs, (1-alphas)+k);
      
      i=0;
      max=0;
      while(i<sims)
      {
        
        if(LDs[i]<=digi)
        {
          if(LDs[i]>=max)
          {
            max = LDs[i];
          }
        }
        i++;
      }
      i=0;
      min=1;
      while(i<sims)
      {
        if(LDs[i]>=tolle)
        {
          if(LDs[i]<=min)
          {
            min=LDs[i];
          }
        }
        i++;
      }
      
      ups[counts]=max;
      lows[counts]=min;
      lengths[counts]=max-min;
      
      counts++;
    }
    
    counts=0;
    i=0;
    min=100;
    while(i<(alphas/0.001))
      
    {
      if(lengths[i]<min)
      {
        min=lengths[i];
        max=counts;
      }
      
      counts++;
      i++;
    }
    
    cilow[j]=lows[(int)max];
    ciup[j]=ups[(int)max];
  }
}
/*
 * qauntile definition for an array
 */
double quantile(int sims,  double LDs[],double alpha)
{
  double bound = alpha*(sims-1);
  if(floorf(bound)!=bound)
  {
    return (LDs[(int) floorf(bound)]*(1-fabs(bound - floorf(bound)))+ LDs[(int) floorf(bound)+1]*(fabs(bound - floorf(bound)))); 
  }
  else
  {
    return LDs[(int) floorf(bound)];
  }
  
  return 0;
}


/*
 * Sorting alg.
 */
void sort( double numbers[], int count)
{
  for(int i = 0; i < count - 1; i++)
  {
    double currentMin = numbers[i];
    int currentMinIndex = i;
    
    for(int j = i + 1; j < count; j++)
    {
      if(currentMin > numbers[j])
      {
        currentMin = numbers[j];
        currentMinIndex = j;
      }
    }
    
    if(currentMinIndex != i){
      numbers[currentMinIndex] = numbers[i];
      numbers[i] = currentMin;
      
    }
  }
}
/*
 * Construct Bootstrap table from genotype freuqencies
 * genofreq - genotype frequency table as intput for multinomial sampler
 * N - sample size
 * genoFreqBS - Genptype count table; sample table drawn from multinomial sampler
 */

void bootstrapGenoFreq( double genoFreq[4][4], int *N, double genoFreqBS[4][4], char **paradigm )
{
  double genoFreqVector[9], af[2]={0,0};
  int genoCountsBSVector[9], i=0, j=0, count=0;
  //double hapfreq[3][3]={0,0,0,0,0,0,0,0,0};
  //int Neff =2* *N;
  //double Dir[4]={1,1,1,1};
  /*  represent genotypic frequencies  as vector */
  for(i=0; i<3; i++)
  {	for(j=0;j<3;j++)
  {
    genoFreqVector[count] = genoFreq[i][j];
    genoFreqBS[i][j]=0;
    genoCountsBSVector[count] = 0;
    count++;
    //Rprintf("%f ", genoFreq[i][j]);
  }	
  //Rprintf("\n");
  }
  
  /* draw a multinomial sample with probabilities  */
  GetRNGstate();
  rmultinom(*N, genoFreqVector, 9, genoCountsBSVector);
  PutRNGstate();
  
  for(int l=0; l<9; l++)
  {
    if( genoCountsBSVector[l]==NAN)
    {bootstrapGenoFreq(genoFreq, N, genoFreqBS, paradigm);} 
  }
  
  //for(i=0;i<9;i++) Rprintf("%d ", genoCountsBSVector[i]);
  //Rprintf("\n");
  
  /* 	- the drawn sample is stored in  'genoCountsBSVector' 
   - convert this vector into a 4x4 matrix of
   genotypic frequencies with marginals 
   */
  count=0;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    {
      genoFreqBS[i][j] = (double)genoCountsBSVector[count]/ (double)(*N);
      count++;
    }
    
    genoFreqBS[3][0] = genoFreqBS[0][0]+genoFreqBS[1][0]+genoFreqBS[2][0];
  genoFreqBS[3][1] = genoFreqBS[0][1]+genoFreqBS[1][1]+genoFreqBS[2][1];
  genoFreqBS[3][2] = genoFreqBS[0][2]+genoFreqBS[1][2]+genoFreqBS[2][2];
  genoFreqBS[0][3] = genoFreqBS[0][0]+genoFreqBS[0][1]+genoFreqBS[0][2];
  genoFreqBS[1][3] = genoFreqBS[1][0]+genoFreqBS[1][1]+genoFreqBS[1][2];
  genoFreqBS[2][3] = genoFreqBS[2][0]+genoFreqBS[2][1]+genoFreqBS[2][2];
  genoFreqBS[3][3] = genoFreqBS[3][0]+genoFreqBS[3][1]+genoFreqBS[3][2];
  
  /* Rprintf(" %f", genoFreqBS[0][0]);
   Rprintf(" %f", genoFreqBS[0][1]);
   Rprintf(" %f", genoFreqBS[0][2]);
   Rprintf(" %f", genoFreqBS[1][0]);
   Rprintf(" %f", genoFreqBS[1][1]);
   Rprintf(" %f", genoFreqBS[1][2]);
   Rprintf(" %f", genoFreqBS[2][0]);
   Rprintf(" %f", genoFreqBS[2][1]);
   Rprintf(" %f", genoFreqBS[2][2]);*/
  
  return;
  //estimateFrequencies(genoFreqBS, hapfreq, &Neff,  Dir, mc, paradigm);
  //Rprintf(" %f", hapfreq[0][0]);
  /* check wether the allele frequencies are valid */
  alleleFreq(genoFreqBS, af);
  if( (af[0] > 0 && af[0] < 1) && ( af[1] > 0 && af[1] < 1 ))
  {
    return;
  }
  //if the allele frequencies are NOT valid, i.e. 0 or 1, the function is called recursively 
  else
  {
    bootstrapGenoFreq(genoFreq, N, genoFreqBS, paradigm);
  }
  
}

/*  CI construction methods: Bootstrap: quantile, se; Jackknife: pseudo values, leave one out
 * 
 * 
 * strategy - which CI method
  ci - computed or not
  mc - how many MC iterations
  alpha - significance level
  clow, cup - lower and upper bounds of CI
  nSim - number of Bootstrap samples
  LDdist - out for internal poperperties (e. g. LDs for Bootstrap samples)
    vars - output for standard error of LDs of CIs
      intervall - which CI method*/

void confidenceGenoInterval( double genoCounts[4][4],double genoFreq[4][4], int *N, char **paradigm,int *nSim, char **LD, int *LDnumb, double *tol, int *digits, double *LDdisti, double *HSweight, double *alpha, char **strategy, double *cilow, double *ciup, double tables[*nSim][9], double Dir[4], double vars[*LDnumb],int *intervall, int *mc)
{
  //set initail conditions
  int i=0, j=0;
  double genoFreqBS[4][4], haploCounts[3][3]={{0,0,0},{0,0,0},{0,0,0}},haploFreq[3][3]={{0,0,0},{0,0,0},{0,0,0}}, LDvaluei[*LDnumb], pooHat[3];
  double tolerance=*tol;
  int digtials=*digits;
  int sims=*nSim;
  /* initialize the matrix of bootstrapped genotype frequencies*/
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      genoFreqBS[i][j]=0;
  
  if((strcmp( strategy[0], "jackknife") == 0))
  {

    double observed_LDs[*LDnumb];
    // set sample size
    int new_hap_N=(2* *N);
    
    //Calculate genotype frequencies
    for(int k=0; k<4; k++)
    {
      for(int l=0; l<4;l++)
      {
        genoFreq[k][l]=genoCounts[k][l]/(*N);
      }
    }
    //estimate Haplotype frequencies with Cardano's method
    estimateHaploFreq(genoFreq, haploFreq,tol, digits, pooHat);
    // calculate Haplotype frequencies back to Haplotype counts
    for(int k=0; k<3; k++)
    {
      for(int l=0; l<3; l++)
      {
        haploCounts[k][l]= haploFreq[k][l]*(2 *(*N));
      }
    }
    //Estimate haplotype frequencies with estimator
    double exponent=0;
    estimateFrequencies(haploCounts,haploFreq, &new_hap_N, Dir, mc, paradigm);
    estimateLD(haploFreq, LD, LDnumb,observed_LDs, HSweight, &exponent);
    
    // Jackknife method: leave one observation out
    int jacktimes=0;
    int element=0;
    int tab=0;
    // for each genotype count table entry
    for(int a=0; a<3; a++)
    {
      for(int b=0; b<3; b++)
      {
        
        for(int j=0;j<*LDnumb;j++)
        {LDvaluei[j] = 0;}
        //Only calculate when genocounts > 0
        jacktimes=(int) round(genoCounts[a][b]);
        if((int) round(genoCounts[a][b]) > 0)
        {
          // substract one observation from population
          genoCounts[a][b]=(int) round(genoCounts[a][b])-1;
          //Resize table with new population size
          genoCounts[3][0] = genoCounts[0][0]+genoCounts[1][0]+genoCounts[2][0];
          genoCounts[3][1] = genoCounts[0][1]+genoCounts[1][1]+genoCounts[2][1];
          genoCounts[3][2] = genoCounts[0][2]+genoCounts[1][2]+genoCounts[2][2];
          genoCounts[0][3] = genoCounts[0][0]+genoCounts[0][1]+genoCounts[0][2];
          genoCounts[1][3] = genoCounts[1][0]+genoCounts[1][1]+genoCounts[1][2];
          genoCounts[2][3] = genoCounts[2][0]+genoCounts[2][1]+genoCounts[2][2];
          genoCounts[3][3] = genoCounts[3][0]+genoCounts[3][1]+genoCounts[3][2];
          
          tables[tab][0]=genoCounts[0][0];
          tables[tab][1]=genoCounts[0][1];
          tables[tab][2]=genoCounts[0][2];
          tables[tab][3]=genoCounts[1][0];
          tables[tab][4]=genoCounts[1][1];
          tables[tab][5]=genoCounts[1][2];
          tables[tab][6]=genoCounts[2][0];
          tables[tab][7]=genoCounts[2][1];
          tables[tab][8]=genoCounts[2][2];
          //Calculation of genotype frequencies
          int N_minus_1=*N-1;
          int hap_n_minus_one=N_minus_1*2;
          for(int k=0; k<4; k++)
          {
            for(int l=0; l<4; l++)
            {
              genoFreq[k][l]= genoCounts[k][l]/N_minus_1;
            }
          }
          //calculations of haplotype freuqncies
          estimateHaploFreq(genoFreq,haploFreq,tol,digits,pooHat);
          //Resize
          for(int k=0; k<3; k++)
          {
            for(int l=0; l<3; l++)
            {
              haploCounts[k][l]= haploFreq[k][l]*hap_n_minus_one;
            }
          }
          
          for(int j=0;j<*LDnumb;j++)
          {
            // estimate haplotype freuqenceis with estimator
            estimateFrequencies(haploCounts, haploFreq, &hap_n_minus_one, Dir, mc, paradigm);
            
            /*estimate LD*/
            double exponent=0;
            estimateLD(haploFreq, LD, LDnumb, LDvaluei, HSweight, &exponent);
          }
          // store calculated LD values in array
          for(int k=element; k<element+jacktimes; k++)
          {
            for(int f=0; f<*LDnumb; f++)
            {
              LDdisti[f*(*N)+k]=LDvaluei[f];
              //Rprintf(" %f", LDvaluei[f]);
            }
          }
          // Add the taken out observation
          genoCounts[a][b]=genoCounts[a][b]+1;
          
          genoCounts[3][0] = genoCounts[0][0]+genoCounts[1][0]+genoCounts[2][0];
          genoCounts[3][1] = genoCounts[0][1]+genoCounts[1][1]+genoCounts[2][1];
          genoCounts[3][2] = genoCounts[0][2]+genoCounts[1][2]+genoCounts[2][2];
          genoCounts[0][3] = genoCounts[0][0]+genoCounts[0][1]+genoCounts[0][2];
          genoCounts[1][3] = genoCounts[1][0]+genoCounts[1][1]+genoCounts[1][2];
          genoCounts[2][3] = genoCounts[2][0]+genoCounts[2][1]+genoCounts[2][2];
          genoCounts[3][3] = genoCounts[3][0]+genoCounts[3][1]+genoCounts[3][2];
          
          
          element+=jacktimes;
          tab++;
        }
      }
    }
    
    for(int e=0;e<*LDnumb;e++)
    {
      //Jackknife pseudo values method 
      if(*intervall==0)
      {
        // for one LD
        double LDs[((*N))];
        for(int h=0; h<(*N); h++)
        {
          LDs[h]=0;
        }
        int intList=0;
        // pseudo values approach
        for(int g=e*(*N); g<e*(*N)+(*N); g++)
        {
          //LDs[intList]=atanh(LDdisti[g]);// - atanh(observed_LDs[e]));
          if((strcmp( paradigm[0], "freq") == 0))
          {
            LDs[intList]=((*N)*observed_LDs[e]) - ((*N-1)*LDdisti[g]);
          }
          else
          {
            LDs[intList] = ((*N)*atanh(observed_LDs[e])) - ((*N-1)*atanh(LDdisti[g]));
          }
          //Rprintf ("%f ", LDdisti[g]);
          //Rprintf(" %f", observed_LDs[e]);
          
          intList++;
        }
        // calculate mean
        double sum_ps=0;
        for(int a=0; a<(*N); a++)
        {
          sum_ps+=LDs[a];
          
          
        }
        sum_ps= sum_ps/(*N);
        
        //calculate standard error
        double var_ps=0;
        
        for(int a=0; a<(*N); a++)
        {
          var_ps+=(pow((LDs[a]-sum_ps),2));
        }

        var_ps=sqrt(var_ps*1/(((*N)-1)*((*N))));
        
        
        vars[e]=tanh(var_ps);

        // Student t quantile
        double stud_t=0;
        double alphafull=1-(*alpha/2);
        double sim_minus_1=(*N)-1;
        int logger=0;
        int tail=1;
        //Rf_qt(alphafull,sim_minus_1);
        stud_t= Rf_qt(alphafull, sim_minus_1,tail,logger);

        //Calculate CI
        long double Dupper=0;
        long double Dlower=0;
        if((strcmp( paradigm[0], "freq") == 0))
        {
          Dupper=observed_LDs[e]+stud_t*var_ps;
          Dlower=observed_LDs[e]-stud_t*var_ps;        
        }
        else
        {
          Dupper=tanh(atanh(observed_LDs[e])+stud_t*var_ps);
          Dlower=tanh(atanh(observed_LDs[e])-stud_t*var_ps);        
        }
        
        if(Dlower < -1){Dlower=-1;}
        if(Dupper > 1){Dupper=1;}
        
        cilow[e]=Dlower;
        ciup[e]=Dupper;
      }
      // Jackknife: leave one out method
      else
      {
        double LDs[((*N))];
        for(int h=0; h<(*N); h++)
        {
          LDs[h]=0;
        }
        int intList=0;
        for(int g=e*(*N); g<e*(*N)+(*N); g++)
        {
          // Fischer Transformation
          if((strcmp( paradigm[0], "freq") == 0))
          {
            LDs[intList]=LDdisti[g];
          }
          else
          {
            LDs[intList]=atanh(LDdisti[g]);
          }
          intList++;
        }
        // Sum
        double sum_ps=0;
        for(int a=0; a<(*N); a++)
        {
          sum_ps+=LDs[a];
        }
        sum_ps= sum_ps/(*N);
        //Rprintf(" %f", sum_ps);
        
        double var_ps=0;
        // standard error
        for(int a=0; a<(*N); a++)
        {
          var_ps+=(pow((LDs[a]-sum_ps),2));
        }
        //Rprintf(" %f", var_ps);
        
        double new_ns=(*N);
        double new_ns_min_one=new_ns-1;
        double ratio=(new_ns_min_one)/new_ns;
        var_ps=sqrt(var_ps*ratio);

        vars[e]=tanh(var_ps);

        // calculate quantile        
        double stud_t=0;
        double alphafull=1-(*alpha/2);
        double sim_minus_1=(*N)-1;
        int logger=0;
        int tail=1;
        //Rf_qt(alphafull,sim_minus_1);
        stud_t= Rf_qt(alphafull, sim_minus_1,tail,logger);
        //Calculate CIs
        
        long double Dupper=0;
        long double Dlower=0;
        if((strcmp( paradigm[0], "freq") == 0))
        {
          Dupper=observed_LDs[e]+stud_t*var_ps;
          Dlower=observed_LDs[e]-stud_t*var_ps;        
        }
        else
        {
          Dupper=tanh(atanh(observed_LDs[e])+stud_t*var_ps);
          Dlower=tanh(atanh(observed_LDs[e])-stud_t*var_ps);        
        }
        
        if(Dlower < -1){Dlower=-1;}
        if(Dupper > 1){Dupper=1;}
        
        cilow[e]=Dlower;//D_prime - 1.644854*sqrt(Var_Zapata);
        ciup[e]=Dupper;
      }
    }
    
  }
  // Bootstrap: quantile, standard error
  int hap_N=2* *N;
  if((strcmp( strategy[0], "bootstrap") == 0))
  {
    double observed_LDs[*LDnumb];
    
    int new_hap_N=(2* *N);
    //Calc genotype freqs
    for(int k=0; k<4; k++)
    {
      for(int l=0; l<4;l++)
      {
        genoFreq[k][l]=genoCounts[k][l]/(*N);
      }
    }
    // est haplo freqs Cardano's method
    estimateHaploFreq(genoFreq, haploFreq,tol, digits, pooHat);
    // back to haplotype counts
    for(int k=0; k<3; k++)
    {
      for(int l=0; l<3; l++)
      {
        haploCounts[k][l]= haploFreq[k][l]*(2 *(*N));
      }
    }
    // est haplo freqs
    estimateFrequencies(haploCounts,haploFreq, &new_hap_N, Dir, mc, paradigm);
    double exponent=0;
    //Est LD
    estimateLD(haploFreq, LD, LDnumb,observed_LDs, HSweight, &exponent);
    
    
    /* loop over the number of simulations */
    for(i=0; i<*nSim; i++ )
    {
      for(int p=0;p<3;p++)
      { 
        for(int q=0;q<3;q++)
        {
          genoFreqBS[p][q]=0;
          haploCounts[p][q]=0;
          haploFreq[p][q]=0;
        }
      }
      
      for(j=0;j<*LDnumb;j++)
      {LDvaluei[j] = 0;}
      
      /* generate a bootstraped 3x3 table of genotypic frequencies */
      bootstrapGenoFreq(genoFreq, N, genoFreqBS, paradigm);
      
      
      /*tables[i][0]=genoFreqBS[0][0]*(*N);
      tables[i][1]=genoFreqBS[0][1]*(*N);
      tables[i][2]=genoFreqBS[0][2]*(*N);
      tables[i][3]=genoFreqBS[1][0]*(*N);
      tables[i][4]=genoFreqBS[1][1]*(*N);
      tables[i][5]=genoFreqBS[1][2]*(*N);
      tables[i][6]=genoFreqBS[2][0]*(*N);
      tables[i][7]=genoFreqBS[2][1]*(*N);
      tables[i][8]=genoFreqBS[2][2]*(*N);*/
      /* estimate haplotype frequencies */
      *tol=tolerance;
      *digits=digtials;
      /*pooHat[0]=0;
      pooHat[1]=0;
      pooHat[2]=0;*/
      estimateHaploFreq(genoFreqBS, haploCounts, tol, digits, pooHat);
      for(int a=0;a<3;a++)
      {
        for(int b=0;b<3;b++)
        {
          haploCounts[a][b]=haploCounts[a][b]*(2* *N);
        }
      }
      
      //Rprintf(" %f", haploCounts[2][2]);
      estimateFrequencies(haploCounts, haploFreq, &hap_N, Dir, mc, paradigm);
      
      /* estimate LD */
      double exponent=0;
      estimateLD(haploFreq, LD, LDnumb, LDvaluei, HSweight, &exponent);
      //Rprintf ("%f ", LDvaluei[0]);
      /* store the result  in a vector */
      for(j=0;j<*LDnumb;j++)
      {
        LDdisti[ j* (*nSim) + i] = LDvaluei[j];
      }
    }
    //Rprintf ("%d ", hap_N);
    for(int e=0;e<*LDnumb;e++)
    {
      // construct quantile
      if(*intervall==0)
      {
        double LDs[*nSim];
        for(int h=0; h<*nSim; h++)
        {
          LDs[h]=0;
        }
        
        int intList=0;
        for(int g=e*(*nSim); g<e*(*nSim)+(*nSim); g++)
        {
          LDs[intList]=LDdisti[g];
          
          //Rprintf ("%d ", g);
          intList++;
        }
        sort(LDs, sims);
        double alphas=*alpha;
        double alpha_u = alphas/2;
        double tolle=quantile(sims, LDs, alpha_u);
        //Rprintf ("%f ", tolle);
        //digits=quantile(LDdist, *nSim);
        double alpha_l = 1-(alphas/2);
        double digi=quantile(sims, LDs, alpha_l);
        //Rprintf ("%f ", digi);
        // char strat=strategy;
        //LDdist=0;
        cilow[e]=tolle;
        ciup[e]=digi;
      }
      // standard error method
      else
      {
        double LDs[(*nSim)];
        //double mean_ld=0;
        for(int h=0; h<*nSim; h++)
        {
          LDs[h]=0;
        }
        // Fischer Tranformation
        int intList=0;
        for(int g=e*(*nSim); g<e*(*nSim)+(*nSim); g++)
        {

          if((strcmp( paradigm[0], "freq") == 0))
          {
            LDs[intList]=LDdisti[g];
          }
          else
          {
            LDs[intList]=atanh(LDdisti[g]);
          }
          intList++;
        }
        //Mean
        double sum_ps=0;
        for(int a=0; a<*nSim; a++)
        {
          sum_ps+= LDs[a];
        }
        sum_ps = sum_ps/(*nSim);
        //Rprintf ("%f ", sum_ps);
        
        double var_ps=0;
        // Standard error
        for(int a=0; a<*nSim; a++)
        {
          var_ps+=(pow((LDs[a]-sum_ps),2))/(*nSim-1);
        }
        var_ps=sqrt(var_ps);
        
        vars[e]=var_ps;
        //var_ps=sqrt(var_ps*(new_hap_N-1)/(new_hap_N));
        
        //quantile
        double stud_t=0;
        double alphafull=1-(*alpha/2);
        double sim_minus_1=*N-1;
        int logger=0;
        int tail=1;
        //Rf_qt(alphafull,sim_minus_1);
        stud_t= Rf_qt(alphafull, sim_minus_1,tail,logger);  
        //Rprintf(" %f", stud_t);
        //double Dupper=  observed_LDs[e] + stud_t * sqrt((var_ps));
        //double Dlower = observed_LDs[e] - stud_t * sqrt((var_ps));
        
        // CI construction
        long double Dupper=0;
        long double Dlower=0;
        if((strcmp( paradigm[0], "freq") == 0))
        {
          Dupper=observed_LDs[e]+stud_t*var_ps;
          Dlower=observed_LDs[e]-stud_t*var_ps;        
        }
        else
        {
          Dupper=tanh(atanh(observed_LDs[e])+stud_t*var_ps);
          Dlower=tanh(atanh(observed_LDs[e])-stud_t*var_ps);        
        }
        
        //double Dupper = tanh(sum_ps + 1.984217 *(var_ps));
        //double Dlower = tanh(sum_ps - 1.984217 *(var_ps));
        
        if(Dlower < -1){Dlower=-1;}
        if(Dupper > 1){Dupper=1;}
        
        cilow[e]=Dlower;//D_prime - 1.644854*sqrt(Var_Zapata);
        ciup[e]=Dupper;
      }
    }
  }

  return;
}


/*Karsten Krug 2008-03-28*/

# include<R.h>
# include<Rmath.h>
# include <Rdefines.h>
# include <Rinternals.h>
# include "my_solve_poly.h"
# include "pwLD.h"
# include "LDmeasures.h"
# include "distribution_statistics.h"

/*#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include "LDmeasures.h"*/

#ifndef oneDivLn2
#define oneDivLn2           1.442695040888963  /* 1/ln(2), scale factor for Mutual Information */
#endif


/* ###############################################################
 * implementation of association maesure HS
 */

double HS( double tab[3][3],  double *HSweight, double *exponent)
{
  double x = log(sqrt((tab[0][0]*tab[1][1])/(tab[0][1]*tab[1][0])));
  double y = log(sqrt((tab[0][0]*tab[0][1])/(tab[1][0]*tab[1][1])));
  double z = log(sqrt((tab[0][0]*tab[1][0])/(tab[0][1]*tab[1][1])));
  long double yule = tanhl(x/2);
  double H_diag = 1+log2(1+exp(x))-(x/(log(2)*(exp(-x)+1)));
  double H = log2(exp(x+y+z)+exp(x)+exp(y)+exp(z))-(((x+y+z)*(exp(x+y+z))+x*exp(x)+y*exp(y)+z*exp(z))/(log(2)*(exp(x+y+z)+exp(x)+exp(y)+exp(z))));
  double exponenten=(exp((*HSweight) *(H_diag-H)));
  double HS = sign(yule)*(pow(fabs(yule),exponenten));
  *exponent=exponenten;
  return(HS);
  
}

/* ###############################################################
 * implementation of Yule's Y (naive plug-in estimator)
 
 */
double Y(  double tab[3][3])
{
  double tmp1=0, tmp2=0;
  tmp1 = tab[0][0] * tab[1][1] ;
  tmp2 = tab[1][0] * tab[0][1];
  
  return(sqrt(tmp1)-sqrt(tmp2))/(sqrt(tmp1)+sqrt(tmp2));
}



/* ###############################################################
 * implementation of Yule's Q (naive plug-in estimator)
 
 */
double Q(  double tab[3][3])
{
  double Qval, D;
  
  /* disequilibrium coefficient */
  D = tab[0][0] - tab[0][2]*tab[2][0];
  
  /* Q */
  Qval = D/(tab[0][0]*tab[1][1] + tab[0][1]*tab[1][0] );
  
  /* make it robust against numerical issues */
  if(Qval > 0) Qval = fmin2(1.0, Qval);
  if(Qval < 0) Qval = fmax2(-1.0, Qval);
  
  return Qval;
}


/* ###############################################################
 * implementation of Lewontin's D'
 *tab is a pointer to a vector of size 4: p00, p01, p10, p11  
 */
double Dprime(  double tab[3][3])
{
  double Dmax=0, D=0;
  
  /* disequilibrium coefficient */
  D = tab[0][0] - tab[0][2]*tab[2][0];
  
  /* Dmax */
  if(D < 0) 
  { Dmax = fmin2( tab[0][2]*tab[2][0], tab[2][1]*tab[1][2]);}
  else if (D > 0) 
  {Dmax = fmin2( tab[0][2]*tab[2][1], tab[1][2]*tab[2][0]);}
  else {Dmax = NA_REAL ;}
  
  if(Dmax==0){return 0;}
  /*D'*/
  if(D > 0 ) {return fmin2(1.0, (D/Dmax) );}
  else if(D < 0) {return fmax2(-1.0, (D/Dmax));}
  else {return 0;}
}

/*###############################################################
 * implementation of r^2
 */
double R( double tab[3][3])
{
  double D=0, r=0;
  
  /* disequilibrium coefficient */
  D = tab[0][0] - tab[0][2]*tab[2][0];
  
  r= D/sqrt(tab[0][2]*tab[1][2]*tab[2][0]*tab[2][1]);
  
  return fmin2(r, 1.0);
}

/*##############################################################
 * 	odds ratio
 */
double oddsratio( double tab[3][3])
{
  double tmp1=0, tmp2=0;
  
  tmp1 = tab[0][0] * tab[1][1] ;
  tmp2 = tab[1][0] * tab[0][1];
  
  return  (tmp1 / tmp2);
}
/*##############################################################
 *	disequilibrium coefficient
 */
double D( double tab[3][3])
{
  return tab[0][0] - tab[0][2]*tab[2][0];
}

/*#############################################################
 *	Mutual Information
 */
double MI( double tab[3][3])
{
  int i=0, j=0;
  double MI=0, tmp=0;
  
  for(i=0; i< 2; i++)
  {	for(j=0; j<2; j++)
  {
    //if(tab[i][j] < 0.0000001) tmp = 0;
    //else tmp = tab[i][j] * log1p
    tmp =tab[i][j] * log1p( ( tab[i][j]/( (tab[i][0]+tab[i][1])*(tab[0][j]+tab[1][j]) ) )  -  1 );
    
    MI += tmp;
  }
  }
  
  return fmin2( 1.0 , MI*oneDivLn2);
}

/*############################################################
 *	Pearson's chi-square goodness of fit test between observed and expected
 haplotype frequencies assuming complete linkage equilibrium
 */
double ChiSquareHaplotypeFrequencies( double O[3][3])
{
  double  N=O[2][2], af1=O[0][2]/N, af2=O[2][0]/N, PooE = (af1*af2), chi2=0; 
  double E[3][3] = { {PooE*N, ((af1-PooE)*N), af1*N}, {(af2-PooE)*N, ((1 - af1 - af2 + PooE)*N), (1- af1)*N}, { af2*N, (1-af2)*N, 1*N} };
  
  int i=0, j=0;
  
  for(i=0; i< 3; i++)
    for(j=0; j<3; j++)
      chi2 += ( (O[i][j] - E[i][j])*(O[i][j] - E[i][j]) ) / E[i][j];
  
  return chi2;
}



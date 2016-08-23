/* Karsten Krug, 2008-04-06 
 * 
 * implementation of some functions from statistics
 * */

// #include <gsl/gsl_sf_gamma.h>
 #include <R.h>
 #include <Rmath.h>
 #include "distribution_statistics.h"

/* probability mass function, multinomial distribution */ 
double dmultinom(int size, unsigned int *n, double *p)
{ 
	double result, B=1, A=1;
	int i;
	unsigned int N=0;
	
	
	for(i=0;i<size;i++)
	{
		N += n[i];
		B *= pow(p[i], n[i]);
		A *= factorial(n[i]);
	}

	result= (factorial(N)/A)*B;

	return result;		
	
}


/*	factorial function */
unsigned int factorial(unsigned int a)
{
	
	unsigned int i, result=1;
	 
	for(i = a; i > 1; i--)
			result*=i;	
	
	return result;
}	


/*#################################################
	random number generation from a Dirichlet distribution

##################################################*/

/*	implementation by myself, translated from the C implementation above */
/* Dirichlet distribution, random generator*/
void rdirichlet( double* alpha, int K, double *theta )
{
	double norm=0.0;
	int i;
	
	for(i = 0; i < K; i++) 
	{
		GetRNGstate();
    		theta[i] = rgamma( alpha[i], 1.0);
		PutRNGstate();
	}

  	for(i = 0; i < K; i++) 
    		norm += theta[i];
  
  	for(i = 0; i < K; i++) 
    		theta[i] /= norm;
	
}


/* probability density function, Dirichlet distribution */
double ddirichlet(double* x, double* alpha, int K)
{
	double f=1, Beta=1, sumAlpha=0;
	int i;

	for(i = 0; i < K; i++)
	{	
		f*=pow(x[i], alpha[i]-1 );
		sumAlpha += alpha[i];
	
		Beta *= gammafn(alpha[i]);	
	}

	Beta = Beta/gammafn(sumAlpha);

	return f/Beta;

}

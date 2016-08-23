/* Karsten Krug, 2008-06-12 

	C code for package  pwLD

*/

# include<R.h>
# include<Rmath.h>
# include <Rdefines.h>
 
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
void estimateLD(double tab[3][3], char **what, int *numb, double *LD)
{
	int i=0;
	
	for(i=0; i<*numb; i++)
	{		
		if( strcmp(what[i], "Q") == 0  )
				LD[i] = Q( tab);
		if( strcmp(what[i], "Dprime") == 0 )
				LD[i] = Dprime(tab);
		if( strcmp(what[i], "r2") == 0 )
				LD[i] = Rsquare(tab);
		if( strcmp(what[i], "OR") == 0 )
				LD[i] = oddsratio(tab);
		if( strcmp(what[i], "D") == 0 )
				LD[i] = D(tab);
		if( strcmp(what[i], "MI") == 0 )
				LD[i] = MI(tab);
		if( strcmp(what[i], "chi2") == 0 )
				LD[i] = ChiSquareHaplotypeFrequencies(tab);


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

	/* estimate frequencies */
	if( strcmp( *paradigm, "freq") == 0)
	{
		for(i=0;i<4;i++)
			for(j=0;j<4;j++)
				 genoFreq[i][j] = genoCounts[i][j] / (*Neff); 
	}
	
	if( strcmp( *paradigm, "bayes") == 0 )	
	{	
		int nTmp=0;
		for(i=0;i<3;i++)
		{	for(j=0;j<3;j++)
			{
				genoFreq[i][j] = (genoCounts[i][j] + Dir[nTmp])/(*Neff + DirSum);
				nTmp++;
			}				
		}

		/* marginals */
		genoFreq[0][3] = genoFreq[0][0] + genoFreq[0][1] + genoFreq[0][2];
		genoFreq[1][3] = genoFreq[1][0] + genoFreq[1][1] + genoFreq[1][2];
		genoFreq[2][3] = genoFreq[2][0] + genoFreq[2][1] + genoFreq[2][2];
		genoFreq[3][0] = genoFreq[0][0] + genoFreq[1][0] + genoFreq[2][0];
		genoFreq[3][1] = genoFreq[0][1] + genoFreq[1][1] + genoFreq[2][1];
		genoFreq[3][2] = genoFreq[0][2] + genoFreq[1][2] + genoFreq[2][2];

		genoFreq[3][3] = genoFreq[0][3] + genoFreq[1][3] + genoFreq[2][3];

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

*/
void LDall(  char **data, int *nrow, int *ncol, char **LD, int *LDnumb,  char
**code, char **paradigm, double *Dir, double *MAF, double *tol, int *digits, double *LDmatPtr)
{
	int i=0, j=0, k=0, m=0, n=0,  N=0, countEntries=0, Nentries=(*nrow)*(*nrow-1)/2;
	char *genoSnp1[*ncol], *genoSnp2[*ncol];
	double genoFreq[4][4], genoCounts[4][4], haploFreq[3][3], af[2]={0,0}, *LDptr=0, LDvalue[*LDnumb], pooHat[3]; 

	LDptr = LDvalue;

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
					if((n < 3) && (m < 3))
						haploFreq[n][m] = NA_REAL;
				}
	
			/* estimate genotypic frequencies*/
			estimateGenoFreq( genoSnp1, genoSnp2, code, ncol, paradigm,  Dir,  genoFreq, genoCounts, &N );



			//estimateGenoFreq( data[i], data[j], code, ncol, paradigm,  Dir,  genoFreq, Nptr );
			
			/* estimate allele frequencies from genotypic
			frequencies*/
			alleleFreq(genoFreq, af);

			/* check wether the allele frequencies are above the MAF
			threshold */
			//if( (af[0] > *MAF) && (af[0] < (1-*MAF)) && (af[1] > *MAF) && (af[1] < (1-*MAF)) )	
			if( (af[0] > 0.0) && (af[0] < 1.0) && (af[1] > 0.0) && (af[1] < 1.0) )	
			{
				/* estimate haplotype frequencies */
				estimateHaploFreq( genoFreq, haploFreq, tol, digits, pooHat );

				/*estimate LD*/
				estimateLD(haploFreq, LD, LDnumb, LDptr);	

				/* store the result */
				for(k=0; k < *LDnumb; k++)
					//LDmatPtr[k][countEntries ] = LDvalue[k];	
					LDmatPtr[ k*Nentries+ countEntries ] = LDvalue[k];

			}
			/* if any allele frequency does not suffice the MAF
			criteria, the respective LD value is set to 'NA' */
			else 
			{	
				for(k=0; k < *LDnumb; k++) 
						LDmatPtr[ k*Nentries + countEntries ] = NA_REAL;
	
						//LDmatPtr[k][countEntries] =NA_REAL;
			}
		
			/* increment the number of entries of the upper
			triangular matrix */
			countEntries++;
		}	
				
	}
}

/*#########################################

	generates a bootstrap sample from genotypic frequencies by
	multinomial sampling from 'genoFreq' with sample size 'N'

	genoFreq		- 4x4 matrix of genotypic frequencies with marginals
	N			- sample size
	genoFreqBS	- 4x4 matrix to store the bootstrapped genotype frequencies
*/
void bootstrapGenoFreq( double genoFreq[4][4], int *N, double genoFreqBS[4][4] )
{
	double genoFreqVector[9], af[2]={0,0};
	int genoCountsBSVector[9], i=0, j=0, count=0;
	
	/*  represent genotypic frequencies  as vector */
	for(i=0; i<3; i++)
	{	for(j=0;j<3;j++)
		{
			genoFreqVector[count] = genoFreq[i][j];
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
			genoFreqBS[i][j] = (double)genoCountsBSVector[count]/(double)(*N);
			count++;
		}

	/* initialize the marginals of the table with zero */	
	for(i=0; i<3; i++) 
	{ 	
		genoFreqBS[i][3]=0;
		genoFreqBS[3][i]=0;
	}
	/* estimate marginal distributions */
	for(i=0;i<3;i++)
	{
		genoFreqBS[0][3] +=genoFreqBS[0][i];
		genoFreqBS[1][3] +=genoFreqBS[1][i];
		genoFreqBS[2][3] +=genoFreqBS[2][i];

		genoFreqBS[3][0] +=genoFreqBS[i][0];
		genoFreqBS[3][1] +=genoFreqBS[i][1];
		genoFreqBS[3][2] +=genoFreqBS[i][2];
	}

	genoFreqBS[3][3] = genoFreqBS[3][2] + genoFreqBS[3][1] + genoFreqBS[3][0];

	/* check wether the allele frequencies are valid */
	alleleFreq(genoFreqBS, af);
	if( (af[0] > 0 && af[0] < 1) && ( af[1] > 0 && af[1] < 1 )) return;
	/* if the allele frequencies are NOT valid, i.e. 0 or 1, the function is called recursively */
	else
		bootstrapGenoFreq(genoFreq, N, genoFreqBS);

}


/*#####################################################

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
void confidenceInterval( double genoFreq[4][4], int *N, int *nSim, char **LD, int *LDnumb, double *tol, int *digits,double *LDdist )
{
	int i=0, j=0;
	double genoFreqBS[4][4], haploFreq[3][3]={{0,0,0},{0,0,0}, {0,0,0}}, LDvalue[*LDnumb], pooHat[3];

	/* initialize the matrix of bootstrapped genotype frequencies*/
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			genoFreqBS[i][j]=0;


	/* loop over the number of simulations */
	for(i=0; i<*nSim; i++ )
	{
		for(j=0;j<*LDnumb;j++) LDvalue[j] = 0;
	
		/* generate a bootstraped 3x3 table of genotypic frequencies */
		bootstrapGenoFreq(genoFreq, N, genoFreqBS);

		/* estimate haplotype frequencies */
		estimateHaploFreq(genoFreqBS, haploFreq, tol, digits, pooHat);

		/* estimate LD */
		estimateLD(haploFreq, LD, LDnumb, LDvalue);

		/* store the result  in a vector */
		for(j=0;j<*LDnumb;j++)
			LDdist[ j* (*nSim) + i] = LDvalue[j];

	}

}

/*############################################

	generates a new set of genotypic frequencies assuming 
	Dirichlet prior

*/
void sampleDirichletGenoFreq( double genoCounts[4][4], double *Dir, double genoFreqRS[4][4] )
{
	double shapeDirichlet[9]={0,0,0,0,0,0,0,0,0}, dirichletGenoFreq[9]={0,0,0,0,0,0,0,0,0}, af[2]={0,0};
	int i=0, j=0, count=0; // N=genoCounts[3][3];
	

	/*  determine the shape parameters of the Dirichlet distribution to sample from, i.e. n_ij + Dir_ij */
	for(i=0; i<3; i++)
	{	for(j=0;j<3;j++)
		{
			shapeDirichlet[count] = genoCounts[i][j] + Dir[count];
			count++;
			//Rprintf("%f ", genoCounts[i][j]);
		}	
		//Rprintf("\n");
	}

	/* sample from Dirichlet( n_ij + Dir_ij )  */
	rdirichlet( shapeDirichlet, 9, dirichletGenoFreq);

	
	//for(i=0;i<9;i++) Rprintf("%f ", dirichletGenoFreq[i]);
	//Rprintf("\n");

	/* 	- the drawn sample is stored in  'dirichletGenoFreq' 
		- convert this vector into a 4x4 matrix of
			genotypic frequencies with marginals 
	*/
 	count=0;
 	for(i=0;i<3;i++)
 		for(j=0;j<3;j++)
 		{
 			genoFreqRS[i][j] = dirichletGenoFreq[count];
 			count++;
 		}

	/* initialize the marginals of the table with zero */	
	for(i=0; i<3; i++) 
	{ 	
		genoFreqRS[i][3]=0;
		genoFreqRS[3][i]=0;
	}
	/* estimate marginal distributions */
	for(i=0;i<3;i++)
	{
		genoFreqRS[0][3] +=genoFreqRS[0][i];
		genoFreqRS[1][3] +=genoFreqRS[1][i];
		genoFreqRS[2][3] +=genoFreqRS[2][i];

		genoFreqRS[3][0] +=genoFreqRS[i][0];
		genoFreqRS[3][1] +=genoFreqRS[i][1];
		genoFreqRS[3][2] +=genoFreqRS[i][2];
	}

	genoFreqRS[3][3] = genoFreqRS[3][2] + genoFreqRS[3][1] + genoFreqRS[3][0];

	/* check wether the allele frequencies are valid */
	alleleFreq(genoFreqRS, af);
	if( (af[0] > 0 && af[0] < 1) && ( af[1] > 0 && af[1] < 1 )) return;
	/* if the allele frequencies are NOT valid, i.e. 0 or 1, the function is called recursively */
 	else
 		sampleDirichletGenoFreq(genoCounts, Dir, genoFreqRS);

}
/*###########################################
	
		estimation of Bayesian Credible Intervals 

*/

void credibleInterval(double genoCounts[4][4], int *nSim, double Dir[9], char **LD, int *LDnumb, double *tol, int *digits,double *LDdist)
{
	int i=0, j=0;
	double genoFreqRS[4][4], haploFreq[3][3]={{0,0,0},{0,0,0}, {0,0,0}}, LDvalue[*LDnumb], pooHat[3];

	/* initialize the matrix of simulated genotype frequencies*/
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			genoFreqRS[i][j]=0;

	for(i=0; i< *nSim;i++)
	{
		/* simulate genotypic frequencies assuming a Dirichlet prior with shape parameter n_ij + Dir_ij */
		sampleDirichletGenoFreq(genoCounts, Dir, genoFreqRS);
	
		
		/* estimate haplotype frequencies */
		estimateHaploFreq(genoFreqRS, haploFreq, tol, digits, pooHat);

		/* estimate LD */
		estimateLD(haploFreq, LD, LDnumb, LDvalue);

		/* store the result  in a vector */
		for(j=0;j<*LDnumb;j++)
			LDdist[ j* (*nSim) + i] = LDvalue[j];
	}

}


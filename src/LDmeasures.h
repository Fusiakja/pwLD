/* Karsten Krug, 2008-03-28
 *   */

#ifndef LDMEASURES_H_
#define LDMEASURES_H_


double HS( double tab[3][3],  double *HSweight, double *exponent);
double Y(  double tab[3][3]);
/* ###############################################################
 * implementation of Yule's Q
 *tab is a pointer to a vector of size 4: p00, p01, p10, p11  
 */
double Q(  double tab[3][3]);

/* ###############################################################
 * implementation of Lewontin's D'
 *tab is a pointer to a vector of size 4: p00, p01, p10, p11  
 */
double Dprime(  double tab[3][3]);

/*###############################################################
 * implementation of r
 */
double R( double tab[3][3]);

/*##############################################################
 * 	odds ratio
 */
double oddsratio( double tab[3][3]);

/*##############################################################
 *	disequilibrium coefficient
 */
double D( double tab[3][3]);


/*#############################################################
 *	Mutual Information
 */
double MI( double tab[3][3]);

/*############################################################
 *	Pearson's chi-square goodness of fit test between observed and expected
 haplotype frequencies assuming complete linkage equilibrium
 */
double ChiSquareHaplotypeFrequencies( double O[3][3]);


#endif /*LDMEASURES_H_*/
